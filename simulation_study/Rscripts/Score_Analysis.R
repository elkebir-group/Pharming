library(tidyverse)
library(yaml)
library(glue)
figpath <- "figures"
xvert <- theme(axis.text.x = element_text(
    angle = 90,  # Rotate labels to be vertical
    hjust = 1,   # Adjust horizontal justification
    vjust = 0.5  # Adjust vertical justification
  ))
compute_acc <- function(df){
  df <- df %>%
    mutate(anc_mut_pairs = anc_mut_TP + anc_mut_FN,
           inc_mut_pairs = inc_mut_TP + anc_mut_FN,
           mut_pairs = mut_TP + mut_FN,
           mut_total = anc_mut_pairs + inc_mut_pairs + mut_pairs,
           anc_mut_wt = anc_mut_pairs/mut_total,
           inc_mut_wt = inc_mut_pairs/mut_total,
           mut_wt = mut_pairs/mut_total,
           mut_sum_score = mut_wt*mut_recall + anc_mut_wt*anc_mut_recall + inc_mut_wt*inc_mut_recall,
           anc_cell_pairs = anc_cell_TP + anc_cell_FN,
           inc_cell_pairs = inc_cell_TP + inc_cell_FN,
           cell_pairs = cell_TP + cell_FN,
           cell_total = anc_cell_pairs + inc_cell_pairs + cell_pairs,
           anc_cell_wt = anc_cell_pairs/cell_total,
           inc_cell_wt = inc_cell_pairs/cell_total,
           cell_wt = cell_pairs/cell_total,
           cell_sum_score = anc_cell_wt*anc_cell_recall + inc_cell_wt*inc_cell_recall + cell_wt*cell_recall
    )
  return(df)
}
#assuming working directory is /scratch/data/leah/Pharming/simulation_study
# metric_names <-  c("Ancestral Pair\nRecall", "Incomparable Pair\nRecall",
#                 "Clustered Pair\nRecall", "Cell\nARI", "CNA Tree\nAccuracy")

metric_labels <- c(" SNV ancestral pair\n recall (APR)",
                  "SNV incomparable pair\nrecall (IPR)" ,
                  "SNV clustered pair\nrecall (CPR)",
                  "SNV accuracy",
                  "cell ancestral pair\n recall (APR)",
                  "cell incomparable pair\nrecall (IPR)" ,
                  "cell clustered pair\nrecall (CPR)",
                  "cell accuracy",
                  "CNA tree accuracy",
                  "SNV tree accuracy",
                  "CNA gentoype similarity",
                  "SNV genotype similarity",
                  "SNV presence similarity",
                  "SNV genotype allele-specific similarity"  )

metric_names <- c("anc_mut_recall", "inc_mut_recall", "mut_recall", "mut_sum_score",
  "anc_cell_recall", "inc_cell_recall", "cell_recall", "cell_sum_score",
  "cna_tree_accuracy", "snv_tree_accuracy",
  "cna_genotype_similarity", "genotype_similarity",
  "genotype_pres_similarity","genotype_allele_spec_similarity"
)
col_to_name_map <- data.frame(name=metric_names,
                              full_name = metric_labels)
header <- read.csv("cpp/metric_header.csv")


read_csv_wrapper <- function(folder, prefix, suffix){
  fname <- file.path(prefix, folder, suffix)
  if(file.exists(fname)){
    df<- read.csv(fname, header=FALSE)
    colnames(df) <- colnames(header)
     df <- df %>%
    mutate(folder = folder) 
  }else{
    print(glue("{fname} does not exist!"))
    df <- data.frame()
  }

  return(df)
}

bpath <- "/scratch/leah/Pharming"
config = yaml.load_file("simulate.yml")
pharming_config <- yaml.load_file( "pharming.yml")
seeds <- 10:14
runs <- expand.grid(cells=config$cells, 
                    snvs=config$snvs,
                     cov= config$cov,
                    mclust = config$mclust,
                    err = config$cerror, 
                    s = seeds, 
                    segs=config$nsegs,
                    dirch = config$dirch) %>%
                    mutate_all(factor)


phert_runs <- runs %>% 
    mutate(folder = glue("s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))
phert <- bind_rows(lapply(phert_runs$folder,  read_csv_wrapper, "phertilizer/sims", "metrics.csv")) %>%
  group_by(folder) %>%
  mutate(clust_method = "CN", method="Phertilizer", dcf_values="UMAP") %>% inner_join(phert_runs)

saveRDS(compute_acc(phert), "data/phertilizer.rds")


pharm_runs <-  expand.grid(
    isegs = pharming_config$ninit_segs,
    tm = pharming_config$ninit_tm,
    order = pharming_config$order,
    topn = pharming_config$topn,
    lamb = pharming_config$lamb,
    clust_method = c("gt", "dcf_clust_gtk")
    ) %>% 
    merge(runs)  %>%
    mutate_all(factor)   %>%
  mutate(folder = glue("{clust_method}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))

# dcf_vals = c("post_dcfs", "dcfs"),
pharm_res_old <-readRDS("data/pharm_res_5-16-2024.rds")

pharm_res <- bind_rows(lapply(pharm_runs$folder,  read_csv_wrapper, "pharming", "metrics.csv")) %>%
mutate(method="Pharming") %>% inner_join(pharm_runs)

dcf_clust_old <- filter(pharm_res_old, clust_method =="dcf_clustering", 
                        dcf_vals=="dcfs", snvs==5000, mclust==5) %>%
  mutate(verison = "old")

dcf_clust_old %>% group_by(folder) %>% count()
dcf_clust_new <-  filter(pharm_res, clust_method =="dcf_clust_gtk", snvs==5000, lamb==1000) %>%
  mutate(verison = "new")

dcf_clust_comp.long <- bind_rows(dcf_clust_old, dcf_clust_new) %>%   
  pivot_longer(cols=c(contains("sum_score"), "cna_tree_accuracy", "snv_tree_accuracy", "genotype_similarity",
                                                                                    "cna_genotype_similarity")) %>%
  inner_join(col_to_name_map)

box_plot_dcf_comp <- ggplot(dcf_clust_comp.long, aes(x=clust_method, y=value)) + 
  facet_grid(full_name ~cov, scales="free") + geom_boxplot() +
  xlab("method") +
  scale_x_discrete(labels=c("CN"="Phertilizer", 
                            "dcf_clust_gtk"= "DCF Clust (GT k)",
                            "gt" = "GT DCFs")) + xvert

line_plot_dcf_comp <- ggplot(dcf_clust_comp.long, aes(x=clust_method, y=value, color=s)) + 
  facet_grid(full_name ~cov) + geom_point() + geom_line(aes(group=s)) +
  xlab("method") +
  scale_x_discrete(labels=c("CN"="Phertilizer", 
                            "dcf_clust_gtk"= "DCF Clust (GT k)",
                            "gt" = "GT DCFs")) + xvert

ggsave(file.path(figpath, "dcf_comparison_m5000_l5_line.pdf"), plot=line_plot_dcf_comp, height=10, width=8)
ggsave(file.path(figpath, "dcf_comparison_m5000_l5_box.pdf"), plot=box_plot_dcf_comp, height=10, width=8)
pharm_res <- compute_acc(pharm_res)
head(pharm_res)
ggplot(pharm_res, aes(x=clust_method, y=mut_sum_score, fill=lamb)) + geom_boxplot() +
  facet_wrap(~cov)
temp <- bind_rows( pharm_res, compute_acc(phert)) %>% filter(snvs==5000, mclust==5, cells==1000,
                                                              clust_method %in% c("gt", "dcf_clust_gtk", "CN"))
temp.long <- temp %>% 
  pivot_longer(cols=c(contains("sum_score"), "cna_tree_accuracy", "snv_tree_accuracy", "genotype_similarity",
                      "cna_genotype_similarity")) %>%
                  inner_join(col_to_name_map)


lamb_comp <- ggplot(temp.long, aes(x=clust_method, y=value, fill=lamb)) + 
  facet_grid(full_name ~cov, scales="free") + geom_boxplot() +
  xlab("method") +
  scale_x_discrete(labels=c("CN"="Phertilizer", 
                            "dcf_clust_gtk"= "DCF Clust (GT k)",
                            "gt" = "GT DCFs")) + xvert
ggsave(file.path(figpath, "lambda_comparison_m5000_l5.pdf"), plot=lamb_comp, height=10, width=8)
pharm_res.long <- pharm_res %>% 
  pivot_longer(c(contains("sum_score"), "cna_tree_accuracy", "genotype_similarity",
                 "genotype_pres_similarity"))


pharm_res <- pharm_res %>% filter(lamb==1000)

pharm_res.long <- pharm_res %>% 
  pivot_longer(cols=c(contains("sum_score"), "cna_tree_accuracy", "snv_tree_accuracy", "genotype_similarity",
                      "cna_genotype_similarity")) %>%
  inner_join(col_to_name_map)

ggplot(pharm_res.long, aes(x=clust_method, y=value)) + 
  facet_grid(full_name ~cov, scales="free") + geom_boxplot() +
  xlab("clustering method")

p <- ggplot(pharm_res.long, aes(x=clust_method, y=value, color=s)) + 
  facet_grid(full_name ~cov, scales="free") + geom_point() +
  geom_line(aes(group=s))
  xlab("clustering method")
  
  ggsave(file.path(figpath, "line_comp_gt_dcf_clustgtk.pdf"), plot=p, height=10, width=8)

ggplot(pharm_res.long, aes(x=clust_method, y=value, fill=lamb)) + 
  facet_grid(full_name ~cov, scales="free") + geom_boxplot() +
  xlab("clustering method")
temp2.long <- temp %>% 
  pivot_longer(c(contains("mut_precision")))

ggplot(temp.long, aes(x=clust_method, y=value, fill=method)) + 
         facet_grid(full_name ~cov, scales="free") + geom_boxplot() +
        xlab("clustering method")

ggplot(temp.long, aes(x=clust_method, y=value, color=s)) + 
  facet_grid(name ~cov, scales="free") + geom_point() + geom_line(aes(group=s)) +
  xlab("clustering method")
ggplot(temp.long %>% filter(method  != "Phertilizer"), aes(x=factor(isegs), y=value, color=s)) + 
  facet_grid(cov~prop_name) + geom_point() + geom_line(aes(group=s)) +
  xlab("initial segments")
ggplot(temp.long, aes(x=factor(isegs), y=value, color=s)) + 
  facet_grid(cov~prop_name, scales="free") + geom_point() + geom_line(aes(group=s)) +
  xlab("initial segments")
colnames(temp)
pharm_res <- compute_acc(pharm_res)
View(pharm_res)
saveRDS(pharm_res2, "data/pharm_res_5-16-2024.rds")
pharm_res2 <- compute_acc(pharm_res)
  
 res <- bind_rows(pharm_res, phert) %>%
   mutate(clust_method = factor(clust_method, levels= c("decifer", "dcf_clustering", "gt", "CN"), labels=c("Decifer", "DCF Clustering", "GT", "CN")),
   dcf_vals = factor(dcf_vals, levels=c("dcfs", "post_dcfs", "UMAP"), labels=c("DCFs", "Post DCFs", "UMAP"))) %>%
   unite(method_name, c(clust_method, dcf_vals, method), sep="+", remove=FALSE)



res <- res %>%
  mutate(anc_mut_pairs = anc_mut_TP + anc_mut_FN,
         inc_mut_pairs = inc_mut_TP + anc_mut_FN,
         mut_pairs = mut_TP + mut_FN,
         mut_total = anc_mut_pairs + inc_mut_pairs + mut_pairs,
         anc_mut_wt = anc_mut_pairs/mut_total,
         inc_mut_wt = inc_mut_pairs/mut_total,
         mut_wt = mut_pairs/mut_total,
         mut_sum_score = mut_wt*mut_recall + anc_mut_wt*anc_mut_recall + inc_mut_wt*inc_mut_recall,
         anc_cell_pairs = anc_cell_TP + anc_cell_FN,
         inc_cell_pairs = inc_cell_TP + inc_cell_FN,
         cell_pairs = cell_TP + cell_FN,
         cell_total = anc_cell_pairs + inc_cell_pairs + cell_pairs,
         anc_cell_wt = anc_cell_pairs/cell_total,
         inc_cell_wt = inc_cell_pairs/cell_total,
         cell_wt = cell_pairs/cell_total,
         cell_sum_score = anc_cell_wt*anc_cell_recall + inc_cell_wt*inc_cell_recall + cell_wt*cell_recall
         )


df.tree.mut<- res %>%
  select(s, cells, snvs, mclust,
                                          cov,
                                          method_name,
                                          method,
                                          anc_mut_recall,
                                          inc_mut_recall,
                                          mut_recall,
                                          mut_sum_score,


         ) %>%
  pivot_longer(c("anc_mut_recall", "inc_mut_recall", "mut_recall", "mut_sum_score")) %>%
  inner_join(col_to_name_map) %>%
  mutate(type="SNV")

col_to_cell_map <- data.frame(name=c("anc_cell_recall", "inc_cell_recall", "cell_recall", "cell_sum_score"),
                              full_name = metric_names)
df.tree.cell<- res %>%
  select(s, cells, snvs, mclust, cov,
        method_name,
         method,
         anc_cell_recall,
         inc_cell_recall,
         cell_recall,
         cell_sum_score) %>%
  pivot_longer(c("anc_cell_recall", "inc_cell_recall", "cell_recall", "cell_sum_score")) %>%
  inner_join(col_to_cell_map) %>%
  mutate(type="cell")


df.tree.met <- bind_rows(df.tree.mut, df.tree.cell) %>%
  mutate(type = factor(type, levels=c("SNV", "cell"), ordered=TRUE),
         full_name = factor(full_name, levels=c( "accuracy",
                                                  "ancestral pair\n recall (APR)",
                                                  "incomparable pair\nrecall (IPR)" ,
                                                  "clustered pair\nrecall (CPR)"), ordered=T))

#########
#
# Find instances where there is a big delta between DCFs methods and GT dcfs 
########

filter(df.tree.met, full_name=="accuracy", snvs==5000, type=='SNV', cov=="0.25", s==10) %>% ggplot(aes(x=s, y=value, color=method_name)) +
  geom_jitter(size=3) + facet_wrap(~cov)


tree.plot <- ggplot(df.tree.met,
       aes(x=cov, y=value, fill=method_name)) +
  geom_boxplot() +
  facet_grid(type~full_name)   + #theme_vertx +
  # meth_col_scale + cont_y_scale +
  # myxlab +
  ylab("") + 
  scale_fill_discrete(name="") +
  theme(legend.position="top")

# res.long <- res %>% 
#    pivot_longer(c(contains("recall"), "cell_ari", "perc_cna_trees_correct"))  

ggsave(file.path(figpath, "tree_metrics_k5.pdf"), plot=tree.plot,width=10, height=10 )   



df.sim<- res %>%
  select(s, cells, snvs, mclust, cov,
         method_name,
         method,
         contains("similarity"),
         contains("accuracy")
         ) %>%
  pivot_longer(c(contains("similarity"), contains("accuracy"))) 

sim_acc_plot <- ggplot(df.sim, aes(x=cov, y=value, fill=method_name)) + geom_boxplot() +
  facet_wrap(~name) +   scale_fill_discrete(name="") +
  theme(legend.position="top")
ggsave(file.path(figpath, "similarity_k5.pdf"), plot=sim_acc_plot,width=10, height=10 )  
metric_mapping <- data.frame(name = unique(res.long$name), 
                name_label = metric_names) %>%
                mutate(name_label = factor(name_label, ordered=TRUE,
                 levels=metric_names))
res.long <- inner_join(res.long, metric_mapping)


  # ggplot(res, aes(x=snvs, y=inc_mut_recall, fill=method_name)) +
  # geom_boxplot() +
  #   facet_wrap(~cov) +
  #   xlab("# of SNVs") + ylab("value") +
  #   scale_fill_discrete(name="") + theme(legend.position ="top")

# for(k in c(5)){
#   p <- res.long %>% filter(mclust ==k, cells==1000) %>%
#     ggplot(aes(x=snvs, y=value, fill=method_name)) +
#   geom_boxplot() +
#     facet_grid(cov~name_label) +
#     xlab("# of SNVs") + ylab("value") +
#     scale_fill_discrete(name="") + theme(legend.position ="top")
#     ggsave(file.path( figpath, sprintf("scores_mclust%d.pdf", k)), plot=p, width=12, height=10)
# }




 

