library(tidyverse)
library(yaml)
library(glue)
figpath <- "figures"
#assuming working directory is /scratch/data/leah/Pharming/simulation_study
# metric_names <-  c("Ancestral Pair\nRecall", "Incomparable Pair\nRecall",
#                 "Clustered Pair\nRecall", "Cell\nARI", "CNA Tree\nAccuracy")

metric_names <- c("ancestral pair\n recall (APR)",
                  "incomparable pair\nrecall (IPR)" ,
                  "clustered pair\nrecall (CPR)",
                  "accuracy")
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

pharm_runs <-  expand.grid(
    isegs = pharming_config$ninit_segs,
    tm = pharming_config$ninit_tm,
    order = pharming_config$order,
    topn = pharming_config$topn,
    lamb = pharming_config$lamb,
    dcf_vals = c("post_dcfs", "dcfs"),
    clust_method = c("decifer", "dcf_clustering", "gt")
    ) %>% 
    mutate_all(factor)  %>%
merge(runs) %>%
  mutate(folder = glue("{clust_method}/{dcf_vals}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))




pharm_res <- bind_rows(lapply(pharm_runs$folder,  read_csv_wrapper, "pharming", "metrics.csv")) %>%
mutate(method="Pharming") %>% inner_join(pharm_runs)



  
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
col_to_name_map <- data.frame(name=c("anc_mut_recall", "inc_mut_recall", "mut_recall", "mut_sum_score"),
                              full_name = metric_names)

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




 

