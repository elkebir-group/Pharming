
library(yaml)
library(glue)
library(matrixStats)
library(tidyverse)
figpath <- "figures/pharming_recomb"

theme_set(theme_grey(base_size = 19))
xvert <- theme(axis.text.x = element_text(
    angle = 90,  # Rotate labels to be vertical
    hjust = 1,   # Adjust horizontal justification
    vjust = 0.5  # Adjust vertical justification
  ))

meth_cols <- c("Phertilizer"= "#339746", "Pharming"= "#974EA3", "Baseline+SCITE"="#FFA500")
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

metric_labels <- c(      "Cell placement\naccuracy",
                         "CNA tree accuracy",
                         "CNA state similarity",
                         "SNV placement\naccuracy",
                         "SNV tree accuracy",
                         "SNV genotype similarity",
                  " SNV ancestral pair\n recall (APR)",
                  "SNV incomparable pair\nrecall (IPR)" ,
                  "SNV clustered pair\nrecall (CPR)",
                  "Cell ancestral pair\n recall (APR)",
                  "Cell incomparable pair\nrecall (IPR)" ,
                  "Cell clustered pair\nrecall (CPR)",
                  "SNV presence similarity",
                  "SNV genotype allele-specific similarity"  )

metric_names <- c( "cell_sum_score",
                   "cna_tree_accuracy", 
                   "cna_genotype_similarity", 
                   "mut_sum_score",
                   "snv_tree_accuracy",
                   "genotype_similarity",
               
                 
                  
                   "anc_mut_recall", "inc_mut_recall", "mut_recall",
  "anc_cell_recall", "inc_cell_recall", "cell_recall",
  "genotype_pres_similarity","genotype_allele_spec_similarity"
)
col_to_name_map <- data.frame(name=metric_names,
                              full_name = factor(metric_labels, levels=metric_labels, ordered = TRUE))
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
read_csv_header <- function(folder ,prefix, suffix){
  fname <- file.path(prefix, folder, suffix)
  if(file.exists(fname)){
    df<- read.csv(fname)

    df <- df %>%
      mutate(folder = folder) 
  }else{
    print(glue("{fname} does not exist!"))
    df <- data.frame()
  }

return(df)
}

read_table_wrapper <- function(folder ,prefix, suffix){
  fname <- file.path(prefix, folder, suffix)
  print(fname)
  if(file.exists(fname)){
    df<- read.table(fname, header=T)
    
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
                    mutate(err = as.character(err)) %>%
                    mutate_all(factor)


phert_runs <- runs %>% 
    mutate(folder = glue("s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))
 
phert <- bind_rows(lapply(phert_runs$folder,  read_csv_wrapper, "phertilizer/sims4", "metrics.csv")) %>%
  group_by(folder) %>%
  mutate(clust_method = "CN", method="Phertilizer") %>% compute_acc() %>% inner_join(phert_runs) 


baseline <- bind_rows(lapply(phert_runs$folder,  read_csv_wrapper, "baseline/sims4", "metrics.csv")) %>%
  group_by(folder) %>%
  mutate(clust_method = "OPTICS", method="Baseline") %>% compute_acc() %>% inner_join(phert_runs) 

# saveRDS(compute_acc(phert), "data/phertilizer.rds")

phert_time <- bind_rows(lapply(phert_runs$folder,  read_table_wrapper, "phertilizer/sims3", "benchmark.log")) %>%
  group_by(folder) %>%
  mutate( method="Phertilizer") %>% rename(seconds=s) %>% inner_join(phert_runs) 




pharm_runs <-  expand.grid(
    isegs = pharming_config$ninit_segs,
    tm = pharming_config$ninit_tm,
    order = pharming_config$order,
    topn = pharming_config$topn,
    lamb = pharming_config$lamb,
    clust_method = c("gt", paste("dcf_clustk", 3:6, sep=""))
    ) %>% 
    merge(runs)  %>%
    mutate_all(factor)   %>%
  mutate(folder = glue("{clust_method}/{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))

pharm_folder <- "pharming_recomb"

pharm_res <- bind_rows(lapply(pharm_runs$folder,  read_csv_wrapper, pharm_folder, "metrics.csv")) %>%
  mutate(method="Pharming") %>% inner_join(pharm_runs) %>% compute_acc()



kfolders <- pharm_runs %>% filter(clust_method != "gt") %>% pull(folder)

model_selection <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "model_selection.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>%
  group_by(folder) %>% filter(row_number()==1)

selectedk <-model_selection %>% ungroup()  %>%
  group_by(cells, snvs, cov, s, mclust, err,  dirch) %>% 
  arrange(BIC) %>% filter(row_number()==1)

bst <- select( selectedk, folder ,BIC, ICL, kval)

pharm_model_select <- filter(pharm_res, clust_method != "gt") %>% inner_join(bst) %>% mutate(clust_method="ICL")


res <- bind_rows(phert %>% mutate(method_name="Phertilizer"),
                 baseline %>% mutate(method_name="Baseline+SCITE"),
                 pharm_model_select %>% mutate(method_name="Pharming"))


# %>%
#   mutate(method_label = sprintf("%s + %s", clust_method, method))

# res <- pharm_res %>% 
#   mutate(method_label = sprintf("%s + %s", clust_method, method))
res.long <- res %>% 
  pivot_longer(cols=c(contains("sum_score"),   "genotype_similarity")) 
res.long.no.phert <-  res %>% 
  pivot_longer(cols=c( "cna_tree_accuracy", "cna_genotype_similarity", "snv_tree_accuracy")) %>%
  filter(method !="Phertilizer")

res.long.filt <- bind_rows(res.long, res.long.no.phert) %>% inner_join(col_to_name_map)
for(e in unique(res.long.filt$err)){
  p <- ggplot(res.long.filt %>% filter(err==e),
              aes(x=cov, y=value, fill=method)) + 
    facet_wrap(~full_name) + geom_boxplot() + 
    scale_fill_manual(name="", values=meth_cols)+
    xlab("coverage")  + ylab("") + theme(legend.position = "top")
  ggsave(file.path(figpath, glue("sim_metrics_mclust5_m5000_e{e}.pdf")), plot=p, width=9, height=5)
}




res.long.filt %>% filter(err=="0.035") %>% group_by( full_name, method) %>% summarize(median=median(value))
res.long.filt %>% filter(err=="0.035", cov==0.01) %>% group_by( full_name, method) %>% summarize(median=median(value))

cell.long <- res %>% 
  pivot_longer(cols=c(contains("cell_sum_score"),   contains("cell_recall"))) %>% 
  inner_join(col_to_name_map)

for(e in unique(cell.long$err)){
  p <- ggplot(cell.long %>% filter(err==e),
              aes(x=cov, y=value, fill=method)) + 
    facet_wrap(~full_name, nrow=1) + geom_boxplot() + 
    scale_fill_manual(name="", values=meth_cols)+
    xlab("coverage")  + ylab("") + theme(legend.position = "top")
  p
  ggsave(file.path(figpath, glue("cell_acc_supp_e{e}.pdf")), plot=p, width=10.5, height=5)
}

snv.long <- res %>% 
  pivot_longer(cols=c(contains("mut_sum_score"),   contains("mut_recall"))) %>% 
  inner_join(col_to_name_map)

for(e in unique(snv.long$err)){
  p <- ggplot(snv.long %>% filter(err==e),
              aes(x=cov, y=value, fill=method)) + 
    facet_wrap(~full_name, nrow=1) + geom_boxplot() + 
    scale_fill_manual(name="", values=meth_cols)+
    xlab("coverage")  + ylab("") + theme(legend.position = "top")
  p
  ggsave(file.path(figpath, glue("mut_acc_supp_e{e}.pdf")), plot=p, width=10.5, height=5)
}


# ggplot(pharm_res, aes(x=clust_method, y=mut_sum_score, fill=lamb)) + geom_boxplot() +
#   facet_wrap(~cov)
# temp <- bind_rows( pharm_res, compute_acc(phert)) %>% filter(snvs==5000, mclust==5, cells==1000,
#                                                               clust_method %in% c("gt", "dcf_clust_gtk", "CN"))


# lamb_comp <- ggplot(temp.long, aes(x=clust_method, y=value, fill=lamb)) + 
#   facet_grid(full_name ~cov, scales="free") + geom_boxplot() +
#   xlab("method") +
#   scale_x_discrete(labels=c("CN"="Phertilizer", 
#                             "dcf_clust_gtk"= "DCF Clust (GT k)",
#                             "gt" = "GT DCFs")) + xvert
# ggsave(file.path(figpath, "lambda_comparison_m5000_l5.pdf"), plot=lamb_comp, height=10, width=8)
# pharm_res.long <- pharm_res %>% 
#   pivot_longer(c(contains("sum_score"), "cna_tree_accuracy", "genotype_similarity",
#                  "genotype_pres_similarity"))


# pharm_res <- pharm_res %>% filter(lamb==1000)
# 
# pharm_res.long <- pharm_res %>% 
#   pivot_longer(cols=c(contains("sum_score"), "cna_tree_accuracy", "snv_tree_accuracy", "genotype_similarity",
#                       "cna_genotype_similarity")) %>%
#   inner_join(col_to_name_map)



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

########Loss Analysis#################

loss_error_report <- bind_rows(lapply(phert_runs$folder,  read_csv_header, "sims4", "loss_error_report.csv")) 


loss_results <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "loss_scores.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", "")))


loss_results <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "loss_scores.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% inner_join(bst) %>%
  filter(solution==0)  %>%
  mutate(F1= 2*((precision*recall)/(precision+recall)),
         cna_f1= 2*((cna_precision*cna_recall)/(cna_precision + cna_recall)))

loss_position <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "fps_leaves.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% inner_join(bst)

bst_folder <- select(bst, folder)

total_inf_loss <- loss_results %>% select(folder, num_inf_lost)

total_fps <- loss_position %>% filter(TP=='False') %>% group_by(folder) %>% summarize(num_fps = n())
fps <- loss_position %>% filter(TP=='False') %>% group_by(folder, is_leaf) %>%
  count() %>% 
  inner_join(total_fps) %>% mutate(perc= n/num_fps)


fps_bst <- fps %>% inner_join(bst_folder)
ggplot(fps  %>% inner_join(pharm_runs), aes(x=cov, y=perc, fill=is_leaf)) + 
  geom_boxplot() + facet_wrap(~err) + ylab("Percent of FPs") + xlab("coverage")
# loss_pos_sum <- loss_position %>% group_by(folder, is_leaf, TP, cov, err) %>% count()   %>% 
#   inner_join(total_inf_loss) %>% mutate(perc = n/num_inf_lost)


inf_lost_snvs <- select(loss_position, snv, TP, folder)
######cmb###
cmb <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "cmb.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% 
  inner_join(bst %>% select(-BIC,-ICL))  



ggplot(cmb, aes(x=cov, y=cmb, fill=factor(within_clade))) + geom_boxplot() +
         facet_wrap(~err) +scale_fill_discrete(name="", labels=c("outside clade", "within clade")) +
ylab("gain CMB")

#####cmb loss####
cmb_loss <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "loss_cmb.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% 
  inner_join(bst %>% select(-BIC,-ICL))


  
  
  
  ggplot(cmb_loss, aes(x=cov, y=cmb, fill=factor(within_clade))) + geom_boxplot() +
  facet_wrap(~err) +scale_fill_discrete(name="", labels=c("within lost clade")) +
  ylab("gain CMB") + ylab("loss CMB (WIHIN CLADE ONLY)")
  
  #####vaf validation####
  vaf_val <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "vaf_validation.csv")) %>%
    inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% 
    inner_join(bst %>% select(-BIC,-ICL))
  
  
  filter(vaf_val, total >= 10, err==0) %>% left_join(inf_lost_snvs) %>% ggplot( aes(x=obs_vaf, fill=TP)) + geom_histogram() +
    geom_vline(aes(xintercept=latent_vaf)) + 
    facet_wrap(~latent_vaf, scales="free") + xlab("observed VAF") +
    scale_fill_discrete("SNV loss TP?")
  
  filter(vaf_val, total >= 10, err==0) %>% left_join(inf_lost_snvs) %>%
    filter(!is.na(TP)) %>%
    ggplot( aes(x=obs_vaf, fill=TP)) + geom_histogram() +
    geom_vline(aes(xintercept=latent_vaf)) + 
    facet_wrap(~latent_vaf, scales="free") + xlab("observed VAF") +
    scale_fill_discrete("SNV loss TP?")

  
#### snv tree likelihood comparison ####
 snv_like <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "snv_likelihood_comp.csv")) %>%
    inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% 
    inner_join(bst %>% select(-BIC,-ENT,-ICL))
  
snv_like_wide <- snv_like %>% select(-tree_edges, -cluster) %>%
  pivot_wider(names_from="source", values_from="cost") %>%
  mutate(diff=inf-gt) %>% left_join(inf_lost_snvs)

ggplot(snv_like_wide, aes(x=err, y = diff, fill=TP)) + geom_boxplot() + facet_wrap(~cov, scales="free_y") +
  xlab("coverage") + ylab("inferred SNV likelihood - gt")+ scale_fill_discrete("SNV loss TP?")

ggplot(snv_like_wide %>% filter(!is.na(TP)), aes(x=err, y = diff, fill=TP)) + geom_boxplot() + facet_wrap(~cov, scales="free_y") +
  xlab("coverage") + ylab("inferred SNV likelihood - gt")+ scale_fill_discrete("SNV loss TP?")

snv_tree_costs <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "snv_tree_assignment_costs.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% inner_join(bst) 


snv_tree_costs <- snv_tree_costs %>% mutate(snvs = as.numeric(as.character(snv))) %>%
  select(snvs, cluster, tree_id, tree_edges, is_lost, cost, assigment, assigned_tree, folder )

snv_tree_min <- snv_tree_costs %>% group_by(folder, snvs, cluster) %>% arrange(cost) %>% filter(row_number()==1) %>% 
  arrange(folder, snvs) %>% mutate(cost= -1*cost)

denom <- snv_tree_min %>% group_by(folder, snvs) %>% summarize(log_denom =logSumExp(cost))

snv_tree_min <- snv_tree_min %>%
  left_join(denom, by=c("folder", "snvs"))


snv_tree_min <- snv_tree_min %>%
  mutate(posterior = ifelse(posterior == 0, .Machine$double.eps, posterior))  # Avoid log(0)

# Step 2: Compute entropy for each group (folder, snvs)
entropy_result <- snv_tree_min %>%
  group_by(folder, snvs) %>%
  summarize(entropy = -sum(posterior * log(posterior))) %>% 
  inner_join(pharm_runs %>% rename(nsnvs= snvs))


ggplot(entropy_result, aes(x=cov, y=entropy, fill=err)) + geom_boxplot()
ggplot(entropy_result %>% filter(err==0, cov==0.01), aes(x=entropy, color=factor(s))) + geom_density() + facet_wrap(~cov, scales="free") +
  scale_color_discrete(name="seed")

snv_losses <- loss_position %>% select(node, snvs =snv, is_leaf, TP, folder)

snv_tree_cost_losses <- inner_join(snv_losses, snv_tree_costs)



snv_tree_fp <- snv_tree_cost_losses %>% filter(TP == "False") %>% group_by(folder, snvs, is_lost) %>% arrange(cost) %>% filter(row_number()==1)

snv_tree_tp <- snv_tree_cost_losses %>% filter(TP == "True") %>% group_by(folder, snvs, is_lost) %>% arrange(cost) %>% filter(row_number()==1)
snv_tree_comp<- snv_tree_cost_losses  %>% group_by(folder, TP, snvs, is_lost) %>% arrange(cost) %>% filter(row_number()==1)

snv_tree_wide  <- snv_tree_comp %>%  select(-tree_id, -tree_edges, -cluster) %>%
  pivot_wider( names_prefix = "lost_", names_from = "is_lost", values_from="cost") %>% rename(snv= snvs) %>%left_join(pharm_runs) %>%
  mutate(diff= lost_False  - lost_True)

# snv_tree_tp_wide  <- snv_tree_fp %>%  select(-tree_id, -tree_edges, -cluster) %>%
#   pivot_wider( names_prefix = "lost_", names_from = "is_lost", values_from="cost") %>% rename(snv= snvs) %>%left_join(pharm_runs) %>%
#   mutate(diff= lost_False  - lost_True)

filt_loss <- bind_rows(lapply(kfolders,  read_csv_header, pharm_folder, "filtered_loss_scores.csv")) %>%
  inner_join(pharm_runs) %>% mutate(kval = as.numeric(str_replace(clust_method, "dcf_clustk", ""))) %>% inner_join(bst) 

filt_loss %>% summarize(perc_no_lost_snvs =sum(num_gt_lost==0)/n())
head(filt_loss)
ggplot(snv_tree_wide, aes(x=err, y=diff, fill=TP)) + geom_boxplot() + 
  facet_wrap(~cov, scales="free_y") + xlab("CNA error rate") +
  ylab("lost SNV likelihood difference") + 
  scale_fill_discrete(name="inferred loss", labels=c("FP", "TP")) 

ggplot(loss_pos_sum, aes(x=is_leaf, y=perc, fill=TP)) + geom_boxplot() + facet_grid(cov~err)

pharm_time <- bind_rows(lapply(pharm_runs$folder,  read_table_wrapper, pharm_folder, "benchmark.log")) %>%
  group_by(folder) %>%
  mutate( method="Pharming") %>% rename(seconds=s) %>% inner_join(pharm_runs) %>% 
  inner_join(bst %>% select(folder))

ggplot(loss_results %>% select(-BIC,-ICL) %>%
         pivot_longer(c("precision", "recall", "F1")), aes(x=cov, y=value, fill=err)) +
  geom_boxplot() + facet_wrap(~name) + xlab("coverage") +
  scale_fill_discrete(name="CNA error")

filt_loss <- mutate(filt_loss, F1= 2*((precision*recall)/(precision+recall))) 
ggplot(filt_loss %>% select(-BIC,-ICL) %>% filter(num_gt_lost>0) %>%
         pivot_longer(c("precision", "recall", "F1")), aes(x=cov, y=value, fill=err)) +
  geom_boxplot() + facet_wrap(~name) + xlab("coverage") +
  scale_fill_discrete(name="CNA error")

ggplot(loss_results %>% select(-BIC,-ENT,-ICL) %>%
         pivot_longer(c("cna_precision", "cna_recall", "cna_f1", "cna_acc")), aes(x=cov, y=value, fill=err)) +
  geom_boxplot() + facet_wrap(~name, nrow=1) + xlab("coverage") +
  scale_fill_discrete(name="CNA error")

time <- bind_rows(phert_time, pharm_time)
time_plot <- ggplot(time, aes(x=method, y=seconds/60, fill=method)) + geom_boxplot() + scale_y_log10() + 
  annotation_logticks(sides="l") + ylab("running time (min)") +
  scale_fill_manual(name="", values=meth_cols) +guides(fill="none")
time %>% group_by(method) %>% summarize(median=median(seconds/60))
time_plot
ggsave(file.path(figpath, glue("running_time.pdf")), plot=time_plot, width=3.5, height=4.5)




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



 

