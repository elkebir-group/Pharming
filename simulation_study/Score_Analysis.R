library(tidyverse)
library(yaml)
library(glue)
figpath <- "simulation_study/figures"
metric_names <-  c("Ancestral Pair\nRecall", "Incomparable Pair\nRecall",
                "Clustered Pair\nRecall", "Cell\nARI", "CNA Tree\nAccuracy")




read_csv_wrapper <- function(folder, prefix, suffix){
  fname <- file.path(bpath, prefix, folder, suffix)
  if(file.exists(fname)){
    df<- read.csv(fname) %>%
    mutate(folder = folder) 
  }else{
    print(glue("{fname} does not exist!"))
    df <- data.frame()
  }

  return(df)
}

bpath <- "/scratch/leah/Pharming/simulation_study"
config = yaml.load_file(file.path(bpath, "simulate.yml"))
pharming_config <- yaml.load_file(file.path(bpath, "pharming.yml"))
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


pharm_runs <-  expand.grid(
    isegs = pharming_config$ninit_segs,
    tm = pharming_config$ninit_tm,
    order = pharming_config$order,
    topn = pharming_config$topn,
    lamb = pharming_config$lamb) %>% 
    mutate_all(factor)  %>%
merge(runs) %>%
  mutate(folder = glue("{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))


phert_runs <- runs %>% 
    mutate(folder = glue("s{s}_m{snvs}_k{segs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))

clust_methods <- c("decifer", "dcf_clust", "sims")
pharm_res <- data.frame()
for(clust in clust_methods){
  temp <- bind_rows(lapply(pharm_runs$folder,  read_csv_wrapper, glue("pharming/{clust}"), "scores.csv")) %>%
  group_by(folder) %>%
  mutate(clust_method = clust )
  pharm_res <- bind_rows(pharm_res, temp)
}

phert <- bind_rows(lapply(phert_runs$folder,  read_csv_wrapper, "phertilizer/sims", "scores.csv")) %>%
  group_by(folder) %>%
  mutate(clust_method = "Phertilizer", segment=as.character(segment)) %>% inner_join(phert_runs)


top_pharm  <- pharm_res %>% group_by(folder, clust_method) %>% top_n(1, -inf_cost) %>% 
  inner_join(pharm_runs) 

top_pharm %>% group_by(folder, clust_method) %>% count()
  
 res <- bind_rows(top_pharm ,phert)%>%
   mutate(clust_method = factor(clust_method, levels= c(clust_methods, "Phertilizer"), 
   labels=c("Decifer + Pharming", "DCF Clustering +Pharming", "Ground Truth DCFs + Pharming", "Phertilizer")))



res.long <- res %>% 
   pivot_longer(c(contains("recall"), "cell_ari", "perc_cna_trees_correct"))  

metric_mapping <- data.frame(name = unique(res.long$name), 
                name_label = metric_names) %>%
                mutate(name_label = factor(name_label, ordered=TRUE,
                 levels=metric_names))


 foo <- res.long %>% filter(mclust ==k, cells==1000)

for(k in c(5)){
  p <- res.long %>% filter(mclust ==k, cells==1000) %>%
  inner_join(metric_mapping) %>%
    ggplot(aes(x=snvs, y=value, fill=clust_method)) +
  geom_boxplot() +
    facet_grid(cov~name_label) +
    xlab("# of SNVs") + ylab("value") +
    scale_fill_discrete(name="") + theme(legend.position ="top")
    ggsave(file.path( figpath, sprintf("scores_mclust%d.pdf", k)), plot=p, width=12, height=10)
}




 

