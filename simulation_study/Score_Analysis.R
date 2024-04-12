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
                    nsegs=config$nsegs,
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
  mutate(folder = glue("{order}/isegs{isegs}_tm{tm}_top{topn}_lamb{lamb}/s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))



pharm_res <- bind_rows(lapply(pharm_runs$folder,  read_csv_wrapper, "pharming/sims", "scores.csv")) %>%
  group_by(folder) %>%
  inner_join(pharm_runs)

top_pharm  <- pharm_res %>% group_by(folder) %>% top_n(1, -inf_cost) 

top_pharm %>% group_by(folder) %>% count()
top_pharm.long <- top_pharm %>% 
   pivot_longer(c(contains("recall"), "cell_ari", "perc_cna_trees_correct"))  

metric_mapping <- data.frame(name = unique(top_pharm.long$name), 
                name_label = metric_names) %>%
                mutate(name_label = factor(name_label, ordered=TRUE,
                 levels=metric_names))


for(k in c(5,7)){
  p <- top_pharm.long %>% filter(mclust ==k) %>%
  inner_join(metric_mapping) %>%
    ggplot(aes(x=snvs, y=value, fill=cells)) +
  geom_boxplot() +
    facet_grid(cov~name_label) +
    xlab("# of SNVs") + ylab("value") 
    ggsave(file.path( figpath, sprintf("scores_mclust%d.pdf", k)), plot=p)
}




 

