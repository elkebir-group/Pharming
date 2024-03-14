library(tidyverse)
library(yaml)
library(glue)
library(mclust)
bpath <- "/Users/leah/Documents/Research/projects/Pharming/simulation_study"
config = yaml.load_file(file.path(bpath, "config.yml"))

runs <- expand.grid(cells=config$cells, snvs=config$snvs, cov=c("0.01", "0.05", "0.25", "1"), 
                    mclust = config$mclust,
                    err = config$cerror, s = 10:14, nsegs=config$nsegs)%>% 
  mutate(folder = glue("s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}"))

"decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/input.tsv"

get_decifer<- function(folder, prefix, suffix){
  dec.out <- read.table(file.path(bpath, prefix, folder, suffix), header=TRUE)
  dec.clust <- select(dec.out, snv=mut_index ,cluster, 
                      inf_dcf = point_estimate_DCF0,  true_cluster_DCF0) %>%
    separate(col=true_cluster_DCF0, c("clust_dcf", "CI"), sep=";") %>%
    mutate(clust_dcf= as.numeric(clust_dcf))
  dec.clust$folder <- folder
  return(dec.clust)
}

read_csv_wrapper <- function(folder, prefix, suffix){
  df<- read.csv(file.path(bpath, prefix, folder, suffix)) %>%
    mutate(folder = folder)
  return(df)
}



res <- bind_rows(lapply(runs$folder,  get_decifer, "decifer", "decifer_output.tsv"))
psi <- bind_rows(lapply(runs$folder,  read_csv_wrapper ,"decifer", "gt_psi.csv"))
delta <- bind_rows(lapply(runs$folder,  read_csv_wrapper ,"decifer", "gt_delta.csv"))

all.res  <- inner_join(psi, delta) %>% filter(dcf < 1) %>%
  inner_join(res)


head(all.res)
ari.df <- all.res %>% group_by(folder) %>% summarize(ari= adjustedRandIndex(node, cluster)) %>% 
  inner_join(runs)

mad.df <- all.res %>% mutate(abs_diff = abs(clust_dcf - dcf) ) %>% 
  group_by(folder) %>% summarize(MAD = mean(abs_diff)) %>% 
  inner_join(runs)

all.res %>% group_by(folder) %>%
  summarize(gt_clust = n_distinct(dcf), inf_clust = n_distinct(clust_dcf)) %>%
  pivot_longer(contains("clust")) %>%
  inner_join(runs)  %>%
  ggplot(aes(x=cov, y=value, fill=name)) + geom_boxplot() +
  xlab("coverage") + ylab("number of clusters")

all.res %>% mutate(abs_diff = abs(clust_dcf - dcf) ) %>% inner_join(runs) %>%
  ggplot( aes(x=cov, y=abs_diff)) + geom_boxplot() + 
  facet_wrap(~nsegs, labeller = "label_both") + xlab("coverage") +
  ylab("absolute difference (gt dcf - clust dcf)")

ggplot(mad.df, aes(x=cov, y=MAD)) + geom_boxplot() + 
  facet_wrap(~nsegs, labeller = "label_both") + xlab("coverage") +
  ylab("mean absolute difference (gt dcf - clust dcf)")

ggplot(ari.df, aes(x=cov, y=ari)) + geom_boxplot() + facet_wrap(~nsegs, labeller = "label_both") +
  xlab("coverage")


