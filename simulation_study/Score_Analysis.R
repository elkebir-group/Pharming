library(tidyverse)
library(yaml)
library(glue)
library(mclust)
bpath <- "/Users/leah/Documents/Research/projects/Pharming/simulation_study"
config = yaml.load_file(file.path(bpath, "config.yml"))

runs <- expand.grid(cells=config$cells, snvs=config$snvs, cov=c("0.01", "0.05",  "0.25", "1"), 
                    mclust = config$mclust,
                    err = config$cerror, s = 11, nsegs=config$nsegs)%>% 
  mutate(folder = glue("s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}"))




read_csv_wrapper <- function(folder, prefix, suffix){
  df<- read.csv(file.path(bpath, prefix, folder, suffix)) %>%
    mutate(folder = folder)
  return(df)
}



res <- bind_rows(lapply(runs$folder,  read_csv_wrapper, "pharming", "scores.csv")) %>%
  group_by(folder) %>%
  inner_join(runs) %>% mutate(cell_assign="greedy")

res_ilp <- bind_rows(lapply(runs$folder,  read_csv_wrapper, "pharming_ilp", "scores.csv")) %>%
  group_by(folder) %>%
  inner_join(runs) %>% mutate(cell_assign="using_dcfs")

res <- bind_rows(res, res_ilp)

  
res %>%  filter(cov %in% c("0.01", "0.05", "0.25")) %>%
   pivot_longer(c(contains("recall"), "cell_ari"))  %>% 
  ggplot(aes(x=cell_assign, y=value)) +
  geom_point() + 
  geom_line(aes(group=cov)) +
  facet_grid(cov~name, scales="free_y") +
  xlab("coverage") + ylab("value") 

ggplot(res, aes(x=gt_cost, y= cost, color=factor(s), shape=factor(nsegs))) + 
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  geom_point() + 
  facet_wrap(~cov, scales="free", labeller = "label_both") + 
  scale_color_discrete(name="sim instance") +
  scale_shape_discrete(name="# of segments") +
  xlab("Ground truth cost") + ylab("Inferred cost")


filter(res)  %>% View()
res %>% pivot_longer(c(contains("recall"), "cell_ari"))  %>% 
  ggplot(aes(x=factor(cov), y=value)) +
  geom_boxplot() + 
  facet_grid(nsegs~name) +
  xlab("coverage") + ylab("value") 

res   %>% 
  ggplot(aes(x=factor(cov), y=cna_mad)) +
  geom_boxplot() + 
  xlab("coverage") + ylab("CNA Mean Total Allele Deviation per Cell") 

