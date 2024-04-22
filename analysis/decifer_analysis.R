library(tidyverse)
library(yaml)
library(glue)

bpath <- "simulation_study"
figs <- file.path(bpath, "figures")
config = yaml.load_file(file.path(bpath, "simulate.yml"))
theme_set(theme_gray(base_size = 20))
seeds <- 10:14
runs <- expand.grid(cells=config$cells, 
                    snvs=config$snvs,
                     cov= config$cov,
                    mclust = config$mclust,
                    err = config$cerror, 
                    s = seeds, 
                    nsegs=config$nsegs,
                    dirch = config$dirch) %>%
                    mutate_all(factor) %>%
                     mutate(folder = glue("s{s}_m{snvs}_k{nsegs}_l{mclust}_d{dirch}/n{cells}_c{cov}_e{err}"))


#"decifer/s{s}_m{snvs}_k{nsegs}_l{mclust}/n{cells}_c{cov}_e{err}/input.tsv"

get_decifer<- function(folder, prefix, suffix){
  fname <- file.path(bpath, prefix, folder, suffix)
  if(file.exists(fname)){
      dec.out <- read.table(fname, header=TRUE)
      dec.clust <- select(dec.out, snv=mut_index ,cluster, 
                          inf_dcf = point_estimate_DCF0,  true_cluster_DCF0) %>%
        separate(col=true_cluster_DCF0, c("clust_dcf", "CI"), sep=";") %>%
        mutate(clust_dcf= as.numeric(clust_dcf))
      dec.clust$folder <- folder
      return(dec.clust)

  }else{
    return(data.frame())
  }

}

read_csv_wrapper <- function(folder, prefix, suffix){
  df<- read.csv(file.path(bpath, prefix, folder, suffix), header=F, col.names=c("gt.dcf")) %>%
    mutate(folder = folder, cluster= row_number())
  return(df)
}

read_scores <- function(folder, prefix, suffix){
  df<- read.csv(file.path(bpath, prefix, folder, suffix)) %>%
    mutate(folder = folder)
  return(df)
}

dcf_clust  <-  bind_rows(lapply(runs$folder, read_scores, "dcf_clustering_v2", "scores.csv")) %>% inner_join(runs)

clust_comp <-  bind_rows(lapply(runs$folder, read_scores, "dcf_clustering_v2", "comp.csv")) %>% inner_join(runs)

numk <- clust_comp %>% group_by(folder) %>% count()

clust_comp <- clust_comp %>% inner_join(numk)

clust_plot <- ggplot(clust_comp, aes(x=gt, y=inf, color=factor(n))) + 
geom_point(size=3) + facet_wrap(~cov, ncol=1) +
geom_abline(linetype="dashed") +
scale_color_discrete(name="k") +xlab("ground truth DCFs") +
ylab("inferred DCFs")
clust_plot
ggsave(file.path(figs, "dcf_clust_scatterplot.pdf"), plot=clust_plot)
dcf_clust.long <- dcf_clust %>% pivot_longer(c("mad", contains("diff"), "perc_cna_tree_correct"))

p <- dcf_clust.long %>% ggplot(aes(x=factor(cov), y=value)) + 
geom_boxplot()+ facet_wrap(~name, scale="free_y") +xlab("coverage")

ggsave(file.path(figs, "dcf_clustering.pdf"), plot=p)


# ggplot(dcf_clust, aes(x=factor(cov), y=mad)) + geom_boxplot() + xlab("coverage")
# ggplot(dcf_clust, aes(x=factor(cov), y=max_diff)) + geom_boxplot() + xlab("coverage")
# ggplot(dcf_clust, aes(x=factor(cov), y=perc_cna_tree_correct)) + geom_boxplot() + xlab("coverage")


res <- bind_rows(lapply(runs$folder,  get_decifer, "decifer", "decifer_output.tsv")) 


inf <- res %>% group_by(folder) %>% select(folder, clust_dcf ) %>% distinct() %>%
     arrange(folder, desc(clust_dcf)) %>% group_by(folder) %>% mutate(cluster = row_number())

num.inf <-inf %>% group_by(folder) %>% summarize(ninf= n())
gt <- bind_rows(lapply(runs$folder, read_csv_wrapper, "sims", "dcfs.txt")) %>% 
       group_by(folder ) %>% arrange(folder, desc(gt.dcf)) %>%
       mutate(cluster= row_number())

num.gt <- gt %>% group_by(folder) %>% summarize(ngt = n())
num.clust <- inner_join(num.inf, num.gt)


nclust.plot  <- num.clust %>% group_by(ngt, ninf) %>% count() %>%
 ggplot(aes(x=factor(ngt), y=factor(ninf), fill=n)) + geom_tile() +
  xlab("number of gt cluster") + ylab("number of decifer inferred clusters") +
  geom_text(aes(label=n), size=15, color="white")
nclust.plot
sp <- function(fname, p, h=6, w=8){
  ggsave(file.path(figs, fname), plot=p, heigh=h, width=w)
}

sp("numclust.pdf", nclust.plot)


correct_num <- filter(num.clust, ngt==ninf) %>% pull(folder)

correct.res <- inf %>% filter(folder %in% correct_num) %>% 
inner_join(filter(gt, folder %in% correct_num)) %>%
 mutate(abs_diff = abs(clust_dcf- gt.dcf)) %>% inner_join(runs)

correct.clust.scat <- ggplot(correct.res, aes(x=clust_dcf, y=gt.dcf, color=factor(cov))) + 
geom_point(size=3) + geom_abline() + theme(legend.position="top") +
scale_color_discrete(name="coverage") + xlab("decifer DCF") + ylab("gt DCF")

sp("correct.clust.scat.pdf", correct.clust.scat, h=6, w=6)






