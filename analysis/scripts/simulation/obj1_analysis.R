library(yaml)
library(glue)
library(tidyverse)
spath <- "../simulation_study"
config <- yaml.load_file(file.path(spath, "config.yml"))
seeds <- 10:14

read_csv_wrapper <- function(path, fname){
  file <- file.path(path, fname)
  if(file.exists(file)){
    return(read_csv(file) %>% mutate(pth =path))
  }else{
    return(data.frame(pth=path))
  }
}

cov_vals <- c("0.01", "0.05", "0.1", "0.25", "1")
lamb_vals <- c(0, 0.1, 1, 5)
all_sims <- expand.grid(seed = seeds, m=config$snvs, k=config$nsegs, mclust = config$mclust,
                         n=config$cells,
                        cov = cov_vals, e=c("0") ) %>%
  mutate(inst_folder =glue("s{seed}_m{m}_k{k}_l{mclust}"),
         cell_folder =glue("n{n}_c{cov}_e{e}"),
         pth = file.path(spath, "input", inst_folder, cell_folder)) %>%
  mutate(n=factor(n), m=factor(m), cov = factor(cov, levels=cov_vals, ordered=T))
         

head(all_sims)
# snv_sims <- all_sims %>% 
#   mutate(n=factor(n), m=factor(m), cov = factor(cov, levels=cov_vals, ordered=T),
#   pth = file.path(spath, "input", inst_folder, cell_folder))
  
snv_dat <- bind_rows(lapply(all_sims$pth, read_csv_wrapper, "snv_assignment.csv")) %>%
  inner_join(all_sims) #%>% filter(!is.na(cost))

snv_dat <- snv_dat %>% filter(!is.na(approach))



snv_dat <- mutate(snv_dat, opt_tree = best_tree & (num_best_trees <= 2))

snv_sum <- snv_dat  %>% group_by(seed, n,m,k, approach, cov) %>%
  summarize(acc = sum(best_tree)/n(), acc2 = sum(opt_tree)/n(), mean=mean(num_best_trees), median=median(num_best_trees))



ggplot(snv_sum, aes(x=n, y=acc, fill=approach)) + 
  geom_boxplot() + facet_wrap(~cov, labeller="label_both") + ylab("accuracy")

ggplot(snv_sum, aes(x=n, y=acc2, fill=approach)) + 
  geom_boxplot() +  facet_wrap(~cov, labeller="label_both") +ylab("accuracy2 (# of best fitting trees <= 2)")

ggplot(snv_dat, aes(x=n, y=num_best_trees, fill=approach)) + 
  geom_boxplot() + facet_wrap(~cov, labell="label_both")


ggplot(snv_dat %>% filter(num_best_trees <= 2), aes(x=n, y=num_best_trees, fill=approach)) + 
  geom_boxplot() + facet_wrap(~cov, labeller="label_both")


snv_acc <- snv_dat %>% group_by( lamb, k, m,n,mclust, cov, e) %>% summarize( accuracy = sum(best_tree)/n(), total=n())

ggplot(snv_acc, aes(x=factor(n), y=accuracy, fill=factor(lamb))) + 
  geom_boxplot() + facet_grid(cov ~e, labeller = "label_both") + 
  xlab("number of cells") +ylab("snv assignment accuracy") +
   scale_fill_brewer(name="lamb", palette = "Set1")
  



all_sims %>% group_by(pth) %>% count() %>% filter(n > 1)
res <- bind_rows(lapply(all_sims$pth, read_csv_wrapper, "obj1_cell_accuracy.csv")) %>%
  inner_join(all_sims)
res.long <- select(res, pth, contains("perc")) %>% pivot_longer(contains("perc")) %>% inner_join(all_sims)
res %>% group_by(pth) %>% count() %>% filter(n> 1)

cna_res <- bind_rows(lapply(filter(all_sims, lamb==0) %>% pull(pth), read_csv_wrapper, "cna_only_cell_accuracy.csv")) %>%
  inner_join(all_sims)


#results  <- res

cna.res.long <-select(cna_res, pth, contains("perc")) %>% pivot_longer(contains("perc")) %>% inner_join(all_sims)
cna.res.long <- select(cna.res.long, -lamb) %>% mutate(lamb = "Inf")

res.long <- res.long %>% mutate(lamb = as.character(lamb))

all.res.long <- bind_rows(cna.res.long, res.long)  %>%
  mutate(lamb = factor(lamb, levels=c("0", "0.1", "1", "5", "Inf"), ordered = T))
library(RColorBrewer)
filter(all.res.long, name=="perc_exact")  %>% ggplot(aes(x=m, y=value, fill=lamb)) + geom_boxplot() + 
  facet_grid(cov~e,  labeller="label_both")  + ylab("cell assignment accuracy") +
  scale_fill_brewer(palette = "Set1")

recall.long <-select(res, pth, contains("recall")) %>% pivot_longer(contains("recall")) %>% inner_join(all_sims)

ggplot(res.long, aes(x=m, y=value, fill=name)) + geom_boxplot() + 
  facet_grid(cov~lamb,  labeller="label_both") +
  scale_fill_discrete(name="") + ylab("percentage of cells")
ggplot(recall.long, aes(x=m, y=value, fill=name)) + 
  geom_boxplot() + facet_wrap(~cov, nrow=1, labeller="label_both") +
  scale_fill_discrete(name="") + ylab("recall")


ggplot(cna.res.long, aes(x=m, y=value, fill=name)) + geom_boxplot() + 
  facet_grid(k~e,  labeller="label_both") +
  scale_fill_discrete(name="") + ylab("percentage of cells")

res <- mutate(res, cost2_ratio = cost2/cost2_gt, cost1_ratio = cost1/cost1_gt, perc_increase = (cost2-cost2_gt)/cost2_gt)

res.filt <- filter(res, !is.na(feature))
ggplot(res %>% filter(cov != 1), aes(x=cost2_ratio, y=incomp_pair_recall, color=n)) + geom_point() + 
  facet_wrap(~m, scales="free")

ggplot(res %>% filter(cov == 0.1), aes(x=cost1_ratio, y=perc_exact)) + geom_point() + 
  facet_wrap(~m, scales="free", labeller="label_both") + geom_smooth(method="lm") 

ggplot(res %>% filter(cov == 0.1), aes(x=cost1_ratio, y=perc_close)) + geom_point() + 
  facet_wrap(~m, scales="free", labeller="label_both") + geom_smooth(method="lm") 

ggplot(res %>% filter(cov == 0.1), aes(x=cost2_ratio, y=perc_exact)) + geom_point() + 
  facet_wrap(~m, scales="free", labeller="label_both") + geom_smooth(method="lm") 

ggplot(res %>% filter(cov == 0.1), aes(x=cost2_ratio, y=perc_close)) + geom_point() + 
  facet_wrap(~m, scales="free", labeller="label_both") + geom_smooth(method="lm") 

cor(res.filt$cost2_ratio, res.filt$perc_close)
