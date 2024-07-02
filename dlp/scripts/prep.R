library(tidyverse)
# vtext <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
theme_set(theme_gray(base_size = 20))


deep_copy_in <- snakemake@input[['copy_numbers']]
var_reads_in <- snakemake@input[['var_reads']]


var_reads_out <- snakemake@output[['var_reads']]
deep_copy_out <- snakemake@output[['copy_numbers']]
sample_segs_out <- snakemake@output[['sampled_segs']]

figpath <- snakemake@params[['figpath']]
# nsegs <- as.numeric(snakemake@params[['nsegs']])
selected_sample <- snakemake@wildcards[['sample']]
# seed <- snakemake@params[['seed']]
thresh <- snakemake@params[['thresh']]
minvar <- snakemake@params[['minvar']]


# set.seed(seed)



df <- read_csv(deep_copy_in)

colnames(df) <- c("cell", "chr", "start", "end", "x", "y")
df <- df %>% separate(col=cell, into=c("sample", "clone", "region", "cell_id"), sep="-", remove=F)  %>%
    mutate(total = x + y) 


seg_mapping <- df %>% select(chr, start, end) %>% distinct() %>% arrange(chr, start, end) %>%
  mutate(segment = row_number()-1)

df.filt <- inner_join(df %>% filter(sample %in% selected_sample) , seg_mapping)


ncells <- n_distinct((df.filt$cell))

cp <- select(df.filt, cell,chr, segment, x, y) %>% unite(col=cn, x, y, sep ="|", remove=FALSE)

num_states_df <- cp %>% group_by(segment, cn) %>%
  summarize(num_cells =n(), prop= num_cells/ncells )

nstates <- num_states_df %>% group_by(segment) %>% summarize(nstates= sum(prop >= thresh))



#write_csv(segments.df, "segments.bin.mapping.csv")

# var.dat <- read.table(var_reads_in, header = F, sep="\t",
#                       col.names=c("chr", "loci", "cell", "base", "var", "total") )

var.dat <- read_csv(var_reads_in)
var.dat  <- var.dat %>% separate(col=cell, into=c("sample", "clone", "region", "cell_id"), sep="-", remove=F)

var.dat.filt <- var.dat %>% filter(sample %in% selected_sample) %>%
 separate(col=mutation, into =c("chr", "loci")) %>% rename(var = varReads, total= totReads)

# var <- var.dat.filt %>% 
#   left_join(seg_mapping %>% mutate(chr=as.character(chr)), by="chr", relationship="many-to-many") %>%
#   filter(loci >= start, loci <= end)

var <- var.dat.filt
print(head(var))

input.dat <- var %>% unite(mutation, chr, loci, sep="_") %>%
     select(segment, mutation, cell, varReads = var, totReads = total)
# snvs <- input.dat %>% select(segment, mutation) %>% distinct()


snvs <- input.dat %>% group_by(segment, mutation) %>% 
  summarize(var= sum(varReads), total=sum(totReads)) %>% 
  filter(var > minvar) %>% 
  select(segment, mutation)

m <- nrow(snvs)
snvs.per.seg <- snvs %>% group_by(segment) %>% count()



seg.dat <- full_join(snvs.per.seg, nstates) %>% rename(m=n) %>%
replace_na(list(m=0))



med_snvs <- median(seg.dat$m)
selected_segs <-seg.dat %>% arrange(segment)
#  slice_sample(n=nsegs, weight_by=m) %>% 

p1 <- ggplot(selected_segs %>% pivot_longer(c("m", "nstates")), aes(x=name, y=value)) + facet_wrap(~name, scales="free") +
geom_boxplot() + xlab(selected_sample) + ylab("count")

ggsave(file.path(figpath, "selected_seg_dist.pdf"), plot=p1)

p2 <- ggplot(selected_segs, aes(x=factor(nstates), y=m)) +  geom_boxplot(aes(x=factor(nstates))) + geom_point() + xlab("number of cn states")
ggsave(file.path(figpath, "snvs_by_num_cn_states.pdf"), plot=p2)
  




input.dat <- var %>% unite(mutation, chr, loci, sep="_") %>%
     select(segment, mutation, cell, varReads = var, totReads = total) %>% 
     inner_join(snvs) #filter(segment %in% selected_segs$segment)
write_csv(input.dat, var_reads_out)

write_csv(selected_segs, sample_segs_out)

copy_prof <- select(cp, segment, cell, copiesX = x, copiesY=y) %>% 
      filter(segment %in% selected_segs$segment)
write_csv(copy_prof, deep_copy_out)


