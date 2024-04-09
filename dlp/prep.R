library(tidyverse)
vtext <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
df <- read_csv("DeepCopyPrediction.csv")
colnames(df) <- c("cell", "chr", "start", "end", "x", "y")
n_distinct((df$cell))

bin_mapping <- df %>% select(chr, start, end) %>% distinct() %>% arrange(chr, start, end) %>%
  mutate(bin = row_number()-1)

df <- inner_join(df, bin_mapping)
cp <- select(df, cell,chr, bin, x, y) %>% unite(col=cn, x, y, sep ="|", remove=FALSE)



num_states_df <- cp %>% group_by(bin, cn) %>%
  summarize(num_cells =n())




compute_dist <- function(bin_ids, pseudo=1e-7){
  bin_data <- filter(cp, bin %in% bin_ids) %>% select(-bin)
  states <- unique(cp$cn)
  norm <- n_distinct(bin_data$cell) * length(bin_ids)
  dat <- left_join(data.frame(cn= states), bin_data,  by="cn") %>%
    mutate(count = ifelse(is.na(cell),0, 1 ))

  dist <- dat %>%
      group_by(cn) %>% 
    summarize(counts = sum(count), freq = sum(count)/norm)
  
  dist$bin = bin_ids
  
  return(dist)
  
  
}

segments.df <- read_csv("segmentation.csv") %>% mutate(segment= row_number()-1) %>%
  right_join(bin_mapping %>% rename(start_loci = start, end_loci=end), relationship="many-to-many") %>%
  filter(bin >= start, bin <=end) 

write_csv(segments.df, "segments.bin.mapping.csv")

var.dat <- read.table("variant_data.filt.tsv", header = F, sep="\t",
                      col.names=c("chr", "loci", "cell", "base", "var", "total") )
head(segments.df)
var <- var.dat %>% 
  left_join(segments.df %>% mutate(chr=as.character(chr)), by="chr", relationship="many-to-many") %>%
  filter(loci >= start_loci, loci <= end_loci)



snvs <- var %>% select(chr, loci, segment, start_loci, end_loci, bin) %>% distinct()

ggplot(snvs, aes(x=factor(segment))) + geom_bar() + xlab("segment") + vtext
  

snvs.per.seg <- snvs %>% group_by(segment) %>% count()
ggplot(snvs.per.seg, aes(x=n)) + geom_histogram(binwidth = 25, fill="white", color="black") +
  xlab("# of SNVs per segment") + ylab("number of segments")

ggplot(snvs.per.seg, aes(x=n, y="0")) + geom_boxplot() + 
  theme(axis.text.y = element_blank(),  axis.ticks.y = element_blank()) +  
  xlab("# of SNVs per segment") + ylab("")
  

ggplot(seg.cp, aes(x=factor(segment), y=cell, fill=cn)) + geom_tile() +  
  facet_wrap(~chr, scales="free_x", labeller = "label_both") + guides(fill="none")  + 
  xlab("segment") +  theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank())

#now look at how many cn states per segment with and without imputation

head(seg.cp)


seg.cp <- left_join(segments.df , cp) 
num.cn <- seg.cp %>% group_by(chr, segment) %>% summarize(num_cn_states=n_distinct(cn))

ggplot(num.cn, aes(x=factor(segment), y=num_cn_states)) + geom_bar(stat="identity") +
  vtext + xlab("segment") + ylab("# of unique CN states per segment")

ggplot(num.cn, aes(x=num_cn_states)) + geom_histogram( fill="white", color="black") +
  xlab("# of unique CN states per segment") + ylab("number of segments")

ggplot(num.cn, aes(x=num_cn_states, y="0")) + geom_boxplot() + 
  theme(axis.text.y = element_blank(),  axis.ticks.y = element_blank()) +  
  xlab("# of unique CN states per segment") + ylab("")

seg.cell.impute <- seg.cp %>% group_by(chr, segment, cell, cn) %>% count()
ggplot(seg.cell.impute, aes(x=cell, fill=cn)) + geom_bar(position="fill") +
  guides(fill="none") + theme(axis.text.x  = element_blank(),
                              axis.ticks.x = element_blank()) +
  facet_wrap(~segment, labeller = "label_both")

seg.cell.state <- seg.cell.impute %>% group_by(chr, segment, cell) %>% 
  arrange(desc(n)) %>% filter(row_number()==1) %>% select(-n) %>% rename(cn_impute=cn)
head(seg.cell.state)

head(seg.cp)
seg.cp.imp <- inner_join(seg.cp, seg.cell.state)

num.cn.imp <- seg.cp.imp %>% group_by(chr, segment) %>% summarize(num_cn_states=n_distinct(cn_impute))

ggplot(num.cn.imp, aes(x=factor(segment), y=num_cn_states)) + geom_bar(stat="identity") +
  vtext + xlab("segment") + ylab("# of unique CN states per segment")



ggplot(num.cn.imp, aes(x=num_cn_states, y="0")) + geom_boxplot() + 
  theme(axis.text.y = element_blank(),  axis.ticks.y = element_blank()) +  
  xlab("# of unique CN states per segment") + ylab("")

seg.cell.impute <- seg.cp %>% group_by(chr, segment, cell, cn) %>% count()
ggplot(seg.cell.impute, aes(x=cell, fill=cn)) + geom_bar(position="fill") +
  guides(fill="none") + theme(axis.text.x  = element_blank(),
                              axis.ticks.x = element_blank()) +
  facet_wrap(~segment, labeller = "label_both")

ggplot(seg.cp.imp, aes(x=factor(segment), fill=cn_impute)) + geom_bar(position="fill") +
  xlab("segment") + guides(fill="none") + vtext

head(seg.cp.imp)
mtcars %>%
  count(am, gear) %>% 
  mutate(freq = prop.table(n), .by = am)

seg.states <- seg.cp.imp %>% select(chr, segment, cn_impute, cell ) %>% distinct() %>% 
  count(segment, cn_impute) %>%
  mutate(prop = prop.table(n), .by =segment) %>% mutate(include = prop > 0.05)

seg.state.counts <- seg.states %>% group_by(segment) %>% summarize(num_imp_filt = sum(include))
ggplot(seg.state.counts, aes(x=num_imp_filt, y="0")) + geom_boxplot() + 
  theme(axis.text.y = element_blank(),  axis.ticks.y = element_blank()) +  
  xlab("# of unique CN states per segment after imputation & filtering ( > 0.05)") + ylab("")
# inner_loop <- function(start){
#   res <- data.frame()
#   for(j in (start+1):nbins-1){
#     arc1 <- start:j
#     #bissect the data from 0 to j and then j+1 to nbins
#     bins1 <- compute_dist(start:j) %>% rename(freq1 = freq)
#     bins2  <- compute_dist((j+1): nbins) %>% rename(freq2 = freq)
#     freq.df <- inner_join(bins1, bins2, by="cn")
#     
#     
#     div <- JSD(freq.df$freq1, freq.df$freq2)
#     res <- res %>% rbind(data.frame(start=i, end=j, jsd=div))
#     
#   }
#   
# }


i <- 0


ggplot(res %>%filter(end-start > 10, end >40, end < 60), aes(x=end, y=jsd)) + geom_point() 
filter(res, end==46)
bins <- unique(cp$bin)

dist.df <- bind_rows(lapply(bins, compute_dist)) 
count.dat <- dist.df %>% select(-freq) %>% 
  pivot_wider(names_from = bin, names_prefix = "bin", values_from = counts)

proc.dat <- dist.df %>% select(-freq) %>% inner_join(bin_mapping %>% select(bin, chr))

chr1.dat <- proc.dat %>%   filter(chr==1) %>% select(-chr) 
write_csv(select(proc.dat, chr, bin, cn, counts), "counts.csv")

chr1.dat.wide <- chr1.dat %>%   
  pivot_wider(names_from = bin, names_prefix = "bin", values_from = counts)
write_csv(chr1.dat.wide, "chr1.csv")
write_csv(count.dat, "counts_by_states.csv")
# dist.wide <- dist.df %>% pivot_wider(id_cols=cn, names_from = bin, names_prefix = "bin", values_from =freq )

# ent_vals <- numeric(ncol(dist.wide)-1)
# for(i in 2:(ncol(dist.wide)-1)){
#   p <- dist.wide[,i]
#   q <- dist.wide[,i+1]
#   ent_vals[i] <- kl_div(p,q)
# }
# summary(ent.df$kldiv)
ggplot(cp %>% filter(chr==1)  %>% inner_join(bin_mapping) ,
       aes(x=bin, y=cell, fill=cn)) + geom_tile()  +
  geom_vline(aes(xintercept=7)) +  geom_vline(aes(xintercept=14))  +
  geom_vline(aes(xintercept=23))# +    geom_vline(aes(xintercept=22))
#[(0, 5), (6, 6), (7, 15), (16, 22), (23, 23), (24, 24), (25, 28)]


# kl_div <- function(p, q){
#   return(sum(p*log2(p/q)))
# }
# 
# JSD <- function(p,q){
#   n <- 0.5 * (p + q)
#   JS <- 0.5 * (sum(p * log2(p / n)) + sum(q * log2(q / n)))
#   return(JS)
# }
# ggplot(cp %>% filter(bin >=100, bin <=125), 
#        aes(x=factor(bin), y=cell, fill=cn)) + geom_tile()


