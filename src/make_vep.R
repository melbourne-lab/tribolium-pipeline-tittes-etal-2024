library(tidyverse)
library(magrittr)
#match references names
anno_key <- read.table("data/annotation_key.txt", header = F) %>%
  set_colnames(c("seqid","alt","chrom"))

read_tsv("data/denovo_rc", col_names = F) %>% 
  select(c(1,2, 5)) %>% 
  set_colnames(c("chrom","pos","snp")) %>% 
  mutate(pos2 = pos) %>% 
  left_join(., anno_key, by = "chrom") %>% 
  select(seqid, pos, pos2, snp) %>% 
  mutate(number = as.integer(1)) %>% 
  write.table(x = ., file = "data/denovo_vep.txt", quote = F, col.names = F, row.names = F)

