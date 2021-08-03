library(tidyverse)

source("helper_functions.R")

#used to be tribolium_AllPops_autosomes_auxModel_summary_betai.out
bf <- read.table("data/baypass2/aux_model_summary_betai.out", header = T) %>%
  full_join(treat_key, by = "COVARIABLE")

#sanity check is the bf dataframe 4 x the number of input loci?
#(nrow(filter(bf, COVARIABLE == 1)) * 4) == nrow(bf)
n <- 4
#used to be AllPops_popdps_4-22_meandp4-22_mincnt0_minor3_majfreq-0.998_autosomes.sites
sites <- read.table("data/baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.txt") %>%
  set_colnames(c("chrom", "pos", "pos2", "ref_alt", "freq")) %>%
  replicate(n, ., simplify = FALSE) %>% #need to repeat rows for bind_cols
  bind_rows()

bfcut <- 10

bf_go <- bind_cols(bf, sites) %>%
  full_join(., chrom_key, by = "chrom")  %>%
  filter(BF.dB. >= bfcut) %>% 
  left_join(., vep, by = "MRK") %>% 
  distinct()

#fast
cat(bf_go %>% filter(treatment == "founder") %>% pull(INFO) %>% info_df() %>% map_chr(~ .x$Feature[1] %>% str_replace("_001", "")), sep = "\n", file = "data/go_founder.txt")
cat(bf_go %>% filter(treatment == "core") %>% pull(INFO) %>% info_df() %>% map_chr(~ .x$Feature[1] %>% str_replace("_001", "")), sep = "\n", file = "data/go_core.txt")
cat(bf_go %>% filter(treatment == "edge") %>% pull(INFO) %>% info_df() %>% map_chr(~ .x$Feature[1] %>% str_replace("_001", "")), sep = "\n", file = "data/go_edge.txt")
cat(bf_go %>% filter(treatment == "shuffled") %>% pull(INFO) %>% info_df() %>% map_chr(~ .x$Feature[1] %>% str_replace("_001", "")), sep = "\n", file = "data/go_shuffled.txt")

#load output from http://geneontology.org/
colz <- c("GO_biological_process_complete", "count", "observed", "expected", "+/-", "fold_enrichment", "p_value", "fdr")
pop <-  c("founder", "core", "edge", "shuffled")
gos <- tibble(pop = pop, file = str_glue("data/go_{pop}_out.txt"))
go_df <- pmap_df(gos, function(pop, file){
read_delim(file, delim = "\t", skip = 12, col_names = colz) %>%
           mutate(population = pop)
})

go_df %>% 
  arrange(population, fdr) %>% 
  group_by(population) %>% 
  dplyr::slice(1:5)

go_df %>% 
  group_by(population) %>% 
  summarise(n())


go_summ <- 
  go_df %>% 
  group_by(GO_biological_process_complete) %>% 
  summarise(populations = paste(population, collapse = " "), mean_q_value = mean(fdr)) %>% 
  mutate(nc = nchar(populations)) %>% 
  arrange(nc) %>% 
  select(-nc)

write_tsv(x = go_summ, path = "data/GO_summary.txt")
  
    
go_summ %>% 
  filter(populations %in% pop) %>%
  group_by(populations) %>% 
  summarise(n())

go_summ %>% 
  filter(populations %in% pop) %>% 
  View()

bind_cols(bf, sites) %>%
  full_join(., chrom_key, by = "chrom")  %>%
  filter(BF.dB. >= 20) %>% 
  left_join(., vep, by = "MRK") %>% 
  distinct() %>% 
  filter(grepl(pattern = "missense", x = INFO)) %>% 
  pull(INFO) %>% 
  info_df() %>% 
  map_chr(~ .x$Feature[1] %>% str_replace("_001", "")) %>% 
  unique()



#https://metazoa.ensembl.org/Tribolium_castaneum/Gene/Literature?db=core;g=TC005411;r=LG8:13518586-13525646;t=TC005411_001
#http://metazoa.ensembl.org/Tribolium_castaneum/Gene/Literature?db=core;g=TC001916;r=LG9:1624726-1732489;t=TC001916_001
#http://metazoa.ensembl.org/Tribolium_castaneum/Gene/Literature?db=core;g=TC016314;r=LG5:983942-992524;t=TC016314_001

