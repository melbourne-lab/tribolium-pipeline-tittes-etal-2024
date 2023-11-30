library(tidyverse)

source("src/helper_functions.R")

#used to be tribolium_AllPops_autosomes_auxModel_summary_betai.out
bf <- read.table("data/baypass2/aux_model_summary_betai.out", header = T) %>%
  full_join(treat_key, by = "COVARIABLE")


#used to be AllPops_popdps_4-22_meandp4-22_mincnt0_minor3_majfreq-0.998_autosomes.sites
sites_raw <- read.table("data/baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35.txt") %>%
  set_colnames(c("chrom", "pos", "pos2", "ref_alt", "freq")) %>% 
  select(c(chrom, pos, pos2, ref_alt, freq)) %>% 
  mutate(MRK = 1:n())


build_go_df <- function(sites_df, bf_df, bf_cut, chrom_key = chrom_key){
  bf_go <- 
    bf_df %>% 
    filter(BF.dB. >= bf_cut) %>% 
    left_join(., sites_df, by = "MRK") %>% 
    full_join(., chrom_key, by = "chrom")  %>%
    left_join(., vep, by = "MRK") %>% 
    distinct()
  return(bf_go)
}

bf_go <- build_go_df(sites_df = sites_raw, bf_df = bf, bf_cut = 10, chrom_key = chrom_key)

bf_trait <- read.table("data/baypass2/trait_aux_model_summary_betai.out", header = T) %>%
  full_join(trait_key, by = "COVARIABLE")

bf_trait_go <- build_go_df(sites_df = sites_raw, bf_df = bf_trait, bf_cut = 10, chrom_key = chrom_key)

unique(bf_trait_go$trait)

get_go <- function(bf_sites_df, feature, feature_level, file_out){
  feature <- enquo(feature)
  genename_vec <- bf_sites_df %>% filter(!!feature == feature_level) %>% pull(INFO) %>% info_df() %>% map_chr(~ .x$Feature[1] %>% str_replace("_001", "") %>% str_replace("\\.", ""))
  cat(genename_vec, sep = "\n", file = file_out)
}

#fast
get_go(bf_go, feature = treatment, feature_level = "founder", file_out = "data/go_founder.txt")
get_go(bf_go, feature = treatment, feature_level = "core", file_out = "data/go_core.txt")
get_go(bf_go, feature = treatment, feature_level = "edge", file_out = "data/go_edge.txt")
get_go(bf_go, feature = treatment, feature_level = "shuffled", file_out = "data/go_shuffled.txt")


#not using post revisions
#get_go(bf_trait_go, feature = trait, feature_level = "growth_rate", file_out = "data/go_growth_rate.txt")
#get_go(bf_trait_go, feature = trait, feature_level = "proportion_dispersed", file_out = "data/go_proportion_dispersed.txt")
#get_go(bf_trait_go, feature = trait, feature_level = "mean_weight", file_out = "data/go_mean_weight.txt")


#load output from http://geneontology.org/

colz <- c("GO_biological_process_complete", "count", "observed", "expected", "+/-", "fold_enrichment", "p_value", "fdr")
coltypes <- cols(GO_biological_process_complete = col_character() , count = col_character(), observed = col_integer(), expected = col_number(), `+/-` = col_character(), fold_enrichment = col_number(), p_value = col_number(), fdr = col_number())
pop <-  c("founder", "core", "edge", "shuffled")
gos <- tibble(pop = pop, file = str_glue("data/go_{pop}_out.txt"))
go_df <- pmap_df(gos, function(pop, file){
  read_delim(file, delim = "\t", skip = 12, col_names = colz, col_types = coltypes) %>%
    mutate(population = pop)
})


#load output from http://geneontology.org/
#pop <-  c("growth_rate", "proportion_dispersed", "mean_weight")
#gos <- tibble(pop = pop, file = str_glue("data/go_{pop}_out.txt"))
#go_trait_df <- pmap_df(gos, function(pop, file){
#  read_delim(file, delim = "\t", skip = 12, col_names = colz, col_types = coltypes) %>%
#    mutate(population = pop)
#})


#full_go_df <- bind_rows(go_trait_df, go_df) %>% rename(population = "category")
full_go_df <- go_df %>% rename(population = "category")
full_go_df %>% 
  arrange(category, fdr) %>% 
  select(GO_biological_process_complete, category) %>% 
  distinct() %>% 
  group_by(category) %>% 
  dplyr::slice(1:10) %>% View()

full_go_df %>% 
  select(GO_biological_process_complete, category) %>% 
  distinct() %>% 
  group_by(category) %>% 
  summarise(n())

go_summ <- 
  full_go_df %>% 
  group_by(GO_biological_process_complete) %>% 
  summarise(categories = paste(category, collapse = " "), mean_q_value = mean(fdr)) %>% 
  mutate(nc = nchar(categories)) %>% 
  arrange(mean_q_value,nc) %>% 
  select(-nc)

write_tsv(x = go_summ, path = "data/GO_summary.txt")

View(go_summ)


go_summ %>% filter(categories == "founder") 

go_summ %>% filter(categories == "core")

go_summ %>% filter(categories == "edge") 

go_summ %>% filter(categories == "shuffled") 


bf %>% 
  filter(BF.dB. >= 20) %>% 
  left_join(., sites_raw, by = "MRK") %>% 
  full_join(., chrom_key, by = "chrom")  %>%
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
