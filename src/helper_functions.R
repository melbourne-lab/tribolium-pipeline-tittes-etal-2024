library(tidyverse)
library(cowplot)
theme_set(theme_classic())
library(ggridges)
library(ggrepel)
library(colorspace)
library(ggforce)
library(patchwork)
library(magrittr)
library(wesanderson)
library(ape)
library(IRanges)
library(GenomicRanges)
library(seqinr)
library(Biostrings)
library(zeallot)
#library(moments)
#library(rmutil)
library(data.table)
#library(flextable)
#library(officer)

my_pall <- c("core" = "#6DDBC1", "edge" = "#00BAFF", "founder" = "#AA52BA", "shuffled" = "#F59985")
my_pch <- c("core" = 20, "edge" = 17, "founder" = 18, "shuffled" = 15)


#handy functions
process_rc_file <-  function(filepath) {
  con = file(filepath, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    #do stuff here
    print(get_genic(line))
  }
  close(con)
}

#load tidy version of Z.csv
load_z <- function(){
  read_csv("data/Z.csv") %>%
    gather(treat, log, -prefix) %>%
    filter(log == 1) %>%
    separate(prefix, into = c("landscape", "gen", "treat_p", "id"), sep = "_") %>%
    arrange(landscape) %>%
    mutate(POP = 1:n())
}

Z <- load_z()

#load tidy version of the sites info file (I originally generated along side the input file for baypass)
load_sites <- function(){
  read_tsv("data/conserved_sites/vep_baypass_r0.98_d3_L3_M30_q0.99_a35_thin100.sites", col_names = F) %>% 
    set_colnames(c("MRK", "chrom", "pos", "pos2", "ref_alt", "freq"))
}


#helpful dataframe for making joins
chroms <- factor(c("CM000276.3", "CM000277.3", "CM000278.3", 
                   "CM000279.2", "CM000280.3",
                   "CM000281.3", "CM000282.3",
                   "CM000283.3", "CM000284.3",
                   "CM000285.3"))

chrom_key <- data_frame(
  chrom = chroms,
  chrom_num = c("X", 1:(length(chroms)-1))
)

treat_key <- data_frame(
  COVARIABLE = c(1, 2, 3, 4),
  treatment = c("founder", "core", "edge", "shuffled")
)


vep <- read_tsv("data/VEP_allsites.vcf", skip = 4, 
                col_names = c("CHROM", "POS", "ID", "REF", 
                              "ALT","QUAL", "FILTER","INFO")) %>% 
  mutate(MRK = 1:n())

info_names <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|TSL|APPRIS|ENSP|CLIN_SIG|SOMATIC|PHENO" %>%
  str_split("\\|") %>% 
  unlist()

#organize info column into data frame (takes a long time so don't do it on all loci at once)
info_df <- function(loci){
  loci %>% #!!!
    map(function(y){
      str_split(y, ",") %>% 
        unlist() %>% 
        map_df(~{
          str_split(.x, "\\|") %>% 
            unlist() %>% 
            rbind() %>% 
            set_colnames(info_names) %>% 
            as_tibble()
        })
    })
}
#por exemplo
info_df(vep$INFO[1:10])
