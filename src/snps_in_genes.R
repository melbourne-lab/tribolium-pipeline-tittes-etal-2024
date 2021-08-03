#load libraries and several R objects
source("helper_functions.R")

ref <- read.fasta("data/GCF_000002335.3_Tcas5.2_genomic.fna", as.string = F)

get_genic <- function(rc_line){
  
  #rc_line <- "CM000276.3\t10249\tN\t2\tG/T"
  #extract relevant fields
  c(chrom, pos, alleles) %<-% str_split(rc_line, "\t")[[1]][c(1,2,5)]
  c(ref_allele, alt_allele) %<-% str_split(alleles, "/")[[1]]
  
  #convert chromosome id and position to be numeric
  chrom <- anno_key$alt[anno_key$chrom == chrom]
  pos <- as.numeric(pos)
  
  #create genomic ranges object
  testrr <- GRanges(chrom, IRanges(pos, pos))
  
  #see if snp position overlaps with any gff features (gff with CDS only)
  hits <- gff_ranges[gff_ranges %over% testrr]
  
  #If the snp is in a coding region, translate the cds with ref and alt alleles at snp
  if(length(hits)){
    
    #BE CAREFUL WITH ALL THIS 1 STUFF -- NOT HOW I SHOULD HANDLE MULTIPLE HITS
    seq_id <- as.character(hits@seqnames)[1]
    start <- start(hits)[1]
    end <- end(hits)[1]
    phase_pos <- as.numeric(as.character(hits$phase[1])) + 1
    strand <- hits@strand@values[1]
    
    gene_seq <- ref[[seq_id]][start:end] %>%
      paste(collapse = "") %>%
      DNAString() %>% 
      (function(x) if(strand == "-") reverseComplement(x) else x) %>%
      subseq(start = phase_pos)
    
    v1 <- gene_seq %>% translate()
    
    v2 <- replaceLetterAt(gene_seq, at = pos - start + 1, letter = alt) %>%
      translate()  
    
    return(ifelse(v1 == v2, "Syn", "NonSyn"))
    
  } else {
    return("Non-genic")
  }
}


#match references names
anno_key <- read.table("data/annotation_key.txt", header = F) %>%
  set_colnames(c("seqid","alt","chrom"))

#annotation file of CDSs only
gff <- read.gff(file = "data/ref_Tcas5.2_CDS.gff3") %>%
  full_join(., anno_key, by = "seqid") %>%
  filter(type == "CDS")


gff_ranges <- GRanges(gff$seqid, IRanges(gff$start, gff$end), gff$strand)
mcols(gff_ranges)$phase <- gff$phase
mcols(gff_ranges)$chrom <- gff$chrom
mcols(gff_ranges)$attributes <- gff$attributes


testrr <- GRanges(anno_key$alt[anno_key$chrom == "CM000276.3"], IRanges(106174, 106175))

testrr <- GRanges("NC_007418.3", 
  IRanges(4255934, 4255934)
  )

#first actual denovo gene
test_str <- readLines(con = "data/denovo_rc", n = 1)

#run through denovo file to see if any snps are genic
process_rc_file("data/denovo_rc")

test_str <- "CM000276.3\t10251\tN\t2\tG/T"
get_genic(test_str)


