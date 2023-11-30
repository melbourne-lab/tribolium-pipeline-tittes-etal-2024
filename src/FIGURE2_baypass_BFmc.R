source("src/helper_functions.R")

#used to be tribolium_AllPops_autosomes_auxModel_summary_betai.out
bf <- read.table("data/baypass2/aux_model_summary_betai.out", header = T) %>%
  full_join(treat_key, by = "COVARIABLE")

#sanity check is the bf dataframe 4 x the number of input loci?
#(nrow(filter(bf, COVARIABLE == 1)) * 4) == nrow(bf)

n <- 4
#used to be AllPops_popdps_4-22_meandp4-22_mincnt0_minor3_majfreq-0.998_autosomes.sites
sites <- read.table("data/baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35.txt") %>%
  set_colnames(c("chrom", "pos", "pos2", "ref_alt", "freq")) %>%
  replicate(n, ., simplify = FALSE) %>% #need to repeat rows for bind_cols
  bind_rows()


bf_full <- bind_cols(bf, sites) %>%
  full_join(., chrom_key, by = "chrom") 

treat_bf <- bf_full %>% 
  filter(treatment == "core",
         BF.dB. > -200) %>%
  mutate(cum_pos = cumsum(as.numeric(pos)))

bf_full <- bf_full %>% mutate(treatment = factor(treatment, levels = rev(c("edge", "core", "shuffled", "founder"))))

#######################
### Base R version ####
#######################

bf_counts <- 
  seq(5, max(bf_full$BF.dB.), length.out = 100) %>% 
  map_df(~{
    bf_full %>%
      filter(BF.dB. >= .x) %>% 
      group_by(COVARIABLE) %>%
      summarise(n_loci = n()) %>% 
      mutate(bf_min = .x)
  }) %>%
  full_join(., treat_key, by = "COVARIABLE")

ES_counts <- 
  #seq(5, max(bf_full$BF.dB.), length.out = 100) %>% 
  seq(0.05, max(abs(bf_full$M_Beta)), length.out = 100) %>% 
  map_df(~{
    bf_full %>%
      filter(abs(M_Beta) >= .x) %>% 
      group_by(COVARIABLE) %>%
      summarise(n_loci = n()) %>% 
      mutate(es_min = .x)
  }) %>%
  full_join(., treat_key, by = "COVARIABLE")

fig2 <- function(){
  #get relative chromosome sizes for scaling panels
  frame_size <- bf_full %>% 
    group_by(chrom_num) %>% 
    summarise(size = max(pos)) %>% 
    mutate(size = size/max(size))
  
  plot_mat <- matrix(1:(nrow(frame_size)*4), nrow = 4, byrow = T)
  #shared Y axis label
  plot_mat <- cbind(1, plot_mat+1)
  #side panels
  panel_b <- matrix(rep(max(plot_mat) + 1, 6), nrow = 2)
  panel_c <- matrix(rep(max(plot_mat) + 2, 6), nrow = 2)
  plot_mat <- cbind(plot_mat, rbind(panel_b, panel_c))
  
  layout(plot_mat, widths = c(0.5, frame_size$size))
  
  ################
  ### PANEL A ####
  ################
  panel_cex <- 1.1
  lab_cex <- 1.2
  par(mar = c(0,0,0,0))
  plot(1,1, type = "n", axes = F, xlab = "", ylab = "")
  text(x = 0.8, y = 1, labels = expression(paste("10 ",log[10]," Bayes Factor")), srt = 90, cex = lab_cex)
  Ts <- c("founder", "shuffled", "core", "edge")
  CN <- sort(unique(bf_full$chrom_num))
  
  pwalk(expand_grid(t=Ts, cn = CN), 
  function(t, cn){
    buffer <- 1e4
    bf_plot <- bf_full %>%
      filter(BF.dB. > -10) %>%
      filter(treatment == t) %>%
      filter(chrom_num == cn)
    
      x_rng <- c(-1e6, max(bf_plot$pos)+1e5) 
    
      TOP <- 2
      BOTTOM <- 4
      LEFT <- 0
      if(t == "edge"){
        if(cn == "1"){
          par(mar = c(BOTTOM, LEFT, 0, 0))
        } else{
          par(mar = c(BOTTOM, 0, 0, 0.2))
        }
      } else if(t == "founder"){
        if(cn == "1"){
          par(mar = c(1, LEFT, TOP, 0))
        } else{
          par(mar = c(1, 0, TOP, 0.2))
        }
      } else{
        if(cn == "1"){
          par(mar = c(1, LEFT, 1, 0))
        } else{
          par(mar = c(1, 0, 1, 0.2))
        }
      }
    
    
      with(bf_plot,
           plot(pos, BF.dB., 
                type = "n", 
                axes = F,
                xlab = "", 
                ylab = "",
                bty = "l",
                ylim = c(0, 65),
                xlim = x_rng)
      )
    
    abline(h = 20, lty = 2)
    with(bf_plot, lines(pos, BF.dB., col = "grey"))
    with(filter(bf_plot, BF.dB. >= 20), points(pos, BF.dB., pch = 19, cex = 0.8, col = my_pall[t]))
    
    if(t == "edge") {
      mtext(text = cn, side = 1, cex = 0.8, line = 0.8)
    }
    if(cn == "1"){
      #strip labels for landscapes
      text(5e5, 62, t, cex = lab_cex, col = my_pall[t], font = 2, adj = 0)
      axis(side = 2, 
           cex.lab = lab_cex, 
           cex.axis = lab_cex, 
           at = seq(-10, 70, length.out = 5))
    }
    if(cn == "1" & t == "founder"){
      mtext(text = "A", side = 3, cex = panel_cex, adj = 0, line = 0.2)
    }
    if(cn == "1" & t == "edge"){
      mtext(text = "Chromosome", side = 1, cex = 0.85, adj = -12, line = 2.2)
    }
  })
  
  
  ################
  ### PANEL B ####
  ################
  par(mar = c(4,5.5, 2, 0.5))

  
  with(bf_counts, 
       plot(bf_min, n_loci, type = "n",bty = "l",
            xlab = expression(paste("10 ", log[10], " Bayes factor")),
            ylab = "Number of loci",
            cex.lab = lab_cex,
            cex.axis = lab_cex
            )
  )
  
  mtext(text = "B", side = 3, cex = panel_cex, adj = 0, line = 0.2)
  
  c("founder", "shuffled", "core", "edge") %>% 
    walk(~{
      with(filter(bf_counts, treatment ==.x),
           lines(bf_min, n_loci, col = my_pall[.x], lwd = 3))
    })
  
  rect(xleft = 5.5, ybottom = 4700, xright = 15, ytop = 5300, border = FALSE, col = "white")
  text(7.75, 5000, "founder", col = my_pall["founder"], font = 2, cex = lab_cex)
  rect(xleft = 5.5, ybottom = 3200, xright = 15, ytop = 3800, border = FALSE, col = "white")
  text(7.75, 3500, "shuffled", col = my_pall["shuffled"], font = 2, cex = lab_cex)
  text(6.5, 1900, "core", col = my_pall["core"], font = 2, cex = lab_cex)
  rect(xleft = 4.5, ybottom = 700, xright = 9, ytop = 1100, border = FALSE, col = "white")
  text(6.5, 800, "edge", col = my_pall["edge"], font = 2, cex = lab_cex)
  
  ################
  ### PANEL C ####
  ################
  
  par(mar = c(4.75,5.5,3,0.5))
  
  with(ES_counts, 
       plot(es_min, n_loci, type = "n",bty = "l",
            xlab = expression(paste("Effect size (", beta[i], ")")),
            ylab = "Number of loci",
            cex.lab = lab_cex,
            cex.axis = lab_cex)
       )
  
  mtext(text = "C", side = 3, cex = panel_cex, adj = 0, line = 0.2)
  c("founder", "shuffled", "core", "edge") %>% 
    walk(~{
      with(filter(ES_counts, treatment ==.x),
           lines(es_min, n_loci, col = my_pall[.x], lwd = 3))
    })
}

png("figures/figure2.png", height = 6, width = 10.5, units = "in", res = 200)
fig2()
dev.off()

pdf("figures/figure2.pdf", height = 6, width = 10.5)
fig2()
dev.off()


#############################
#  MAKE OUTLIER DATA FRAME ##
#############################

bfcut <- 20
outlier_df <- bf_full %>%
  filter(BF.dB. >= bfcut) %>% 
  left_join(., vep, by = "MRK")

write_csv(outlier_df, path = "data_out/outlier_df.csv")

write_csv(outlier_df[,!grepl("\\.\\.\\.", colnames(outlier_df))], path = "data_out/suppmat_outlier_df.csv")


outlier_df %>% 
  group_by(treatment) %>% 
  summarize(n()) %>% 
  ungroup()

outlier_info <- info_df(outlier_df$INFO) 

outlier_info %>% 
  map(~{
    .x$Consequence
  }) %>% 
  flatten() %>% 
  unlist() %>% 
  table()
outlier_df[grep("missense", outlier_df$INFO), ] 


outlier_consequences <- 
  info_df(outlier_df$INFO) %>%  
  map(~ .x$Consequence) %>%  
  flatten() %>%  
  unlist() %>%
  table() %>% 
  prop.table() %>% 
  round(4) * 100


info_DF <- info_df(outlier_df$INFO) %>%  
  map_df(~ .x)

##########################
### SUPPLAMENTAL TABLE ###
##########################
write_tsv(info_DF, "data_out/outlier.tsv")
View(info_DF)

filter(info_DF, Consequence == "missense_variant") %>% 
  pull(Gene) %>% 
  unique()

consequences <- 
  names(outlier_consequences) %>% 
  str_remove("\\&intron|_variant|_gene") %>% 
  str_remove("_variant") %>% 
  str_remove("variant") %>% 
  str_replace("_", " ") %>%
  str_replace("_", " ") %>%
  str_remove("&intron")

category_table <- 
  data_frame(
    category = consequences,
    percent = as.vector(outlier_consequences)
  ) 

print(category_table)
write_csv(category_table, "data_out/outlier_categories.csv")


#doc <- body_add_flextable(doc, value = category_table, pos = "after")
#print(doc, "data/tables.docx", pos = "after")

