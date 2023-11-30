source("src/helper_functions.R")

# Principal components on cov and corr --------------------
pc_viz <- function(pca, type = "?", pt_color = NULL){
  
  var_exp <- ((pca$sdev^2)/sum(pca$sdev^2))
  print(paste0("variance explained, pc1: ", var_exp[1], "\n pc2:", var_exp[2]))
  
  var_exp1 <- round(var_exp[1]*100, 2)
  var_exp2 <- round(var_exp[2]*100, 2)
  
  plot(pca$x[,1], pca$x[,2], type = "n", bty = "l",
       #pch = substr(Z$treat, 1, 1),
       #pch = 19,
       #col = factor(Z$treat),
       #cex = cexx,
       xlab = paste0("PC1 (", var_exp1, " %)"),
       ylab = paste0("PC2 (", var_exp2, " %)")
  )
  mtext(type, 3, line = 0.5)
  
  # Visualize founder and descendants --------------------
  
  found_id %>% walk( ~ {
    
    filter(Z, landscape == .x) %>%
      arrange(gen) %>%
      select(POP) %>%
      pull(POP)
    
    rep_index <- filter(Z, landscape == .x) %>%
      arrange(gen) %>%
      select(POP) %>%
      pull(POP)
    
    if(length(rep_index) > 1){
      arrows(x0 = pca$x[rep_index, ord][1,1], 
             y0 = pca$x[rep_index, ord][1,2], 
             x1 = pca$x[rep_index, ord][2,1], 
             y1 = pca$x[rep_index, ord][2,2],
             col = alph_blk, length = arr_len, lwd = 1.5
      )
    }
    
    if(length(rep_index) > 2){
      arrows(x0 = pca$x[rep_index, ord][1,1], 
             y0 = pca$x[rep_index, ord][1,2], 
             x1 = pca$x[rep_index, ord][3,1], 
             y1 = pca$x[rep_index, ord][3,2],
             col = alph_blk, length = arr_len, lwd = 1.5
      )
    }
  })
  
  
  if(is.null(pt_color)){
    points(pca$x[,1], pca$x[,2],
           pch = Z$pch,
           col = factor(Z$treat))
  }else{
    points(pca$x[,1], pca$x[,2],
           pch = Z$pch,
           col = pt_color)
  }  
}



set.seed(3456)

# Load data --------------------
pch_match <- tibble(
  treat = c("core","edge","founder", "shuffled"),
  pch = c(20, 17, 18, 15)
)

Z <- load_z() %>% 
  full_join(., pch_match, by = "treat") %>% 
  mutate(temporal_block = case_when(
    as.numeric(landscape) %in% 21:40 ~ 1,
    as.numeric(landscape) %in% 41:60 ~ 2,
    as.numeric(landscape) %in% 61:80 ~ 3
  )
  )

#exploring PCs of the covariance matrix
#population genetic covariance matrix
cov_mat <- read.table("data/baypass2/core_model_mat_omega.out")

cvmat <- cov2cor(as.matrix(cov_mat))
diag(cvmat) <- NA
cvmat %>%
  data.frame() %>% 
  mutate(treat = Z$treat) %>% 
  pivot_longer(-treat) %>% 
  ggplot(aes(treat, value)) +
  geom_sina()


diag(as.matrix(cov_mat)) %>% 
  tibble(var = .) %>% 
  mutate(treat = Z$treat) %>% 
  group_by(treat) %>% 
  summarise(mean(var, na.rm = TRUE), sd(var, na.rm = TRUE))


svdx <- svd(cov2cor(as.matrix(cov_mat)))
plot(svdx$u[,1], svdx$u[,2], col = factor(Z$pch))

pca_corr <- prcomp(cov2cor(as.matrix(cov_mat)), scale. = T, center = T)
pca_cov <- prcomp(cov_mat, center = T, scale. = T)
plot(pca_corr$x[,1], pca_corr$x[,2], col = factor(Z$pch))
plot(pca_cov$x[,1], pca_cov$x[,2], col = factor(Z$pch))
plot(svdx$u[,1], pca_cov$x[,1])
plot(svdx$u[,2], pca_corr$x[,2])


# setup --------------------

#grab and sort the indices give the landscape replicate (pull converts to vector)
found_id <- Z %>% 
  filter(gen == 0) %>%
  pull(landscape) %>% 
  as.numeric() %>% 
  sort() %>%
  unique()


#helper variables for viz
alph_blk <- alpha("black", 0.1)
arr_len <- 0
ord <- 1:2
cexx <- 1.5
palette(my_pall)

#PCA of allele frequencies
allele_df_full <- data.table::fread("data/baypass2/core_model_summary_yij_pij.out")
outlier_df <- read_csv("data_out/outlier_df.csv")
#exclude outliers
allele_df <- filter(allele_df_full, !MRK %in% outlier_df$MRK)
rm(allele_df_full)
gc()

#grab random 20K sites
mrks <- sort(sample(allele_df$MRK, 20000))
allele_sample <- allele_df %>% 
  filter(MRK %in% mrks) %>% 
  select(POP, MRK, M_P) %>% 
  pivot_wider(id_cols = POP, names_from = MRK, values_from = M_P)

allele_join <- 
  full_join(
    allele_sample,
    Z,
    by = "POP"
  )
allele_mat <- allele_sample %>% select(-POP)


PC <- prcomp(allele_mat, center = T, scale. = T)
pc_viz(PC, type = "PCA of allele frequencies")


pc_df <- bind_cols(POP = Z$POP, PC$x)
write_csv(pc_df, "data_out/core_freq_PCA.csv")

full_pcdf <- full_join(pc_df, Z, by = "POP")

#quantify distance of populations from their ancestral pop
PC_dists <- 
  found_id %>% 
  map_df(~{
    pc_c <- filter(full_pcdf, landscape == .x) %>% 
      select(treat, starts_with("PC"), -pch)
    pc_F <- filter(pc_c, treat == "founder")
    pc_D <- filter(pc_c, treat != "founder")
    
    popid <- t(combn(pc_c$treat, m = 2))
    colnames(popid) <- c("pop1", "pop2")
    
    tibble(
      tibble(data.frame(popid)),
      distance = as.vector(dist(pc_c[,-1]))
    )
  })


PC_dists %>% 
  filter(pop1 == "founder") %>% 
  group_by(pop1, pop2) %>% 
  summarise(md = mean(distance)) %>% 
  arrange(desc(md)) %>% 
  mutate(mx_md = max(md)/md) %>% 
  pull(mx_md) %>% 
  mean()
  

PC_dists <- 
  found_id %>% 
  map_df(~{
    pc_c <- filter(full_pcdf, landscape == .x) %>% 
      select(treat, starts_with("PC"), -pch)
    pc_F <- filter(pc_c, treat == "founder")
    pc_D <- filter(pc_c, treat != "founder")
    
    popid <- t(combn(pc_c$treat, m = 2))
    colnames(popid) <- c("pop1", "pop2")
    
    tibble(
      tibble(data.frame(popid)),
      distance = as.vector(dist(pc_c[,-1]))
    )
  })



filter(full_pcdf, treat == "core") %>% 
  group_by(temporal_block) %>% 
  group_modify(~{
    distmat <- select(., starts_with("PC")) %>% 
      dist() %>%
      as.matrix()
    distmat[upper.tri(distmat, diag = TRUE)] <- NA
    
    distmat %>% 
      data.frame() %>%
      mutate(pops = colnames(.)) %>% 
      pivot_longer(cols = -pops, names_to = "pop", values_to = "distance") %>%
      mutate(treat = "temp") %>% 
      select(distance, treat) %>% 
      drop_na() %>% 
      arrange(distance)

  })
  


pc_rep_dists <- c("edge", "core", "shuffled", "founder") %>% 
  map_df(function(x){

    filter(full_pcdf, treat == x) %>% 
      group_by(temporal_block) %>% 
      group_modify(~{
        distmat <- select(., starts_with("PC")) %>% 
          dist() %>%
          as.matrix()
        distmat[upper.tri(distmat, diag = TRUE)] <- NA
        
        distmat %>% 
          data.frame() %>%
          mutate(pops = colnames(.)) %>% 
          pivot_longer(cols = -pops, names_to = "pop", values_to = "distance") %>%
          mutate(treat = x) %>% 
          select(distance, treat) %>% 
          drop_na() %>% 
          arrange(distance)
        
      })
  }) 


pc_rep_dists %>% 
  group_by(treat) %>% 
  summarise(md = mean(distance)) %>% 
  arrange(desc(md)) %>% 
  mutate(mx_md = max(md)/md) %>% 
  pull(mx_md) %>% 
  mean()


fig1 <- function(){
  layout(mat = matrix(c(1, 1, 2, 3, 4, 4), nrow = 2, byrow = FALSE), widths = c(3, 2, 1))
  par(mar = c(5,4,4,0)*1.1)
  pc_viz(PC, type = "")
  mtext(text = "A", side = 3, line = -3, outer = T, adj = 0.1)
  par(mar = c(5,4,4,0)*1.1)
  par(bty="l")
  boxplot(PC_dists$distance ~ PC_dists$pop2, col = "white", 
          lwd = 1.5, pch = NA, staplewex = 0, lty = 1, axes = T, 
          bty = "l",
          xlab = "",
          ylab = "PC distance from founder"
  )
  points(jitter(as.numeric(factor(PC_dists$pop2)), factor = 1), 
         PC_dists$distance, 
         cex = 0.7, 
         col = my_pall[PC_dists$pop2],
         pch = my_pch[PC_dists$pop2])
  mtext(text = "B", side = 3, line = -3, outer = T, adj = 0.6)
  
  #C
  par(mar = c(5,4,4,0)*1.1)
  
  boxplot(pc_rep_dists$distance ~ pc_rep_dists$treat, col = NA, 
          lwd = 1.5, pch = NA, staplewex = 0, lty = 1, axes = T, 
          bty = "l",
          xlab = "Population",
          ylab = "PC distance between replicates")
  

  points(jitter(as.numeric(factor(pc_rep_dists$treat)), factor = 1), 
         pc_rep_dists$distance, 
         cex = 0.7, 
         col = my_pall[pc_rep_dists$treat],
         pch = my_pch[pc_rep_dists$treat])
 
  boxplot(pc_rep_dists$distance ~ pc_rep_dists$treat, col = NA, add = TRUE,
          lwd = 1.5, pch = NA, staplewex = 0, lty = 1, axes = F)
  
  mtext(text = "C", side = 3, line = -22, outer = T, adj = 0.6)
  
  par(mar = c(0,0,0,0))
  plot(c(0,10), c(0,10), axes = F, type = "n", xlab = "", ylab = "", xlim = c(-0.2,1))
  
  legend(0, 6, c("core", "edge", "founder", "shuffled"),
         pch = pch_match$pch, col = my_pall,
         bty = "n", cex = 1.2)
  
}

fig1()


####################
####################
#### FIGURE 1 ######
####################
####################



pdf("figures/figure1.pdf", width = 8, height = 5)
fig1()
dev.off()

png("figures/figure1.png", width = 8, height = 5, units = "in", res = 200)
fig1()
dev.off()



pdf("figures/pca_temporalblock.pdf", width = 6, height = 5)
layout(mat = matrix(c(1,2), nrow = 1), widths = c(2.5, 0.85))
par(mar = c(5,4,4,0)*1.1)
pc_viz(PC, type = "PCA of allele frequencies", pt_color = Z$temporal_block)
par(mar = c(0,0,0,0))
plot(c(0,10), c(0,10), axes = F, type = "n", xlab = "", ylab = "")
legend(0, 6, c("block 1", "block 2", "block 3"),
       pch = 15, col = my_pall,
       bty = "n", cex = 0.65)
legend(0, 5, c("core", "edge", "founders", "shuffled"),
       pch = pch_match$pch, col = "black",
       bty = "n", cex = 0.65)
dev.off()

pc_df <- bind_cols(POP = Z$POP, PC$x)
write_csv(pc_df, "data/core_freq_PCA.csv")

full_pcdf <- full_join(pc_df, Z, by = "POP") 


pdf("figures/pca_corr.pdf", width = 6, height = 5)
#cairo_pdf(filename = "~/Dropbox/dissertation/dissertation/figs2/pca.pdf", width = 6, height = 5)
#par(mfrow = c(1,2))
layout(mat = matrix(c(1,2), nrow = 1), widths = c(2.5, 0.85))
par(mar = c(5,4,4,0)*1.1)
pc_viz(pca_corr, "")
par(mar = c(0,0,0,0))
plot(c(0,10), c(0,10), axes = F, type = "n", xlab = "", ylab = "")
legend(0, 6, c("core", "edge", "founders", "shuffled"),
       pch = pch_match$pch, col = my_pall,
       bty = "n", cex = 0.65)
#pc_viz(pca_corr, "correlation") 
dev.off()


stats::heatmap(as.matrix(cov_mat), 
               labRow = Z$treat, 
               labCol = NA, 
               Colv = NA, Rowv = NA)


