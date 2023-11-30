source("src/helper_functions.R")
source("src/make_traitdata.R")

baypass_majmin <-  load_sites() %>% 
  mutate(type = ifelse(freq < 0.5, "min", "maj"))

mean(baypass_majmin$freq > 0.5)
hist(baypass_majmin$freq)

cons_pij <- read.table("data/conserved_sites/aux_model_summary_yij_pij_CONSERVED.out", header = T) %>% 
  as_tibble() %>% 
  full_join(., Z, by = "POP") %>% 
  mutate(edge = ifelse(treat == "edge", "edge", "not edge")) %>% 
  left_join(., baypass_majmin, by = "MRK") %>% 
  left_join(., vep, by = "MRK")

#save markers that are both outliers and conserved loci
read_csv("data_out/outlier_df.csv") %>% 
  filter(MRK %in% cons_pij$MRK) %>% 
  write_csv("data_out/outlier_conserved_df.csv")

info_df(cons_pij$INFO[1:3])

outlier_info <- info_df(unique(cons_pij$INFO))

outlier_info %>% 
  map(~{
    .x$Consequence
  }) %>% 
  flatten() %>% 
  unlist() %>% 
  table()

#outlier_df[grep("missense", outlier_df$INFO), ] 

outlier_consequences <- 
  outlier_info %>%  
  map(~ .x$Consequence) %>%  
  flatten() %>%  
  unlist() %>%
  table() %>% 
  prop.table() %>% 
  round(4) * 100

consequences <- 
  names(outlier_consequences) %>% 
  str_remove("\\&intron|_variant|_gene") %>% 
  str_remove("_variant") %>% 
  str_remove("variant") %>% 
  str_replace_all("_", " ") %>%
  str_replace("_", " ") %>%
  str_remove("&intron")

conserved_vep <- 
data_frame(
  category = consequences,
  percent = as.vector(outlier_consequences)
)

print(conserved_vep)
write_csv(conserved_vep, "data_out/conserved_categories.csv")

cons_pij %>% 
  select(MRK, INFO) %>% 
  unique() %>%
  grep(pattern = "miss", x = .$INFO) 

#15 / 1836
length(grep(pattern = "miss", x = vep$INFO)) / nrow(vep)
length(grep(pattern = "intron", x = vep$INFO)) / nrow(vep)

grepl(pattern = "miss", x = vep$INFO)[20]
hist(cons_pij$M_P)
hist(cons_pij$freq)

#Capture maj/min ref/alt orientation of each locus
polar_df <- cons_pij %>% 
  group_by(MRK) %>% 
  summarise(mean_f = mean(M_P)) %>% 
  left_join(., baypass_majmin, by = "MRK") %>%
  ungroup() %>% 
  mutate(polar_freq = ifelse(mean_f < 0.5 & type == "maj" | mean_f >= 0.5 & type == "min", 1 - mean_f, mean_f),
         polar = ifelse(mean_f < 0.5 & type == "maj" | mean_f >= 0.5 & type == "min", FALSE, TRUE)) %>%
  select(MRK, polar_freq, polar) %>% 
  right_join(., cons_pij, by = "MRK") %>% 
  mutate(polar_M_P = ifelse(polar, M_P, 1 - M_P))

hist(polar_df$polar_M_P)
hist(polar_df$freq)

plot(polar_df$freq, polar_df$polar_M_P, pch = ".")
plot(polar_df$freq, polar_df$polar_freq, pch = ".")  
plot(polar_df$polar_freq, polar_df$freq, pch = ".")


temp <- cons_pij %>% 
  group_by(MRK) %>% 
  summarise(mean_f = mean(M_P)) %>% 
  left_join(., baypass_majmin, by = "MRK") %>%
  ungroup()
#baypass chooses allele at random, have to get ref allele according to raw frequencies
plot(temp$freq, temp$mean_f, pch = ".")


diff_df <- polar_df %>% 
  group_by(MRK) %>% 
  select(polar_M_P, landscape, treat) %>% 
  spread(treat, polar_M_P) %>% 
  mutate(
    edge_d =  founder - edge,
    core_d =  founder - core,
    shuff_d = founder - shuffled,
    edge_d_tr =  log(((founder - edge)/(founder*(1-founder)))^2),
    core_d_tr =  log(((founder - core)/(founder*(1-founder)))^2),
    shuff_d_tr = log(((founder - shuffled)/(founder*(1-founder)))^2)
  )

den_shuff <- density(na.omit(diff_df$shuff_d))
den_core <- density(na.omit(diff_df$core_d))
den_edge <- density(na.omit(diff_df$edge_d))

qs <- seq(0,1, length.out = 100)
q_shuff <- quantile(na.omit(diff_df$shuff_d), probs = qs)
q_core <- quantile(na.omit(diff_df$core_d), probs = qs)
q_edge <- quantile(na.omit(diff_df$edge_d), probs = qs)

wds <- 3

#core, edge, founder, shuffled
# "#F1BB7B", "#FD6467", "#5B1A18", "#D67236"
#plot(1:4, pch = 19, col = c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236"))
plot(den_shuff$x, den_shuff$y, type = "l", bty = "l", 
     xlim = range(den_shuff$x, den_core$x, den_edge$x),
     ylim = range(den_shuff$y, den_core$y, den_edge$y),
     ylab = "Density", 
     col = my_pall[4],
     xlab = expression(Delta[p]),
     lwd = wds,
     cex.lab = 1.2)


###################
#### Figure 5 #####
###################

plot(den_shuff$x, den_shuff$y, type = "l", bty = "l", 
     xlim = range(den_shuff$x, den_core$x, den_edge$x),
     ylim = range(den_shuff$y, den_core$y, den_edge$y))
     

fig5 <- function(){
  cx <- 1.65
  layout(mat = matrix(c(1,2,3), nrow = 1), widths = c(2.5, 2.5, 0.85), heights = c(1, 1, 1))
  par(mar = c(5,4,4,0)*1.1)
  plot(den_shuff$x, den_shuff$y, type = "l", bty = "l", 
       xlim = range(den_shuff$x, den_core$x, den_edge$x),
       ylim = range(den_shuff$y, den_core$y, den_edge$y),
       ylab = "Density", 
       col = my_pall[4],
       xlab = expression(Delta[p]),
       lwd = wds,
       cex.lab = cx)
  lines(den_core, col = my_pall[1], lty = 2, lwd = wds)
  lines(den_edge, col = my_pall[2],lty = 1, lwd = wds)
  mtext(text = "A", side = 3, adj = 0, line = 1, cex = 1.2)
  
  par(mar = c(5,5,4,0)*1.1)
  plot(q_shuff, qs, xlim = c(-0.3, 0.3), type = "l", bty = "l", axes = T,
       col = my_pall[4], lwd = wds, 
       xlab =  expression(Delta[p]),
       ylab = "Cumulative  probability",
       cex.lab = cx)
  lines(q_core, qs, 
        col = my_pall[1], lwd = wds,
        lty = 2)
  lines(q_edge, qs, 
        col = my_pall[2], lwd = wds) 
  mtext(text = "B", side = 3, adj = 0, line = 1, cex = 1.2)
  
  par(mar = c(0,0,0,0))
  plot(c(0,10), c(0,10), axes = F, type = "n", xlab = "", ylab = "")
  legend(0, 6, c("shuffled", "core", "edge"),
         lwd = wds, col = my_pall[c(4, 1, 2)],
         bty = "n", cex = cx - 0.15
  )
  
}


pdf("figures/figure5.pdf", height = 5, width = 9)
fig5()
dev.off()

png("figures/figure5.png", units = "in", width = 9, height = 5, res = 200)
fig5()
dev.off()


qq_col <- function(x, n = 100){
  qx <- seq(0, 1, length.out = n)
  qy <- quantile(x, prob = qx, na.rm = T)
  return(cbind(qx, qy))
}

edge_q <- qq_col(diff_df$edge_d)
core_q <- qq_col(diff_df$core_d)
shuff_q <- qq_col(diff_df$shuff_d)

plot(edge_q[,2], core_q[,2], type = "l", 
     xlim = c(-0.8, 0.8), 
     ylim = c(-0.8, 0.8), 
     col = "blue")
lines(edge_q[,2], shuff_q[,2], col = "red")
lines(core_q[,2], shuff_q[,2])
abline(0,1, lty = 2, col = "grey80")


print("KS TESTS FOR CONSERVED SITES")
ks.test(diff_df$core_d, diff_df$edge_d)
ks.test(diff_df$shuff_d, diff_df$edge_d)
ks.test(diff_df$shuff_d, diff_df$core_d)



#######
### TEST AGAIN AT HIGH IMAPACT SITES
#######
length(unique(cons_pij$MRK))

rm(cons_pij)
gc()

effect <- c("MODERATE", "HIGH")
effect <- c("HIGH")
vep_high <- vep %>% 
  filter(grepl(paste(effect, collapse="|"), INFO))

pij_df <- fread("data/baypass2/aux_model_summary_yij_pij.out")

r_mrk <- unique(sort(vep_high$MRK))
# #randomly chosen markers
cons_pij_dfs <- pij_df %>%
  as_tibble() %>%
  mutate(group = round(MRK/1e5, 0)) %>% 
  group_by(group) %>% 
  group_split()

cons_pij <- cons_pij_dfs %>% 
  map_df(~{
    full_join(.x, Z, by = "POP") %>%
      filter(MRK %in% r_mrk) %>% 
      mutate(edge = ifelse(treat == "edge", "edge", "not edge")) %>%
      left_join(., baypass_majmin, by = "MRK") %>% 
      right_join(., vep_high, by = "MRK")  
  })
  
rm(cons_pij_dfs)
gc()

outlier_info <- info_df(unique(cons_pij$INFO))

outlier_info %>% 
  map(~{
    .x$Consequence
  }) %>% 
  flatten() %>% 
  unlist() %>% 
  table()

outlier_consequences <- 
  outlier_info %>%  
  map(~ .x$Consequence) %>%  
  flatten() %>%  
  unlist() %>%
  table() %>% 
  prop.table() %>% 
  round(4) * 100

consequences <- 
  names(outlier_consequences) %>% 
  str_remove("\\&intron|_variant|_gene") %>% 
  str_remove("_variant") %>% 
  str_remove("variant") %>% 
  str_replace_all("_", " ") %>%
  str_replace("_", " ") %>%
  str_remove("&intron")

conserved_vep <- 
  data_frame(
    category = consequences,
    percent = as.vector(outlier_consequences)
  )

print(conserved_vep)
write_csv(conserved_vep, "data_out/highimpact_categories.csv")


baypass_majmin
#Capture maj/min ref/alt orientation of each locus
polar_df <- cons_pij %>% 
  group_by(MRK) %>% 
  summarise(mean_f = mean(M_P, na.rm = TRUE)) %>% 
  left_join(., baypass_majmin, by = "MRK") %>%
  ungroup() %>% 
  mutate(polar_freq = ifelse(mean_f < 0.5 & type == "maj" | mean_f >= 0.5 & type == "min", 1 - mean_f, mean_f),
         polar = ifelse(mean_f < 0.5 & type == "maj" | mean_f >= 0.5 & type == "min", FALSE, TRUE)) %>%
  select(MRK, polar_freq, polar) %>% 
  right_join(., cons_pij, by = "MRK") %>% 
  mutate(polar_M_P = ifelse(polar, M_P, 1 - M_P))

temp <- cons_pij %>% 
  group_by(MRK) %>% 
  summarise(mean_f = mean(M_P)) %>% 
  left_join(., baypass_majmin, by = "MRK") %>%
  ungroup()
#baypass chooses allele at random, have to get ref allele according to raw frequencies

diff_df <- polar_df %>% 
  group_by(MRK) %>% 
  select(MRK, polar_M_P, landscape, treat) %>% 
  drop_na() %>% 
  pivot_wider(id_cols = c(landscape, MRK), names_from = treat, values_from = polar_M_P) %>%
  mutate(
    edge_d =  founder - edge,
    core_d =  founder - core,
    shuff_d = founder - shuffled,
    edge_d_tr =  log(((founder - edge)/(founder*(1-founder)))^2),
    core_d_tr =  log(((founder - core)/(founder*(1-founder)))^2),
    shuff_d_tr = log(((founder - shuffled)/(founder*(1-founder)))^2)
  )


length(unique(diff_df$MRK))

den_shuff <- density(na.omit(diff_df$shuff_d))
den_core <- density(na.omit(diff_df$core_d))
den_edge <- density(na.omit(diff_df$edge_d))

qs <- seq(0,1, length.out = 100)
q_shuff <- quantile(na.omit(diff_df$shuff_d), probs = qs)
q_core <- quantile(na.omit(diff_df$core_d), probs = qs)
q_edge <- quantile(na.omit(diff_df$edge_d), probs = qs)

qq_col <- function(x, n = 100){
  qx <- seq(0, 1, length.out = n)
  qy <- quantile(x, prob = qx, na.rm = T)
  return(cbind(qx, qy))
}

edge_q <- qq_col(diff_df$edge_d)
core_q <- qq_col(diff_df$core_d)
shuff_q <- qq_col(diff_df$shuff_d)

print("KS TESTS FOR VEP HIGH IMPACT")
print(ks.test(diff_df$core_d, diff_df$edge_d))
print(ks.test(diff_df$shuff_d, diff_df$edge_d))
print(ks.test(diff_df$shuff_d, diff_df$core_d))
fig5()
