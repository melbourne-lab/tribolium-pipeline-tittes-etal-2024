#COMBINE TRAIT DATA WITH ALLELE FREQUENCY DATA
#--------
source("src/helper_functions.R")
source("src/make_traitdata.R")
library(glmnet)

set.seed(456745)
Z <- load_z()
vep_c <- vep %>% mutate(MRK = as.character(MRK))
outlier_df <- read_csv("data_out/outlier_df.csv")
pc_df <- read_csv("data_out/core_freq_PCA.csv") %>% 
  select(1:9)


#sites <- read.table("data/baypass2/vep_baypass_r0.98_d3_L3_M30_q0.99_a35.txt") %>%
#  ungroup() %>% 
#  set_colnames(c("chrom", "pos", "pos2", "ref_alt", "freq")) %>%
#  mutate(MRK = 1:n())

all_freqs <- fread("data/baypass2/aux_model_summary_yij_pij.out")


gwa_df <- read.table("data/baypass2/aux_model_summary_betai.out", header = T) %>% 
  filter(BF.dB. >= 30) %>% 
  right_join(all_freqs, ., by = "MRK") %>% 
  distinct() %>% 
  drop_na()

# mx_patch <- trait_df %>% 
#   group_by(POP) %>% 
#   #summarise(mx_patch = mean(total_weight_census)) %>% 
#   summarise(mx_patch = mean(Patch)) %>% 
#   arrange(POP) %>% 
#   ungroup() %>% 
#   full_join(., Z, by = "POP") %>% 
#   drop_na() %>% 
#   select(POP, mx_patch)

trait_df %>%
  ggplot(aes(treat, prop_disp)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0) 

trait_df %>% 
  ggplot(aes(treat, prop_10)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0) 

trait_df %>% 
  ggplot(aes(treat, prop_40)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0) 

trait_df %>% 
  ggplot(aes(treat, mean_gr)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0)

trait_df %>% 
  ggplot(aes(treat, mean_weight)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0)


trait_lasso <- function(in_df, trait_var, nfolds, hist_resid = FALSE, return_mod = FALSE, ...){
  #in_df <- trait_df; trait_var <- "prop_disp"; maxit = 10e8; permute = FALSE; hist_resid = TRUE
  
  trait_var <- enquo(trait_var)
  
  mx_patch_df <- 
    in_df %>% 
    ungroup() %>%  
    select( !!trait_var, POP, treat) %>% 
    dplyr::rename(mx_patch = !!trait_var) %>% 
    drop_na()

  #fubar here
  gwa_train <-
    gwa_df %>% 
    select(MRK, POP, M_P) %>%
    distinct() %>% 
    pivot_wider(names_from = c(MRK), values_from = M_P) %>% 
    full_join(., mx_patch_df, by = "POP", na_matches = "never") %>% 
    left_join(., pc_df, by = "POP", na_matches = "never") %>% 
    ungroup() %>% 
    drop_na()
  

  X <- gwa_train %>% 
    select(-POP) %>% 
    model.matrix(mx_patch ~ . -1, data = .) %>% 
    scale()
  

  if(missing(nfolds)){
    nfolds <- nrow(gwa_train)
  }
  
  patch_cv <- cv.glmnet(X, gwa_train$mx_patch, nfolds = nfolds, grouped = F, ...)

  out_df_short <- outlier_df %>% 
    mutate(MRK = as.character(MRK)) %>% 
    select(MRK, treatment)
  
  coef_df <- coef(patch_cv) %>% 
    as.matrix() %>% 
    t() %>% 
    as_tibble() %>% 
    gather("MRK", "coef") %>%
    filter(abs(coef) > 0.0) %>%
    mutate(MRK = sub(pattern = "\\`", replacement = "", x = MRK)) %>% 
    mutate(MRK = sub(pattern = "\\`", replacement = "", x = MRK)) %>% 
    filter(!MRK %in% c("(Intercept)") ) %>%
    left_join(., vep_c, by = "MRK", na_matches = "never") %>% 
    left_join(., out_df_short, by = "MRK", na_matches = "never")
  
  #plot(gwa_train$mx_patch, predict(patch_cv, newx = X))
  #points(gwa_test$mx_patch, predict(patch_cv, newx = X_test), col = "red", pch = 19)
  #abline(0,1)
  
  #cor(gwa_train$mx_patch, predict(patch_cv, newx = X)) ^ 2
  #r_squared <- summary(lm(predict(patch_cv, newx = X_train) ~ gwa_test$mx_patch))$adj.r.squared
  #print(r_squared)
  
  #pred_lasso <- lm(predict(patch_cv, newx = X) ~ gwa_train$mx_patch)
  #r_squared <- summary(pred_lasso)$adj.r.squared
  
  #conservative estimate of R_squared based on deviance
  
  r_squared <- patch_cv$glmnet.fit$dev.ratio[which(patch_cv$glmnet.fit$lambda == patch_cv$lambda.1se)]
  
  if(hist_resid){
    pred <- predict(patch_cv, newx = X)
    resid <- pred - gwa_train$mx_patch
    hist(resid)
    
    plot(pred, gwa_train$mx_patch)
    abline(0,1)
    
  }
  
  
  if(return_mod){
    return(patch_cv)
  } else{
    return(coef_df %>% mutate(r_squared = rep(r_squared, n())))  
  }
  
}




lasso_n <- function(in_df2, trait_var, n_runs = 100, n_permutes = 1, prop_in = 0.95, nfolds = 10, permute = FALSE, permute_within = FALSE){

  trait_var <- quo_name(enquo(trait_var))
  
  run_n <- seq_len(n_permutes) %>% 
    map_df(function(p){
      if(permute){
        if(permute_within){
          perm_df <- in_df2 %>% group_by(treat) %>% mutate_at(vars(!!!trait_var), function(x) sample(x, n(), replace = FALSE)) %>% ungroup()
        } else {
          perm_df <- in_df2 %>%  mutate_at(vars(!!!trait_var), function(x) sample(x, n(), replace = FALSE)) %>% ungroup()
        }
      } else {
        perm_df <- in_df2
      }
        seq_len(n_runs) %>% 
        map_df(~ {
          trait_lasso(in_df = perm_df, trait_var = trait_var, nfolds = nfolds) %>% 
            mutate(permute_index = p)
        })
    }) 
  
  mark_in <- names(table(run_n$MRK)[table(run_n$MRK) >= n_runs*prop_in])
  run_n %>% 
    filter(MRK %in% mark_in) %>% 
    group_by(MRK, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, treatment) %>%
    summarise(mean_coef = mean(coef),
              mean_r_squared = mean(r_squared),
              permute_index = paste(unique(permute_index), collapse = " "),
              q_low = quantile(coef, 0.25),
              q_high = quantile(coef, 0.75)) %>% 
    ungroup()
}
  
  

set.seed(45456)
trait_lasso(trait_df, trait_var = prop_disp, maxit = 10e8, permute = FALSE, hist_resid = TRUE)
trait_lasso(trait_df, trait_var = mean_gr, maxit = 10e8, permute = FALSE, hist_resid = TRUE)
trait_lasso(trait_df, trait_var = mean_weight, maxit = 10e8, permute = FALSE, hist_resid = TRUE) %>% View()

OUT <- trait_lasso(trait_df, trait_var = mean_gr, maxit = 10e8, hist_resid = TRUE, return_mod = TRUE)

# plot(OUT)
# points(log(0.01766), 0.01624)
# points(log(0.04275), 0.01936, pch = 19)
# coef(OUT) %>% 
#   as.matrix() %>% 
#   t() %>% 
#   as_tibble() %>% 
#   gather("MRK", "coef") %>%
#   filter(coef > 0)
# plot(log(OUT$glmnet.fit$lambda), OUT$glmnet.fit$dev.ratio)


#lasso_disp40 <- lasso_n(trait_var = prop_40, prop_in = 0.95) 
#lasso_disp10 <- lasso_n(trait_var = prop_10, prop_in = 0.95)

lasso_disp <- lasso_n(trait_df, trait_var = prop_disp, n_runs = 100, prop_in = 0.95, permute = FALSE)
lasso_disp
lasso_gr <- lasso_n(trait_df, trait_var = mean_gr, n_runs = 100, prop_in = 0.95, permute = FALSE)
lasso_gr
lasso_weight <- lasso_n(trait_df, trait_var = mean_weight, n_runs = 100, prop_in = 0.95, permute = FALSE)
lasso_weight


#SLOW!!!
lasso_gr_permuted_within <- lasso_n(trait_df, trait_var = mean_gr, n_runs = 100,  n_permutes = 100, prop_in = 0.95, permute = TRUE, permute_within = TRUE)
write_csv(lasso_gr_permuted_within, "data_out/lasso_gr_permuted_within.csv")
lasso_weight_permuted_within <- lasso_n(trait_df, trait_var = mean_weight,  n_permutes = 100, prop_in = 0.95, permute = TRUE, permute_within = TRUE)
write_csv(lasso_weight_permuted_within, "data_out/lasso_weight_permuted_within.csv")
lasso_gr_permuted_across <- lasso_n(trait_df, trait_var = mean_gr, n_runs = 100, n_permutes = 100, prop_in = 0.95, permute = TRUE)
write_csv(lasso_gr_permuted_across, "data_out/lasso_gr_permuted_across.csv")
lasso_weight_permuted_across <- lasso_n(trait_df, trait_var = mean_weight,  n_permutes = 100, prop_in = 0.95, permute = TRUE)
write_csv(lasso_weight_permuted_across, "data_out/lasso_weight_permuted_across.csv")

##########
#PERMUTATIONS
##########

lasso_gr_permuted_within <- read_csv("data_out/lasso_gr_permuted_within.csv")
lasso_weight_permuted_within <- read_csv("data_out/lasso_weight_permuted_within.csv")
lasso_gr_permuted_across <- read_csv("data_out/lasso_gr_permuted_across.csv")
lasso_weight_permuted_across <- read_csv("data_out/lasso_weight_permuted_across.csv")



get_nulldist <- function(permuted){
  1:100 %>% map_dbl(~{
    p_df <- filter(permuted, grepl(.x, permute_index))
    if(nrow(p_df) == 0){
      i <- 0
    } else{
      i <- quantile(p_df$mean_coef, 1)
    }
    i
  })
}


#simpler helper functions to asses permutations
plot_coef <- function(permuted, empirical, q = 0.975, ...){
  p_dist <- get_nulldist(permuted)
  hist(p_dist, xlim = c(-max(c(p_dist, empirical))*1.2, max(c(p_dist, empirical))*1.2), ...)
  abline(v = empirical, col = "blue")
}

         
idx_table <- function(idx_column){
  paste(idx_column, collapse  = " ") %>% 
    str_split(" ", simplify = T) %>% 
    as.numeric() %>% 
    table()
}

prop_gte <- function(permuted, empirical){
  p_dist <- get_nulldist(permuted)
  empirical %>% 
    map_dbl(~ mean(abs(.x) >= abs(p_dist)))
}

plot_coef(lasso_weight_permuted_within, lasso_weight$mean_coef)
plot_coef(lasso_weight_permuted_across, lasso_weight$mean_coef)
plot_coef(lasso_gr_permuted_within, lasso_gr$mean_coef)
plot_coef(lasso_gr_permuted_across, lasso_gr$mean_coef)

prop_gte(lasso_weight_permuted_across, lasso_weight$mean_coef)
prop_gte(lasso_gr_permuted_across, lasso_gr$mean_coef)

prop_gte(lasso_weight_permuted_within, lasso_weight$mean_coef)
prop_gte(lasso_gr_permuted_within, lasso_gr$mean_coef)

(props <- c(prop_gte(lasso_gr_permuted_within, lasso_gr$mean_coef),
           prop_gte(lasso_weight_permuted_within, lasso_weight$mean_coef))
)

idx_table(lasso_weight_permuted_across$permute_index) %>% length()
idx_table(lasso_weight_permuted_within$permute_index) %>% length()

idx_table(lasso_gr_permuted_across$permute_index) %>% length()
idx_table(lasso_gr_permuted_within$permute_index) %>% length()

idx_table(lasso_gr_permuted_across$permute_index) %>% quantile(c(0.5, 0.975))
idx_table(lasso_gr_permuted_within$permute_index) %>% quantile(c(0.5, 0.975))

idx_table(lasso_weight_permuted_across$permute_index) %>% quantile(c(0.5, 0.975))
idx_table(lasso_weight_permuted_within$permute_index) %>% quantile(c(0.5, 0.975))

table(lasso_weight_permuted_across$treatment)
table(lasso_gr_permuted_across$treatment)

table(lasso_weight_permuted_within$treatment)
table(lasso_gr_permuted_within$treatment)


wt_df <- 
  lasso_weight %>% 
  ungroup() %>% 
  mutate(trait = rep("weight", n())) 

gr_df <- 
  lasso_gr %>% 
  ungroup() %>% 
  mutate(trait = rep("growth rate", n()))

out_trait_df <- 
  bind_rows(
    gr_df,
    wt_df
  ) %>% 
  select(CHROM, POS, REF, ALT, INFO, treatment, trait, mean_coef)

outlier_table <- 
  out_trait_df %>% 
  mutate("gene" = info_df(out_trait_df$INFO) %>% map_chr(~ .x$Feature[1] %>% str_replace("_001", "")),
         "con" = info_df(out_trait_df$INFO) %>% map_chr(~ .x$Consequence[1] %>% str_replace("_", " ") %>% str_replace("_", " ") %>% str_replace("variant", "") %>% str_replace("gene", ""))) %>% 
  left_join(., chrom_key, by = c("CHROM" = "chrom")) %>% 
  select(chrom_num, POS, REF, ALT, treatment, trait, mean_coef, gene, con) %>% 
  set_colnames(c("Chrom.", "Position", "Ref.", "Alt.", "Treatment", "Trait", "Coeff.", "Gene", "category")) %>%
  mutate(props_gt = props)



###############
### TABLE 3 ###
###############


print(outlier_table)

write_csv(outlier_table, "data_out/trait_outlier_table.csv")
mrk <- as.integer(wt_df$MRK)

mrk_df <- all_freqs %>% 
  filter(MRK == mrk) %>% 
  full_join(., Z, by = "POP") %>% 
  full_join(., trait_df, by = "POP") %>% 
  drop_na()

mrk_df %>%   
  ggplot(aes(treat.x, M_P)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0)


#Predictor residual plotting function
make_avdf <- function(assoc_df, av_obj, org_df){
  assoc_df <- mutate(assoc_df, as.numeric(MRK))
  MRK <- assoc_df$MRK
  seq_along(MRK) %>% 
    map_df(~{
      data.frame(av_obj[[.x]]) %>%
        set_colnames(c("allele_residual", "trait_residual")) %>% 
        mutate(MRK = MRK[.x]) %>% 
        bind_cols(., org_df)
    }) %>% 
    full_join(., assoc_df) %>% 
    mutate(mrk_treat = paste0(treatment, " (", CHROM, " : ", POS, ")"))
}


#Predictor residual plot for weight outliers ----
weight_sites_df <- 
  all_freqs %>% 
  filter(MRK %in% wt_df$MRK) %>% 
  full_join(., trait_df, by = "POP")

wt_model_df <-   
  weight_sites_df %>% 
  select(MRK, POP, M_P, mean_weight) %>%
  drop_na() %>% 
  arrange(POP) %>% 
  spread(MRK, M_P) %>%
  left_join(., Z, by = "POP") 

#added variable plots for weight
wt_model <- lm(mean_weight ~  `147616` + `8103`, data = wt_model_df)
summary(wt_model)
av_wt <- car::avPlots(wt_model)

wt_av_df <-  make_avdf(assoc_df = wt_df, av_obj = av_wt, org_df = wt_model_df) 

wt_fig_out <- wt_av_df %>% 
  ggplot(aes(allele_residual, trait_residual, colour = treat)) + 
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F, colour = "grey50") +
  facet_wrap(~mrk_treat, ncol = 1) +
  guides(colour = guide_legend(""), shape = guide_legend("")) +
  scale_colour_manual(values = my_pall[c(1,2,4)]) +
  ylab("Body weight (mg) (residuals)") + 
  xlab("Allele frequency (residuals)")
ggsave("figures/figure3.png", wt_fig_out, width = 7, height = 8, units = "in", dpi = 200)
ggsave("figures/figure3.pdf", wt_fig_out, width = 7, height = 8) 



#----

#Predictor residual plot for population growth rate outlier at edge----
gr_sites_df <- 
  all_freqs %>% 
  filter(MRK %in% gr_df$MRK) %>% 
  full_join(., trait_df, by = "POP")

gr_model_df <-   
gr_sites_df %>% 
  select(MRK, POP, M_P, mean_gr) %>%
  drop_na() %>% 
  arrange(POP) %>% 
  spread(MRK, M_P) %>%
  left_join(., Z, by = "POP") 
  
#added variable plot
select(gr_df, MRK, POS)
gr_mod <- lm(mean_gr ~  `58565` + `36490` + `158827` + `84435`, data = gr_model_df)
summary(gr_mod)

#added variable plots for growth rate
av_gr <- car::avPlots(gr_mod)

gr_av_df <- make_avdf(gr_df, av_gr, gr_model_df)

gr_av_plot <- gr_av_df %>% 
  ggplot(aes(allele_residual, trait_residual, colour = treat)) + 
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = F, colour = "grey50") +
  facet_wrap(~mrk_treat, ncol = 1) +
  guides(colour = guide_legend(""), shape = guide_legend("")) +
  scale_colour_manual(values = my_pall[c(1,2,4)]) 


edge_coeff <- filter(outlier_table, Treatment == "edge") %>% pull(Coeff.)
#summary(lm(gr_model_df$mean_gr ~ edge_resid))

edge_mrk <- 
  all_freqs %>% 
  right_join(., mutate(gr_df, MRK = as.numeric(MRK))) %>% 
  full_join(., Z, by = "POP") %>% 
  mutate(mrk_treat = paste0(treatment, " (", CHROM, " : ", POS, ")"))

select(edge_mrk, POP, MRK, M_P, treat) %>% 
  group_by(treat) %>% 
  summarise(mean(M_P))

edge_plot <- 
  edge_mrk %>% 
  ggplot(aes(treat, M_P, colour = treat, shape = treat)) +
  geom_boxplot(fill = "white", colour = "black", outlier.shape = NA) +
  geom_jitter(width = 0.1, height = 0) +
  ylab("Allele frequency") +
  xlab("Population") +
  facet_wrap(~mrk_treat, ncol = 1) +
  scale_colour_manual(values = my_pall) +
  scale_shape_manual(values=c(20, 17, 18, 15)) 
  #ylim(0,1) +
  #theme(legend.position = "n")

fig_out <- gr_av_plot + 
  ylab("Log population growth rate (residuals)")  + 
  xlab("Allele frequency (residuals)") +
  theme(legend.position = "n") + 
  edge_plot

ggsave("figures/figure4.png", fig_out, width = 8, height = 8, units = "in", dpi = 200)
ggsave("figures/figure4.pdf", fig_out, width = 8, height = 8) 


disp_edge <- trait_df %>% 
  select(POP, prop_disp) %>%
  drop_na() %>% 
  arrange(POP) %>% 
  left_join(.,edge_mrk, by = "POP")

plot(disp_edge$M_P, disp_edge$prop_disp)
plot(trait_df$mean_weight, trait_df$mean_gr)
plot(trait_df$prop_disp, trait_df$mean_gr)


##PLAYGROUND

mrk_df <- all_freqs %>% 
  filter(MRK == 14581) %>% 
  full_join(., Z, by = "POP") %>% 
  full_join(., trait_df, by = "POP") %>% 
  full_join(., weight_df, na_matches = "never")

mrk_df %>%
  select(M_P, mean_gr, treat.x) %>% 
  drop_na() %>% 
  ggplot(aes(treat.x, M_P)) +
  geom_boxplot() 

mrk_df %>%
  select(M_P, mean_gr, treat.x) %>% 
  drop_na() %>% 
  #mutate(mean_gr = sample(mean_gr)) %>%   
  ggplot(aes(M_P, mean_gr, color = treat.x)) +
  geom_point() +
  geom_smooth(aes(M_P, mean_gr), inherit.aes = F, method = "lm")  
