#example for Topher
library(tidyverse)

#load patch info, convert to long form, seperate prefix column into its parts
Z <- read_csv("data/Z.csv") %>%
  gather(treat, log, -prefix) %>%
  filter(log == 1) %>%
  separate(prefix, into = c("landscape", "gen", "treat_p", "id"), sep = "_") %>%
  arrange(landscape) %>%
  mutate(POP = 1:n())

#S is shuffled, C is structured
dispersal_df <- read_csv("data/F2_dispersal.csv") %>% 
  mutate(treat = ifelse(Treatment == "S","shuffled", Location)) %>% 
  mutate(landscape = as.character(Landscape)) %>% 
  select(-Landscape)

census_df <- read_csv("data/F2_census.csv") %>% 
  mutate(treat = ifelse(Treatment == "S","shuffled", Location)) %>% 
  mutate(landscape = as.character(Landscape)) %>%
  mutate(treat = as.character(treat)) %>%
  select(-Landscape)

growthrate_df <- 
  census_df %>% 
  mutate(
    growth_rate = ifelse(!is.na(Count), 
                         Count / Density, 
                         total_weight/(weight_50 / 50) / Density)
  ) %>% 
  group_by(landscape, treat) %>% 
  summarise(mean_gr = mean(growth_rate)) %>% 
  left_join(., Z, by = c("landscape" = "landscape", "treat" = "treat")) 


#to get average beetle weight
#if weight_50 is NA -> total_weight_census / count_census
#if weight_50 is NOT NA -> weight_50 / 50


weight_df <- 
  census_df %>% 
  mutate(
    mean_weight = ifelse(is.na(weight_50), 
                         total_weight / Count,
                         weight_50/50)
  ) %>% 
  group_by(treat, landscape) %>% 
  summarise(mean_weight = mean(mean_weight)) 


weight_df %>% 
  ggplot(aes(treat, mean_weight)) +
  geom_boxplot()



prop_disp_df <- 
  dispersal_df %>% 
  filter(Patch == 1, !is.na(landscape)) %>% 
  mutate(disp = (Density - Count)/Density) %>% 
  mutate(disp = ifelse (disp < 0, 0, disp),
         land_treat_loc = paste(landscape, Treatment, Location, sep = "_")) %>% 
  select(land_treat_loc, treat, disp, Density) %>% 
  spread(Density, disp) %>% 
  #mutate(prop_disp = (`40` + `10`)/2) %>% 
  mutate(prop_disp = (`40` - `10`)) %>%
  separate(land_treat_loc, c("landscape", "Treatment", "Location"), sep = "_") %>% 
  right_join(., Z, by = c("landscape" = "landscape", "treat" = "treat")) %>% 
  group_by(landscape, treat) %>% 
  summarise(prop_disp = mean(prop_disp),
            prop_40 = mean(`40`),
            prop_10 = mean(`10`)) %>% 
  drop_na() 

trait_df <- 
  full_join(growthrate_df, prop_disp_df, by = c("landscape" = "landscape", "treat" = "treat")) %>% 
  full_join(., weight_df, by = c("landscape" = "landscape", "treat" = "treat")) %>% 
  mutate(mean_weight = mean_weight*1000,
         mean_gr = log(mean_gr),
         prop_disp = qlogis(prop_disp),
         prop_40 = qlogis(prop_40),
         prop_10 = qlogis(prop_10)
         )

nrow(growthrate_df)
nrow(weight_df)
nrow(prop_disp_df)
nrow(trait_df)


my_scale <- function(x){
  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
}
final_trait_df <- select(trait_df, -prop_40, -prop_10) %>% 
  drop_na() %>% 
  ungroup() %>% 
  mutate(mean_gr_sc = my_scale(mean_gr),
    prop_disp_sc = my_scale(prop_disp),
    mean_weight_sc = my_scale(mean_weight)
  )


cov_mat <- read.table("data/baypass2/core_model_mat_omega.out")
#covariates formatted for baypass
trait_covariates <- select(final_trait_df, c(mean_gr_sc, prop_disp_sc, mean_weight_sc)) %>% 
  t()
write.table(trait_covariates, "data_out/trait_covariates.txt", col.names = FALSE, row.names = FALSE)
cov_mat_sub <- cov_mat[final_trait_df$POP, final_trait_df$POP]

write.table(cov_mat_sub, "data_out/core_model_mat_omega_traitsubset.out", col.names = FALSE, row.names = FALSE)

final_trait_df$POP
columns <- sort(c(seq(1, 2*max(Z$POP), 2)[final_trait_df$POP],
  seq(2, 2*max(Z$POP), 2)[final_trait_df$POP]
  ))

paste(columns+5, collapse = ",")


im <- scale(matrix(rnorm(1000, mean = 0, sd = 0.1), ncol = 10))
plot(cov(im), (t(im) %*% im)/nrow(im))
abline(0, 1)
