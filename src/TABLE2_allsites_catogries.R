library(tidyverse)
library(magrittr)

categories <- c(
  "synonymous",
  "missense",
  "stop gained",
  "stop lost",
  "stop retained",
  "start lost",
  "splice region&synonymous",
  "start lost&splice region",
  "3 prime UTR",
  "5 prime UTR",
  "intron",
  "splice acceptor",
  "splice donor",
  "splice region",
  "upstream",
  "downstream",
  "intergenic"
  )

impact_df <- full_join(
  read_csv("data_out/outlier_categories.csv"),
  read_csv("data_out/conserved_categories.csv"), 
  by = "category") %>% 
  full_join(., 
          read_csv("data_out/highimpact_categories.csv"),
          by = "category") %>% 
  set_colnames(c("category", "outlier percent", "conserved percent" , "high impact percent")) %>% 
  replace_na(list("outlier percent" = 0, "conserved percent" = 0, "high impact percent" = 0)) %>% 
  mutate(category = factor(category, levels = categories)) %>% 
  arrange(factor(category, levels = categories))


print(impact_df)

View(impact_df)
write_csv(impact_df, "data_out/allsites_categories.csv")

