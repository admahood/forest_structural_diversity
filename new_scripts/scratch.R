library(tidyverse)
d<- read_csv('data/data_with_climate_norms.csv')


ggplot(d, aes(x=rumple , y=vai, color = site, size = rer)) +
  geom_point()#+
  scale_y_log10() #+ facet_wrap(~site)
ggplot(d, aes(x=entropy, y=rer, color = site)) +
  geom_point()

d |> dplyr::select(9:23, rer, i_cat, rel_cover_exotic) %>%
  pivot_longer(cols = names(.)[1:14]) 


d |> filter(days_since_change < 10000) |> pull(plotID) |> unique() |> length()
  # ggplot(aes(x = days_since_change/365)) + 
  # geom_histogram()

d |> mutate(days_since_change = ifelse(days_since_change == 99999, NA, days_since_change)) |>
  nrow()
