# make a table of basic information about invasion rates
library(tidyverse)

# invasion table ===============================================================

dd <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
  dplyr::mutate(latitude = latitude.x,
                longitude = longitude.x,
                folded_aspect = topomicro::folded_aspect(aspect),
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect,
                rer = (nspp_exotic/nspp_total)*100) %>%
  dplyr::select(-ends_with(".y"), -ends_with(".x"))


glimpse(dd)
unique(dd$invaded)

invasion_by_site <- dd |>
  mutate(invaded = ifelse(invaded == 'invaded', 1, 0)) |>
  group_by(siteID, year) |>
  summarise(n = n(),
            nspp = mean(nspp_exotic),
            i = sum(invaded),
            ir = mean(rer),
            ii = mean(rel_cover_exotic)) |>
  ungroup()  |>
  mutate(i = i/n) |>
  group_by(siteID) |>
  summarise(mean_exotic_spp = mean(nspp),
            percent_invaded_plots = mean(i) * 100,
            mean_ir = mean(ir),
            mean_ii = mean(ii) * 100) |>
  ungroup()
write_csv(invasion_by_site |> mutate_if(is.numeric, round, 2),
          'output/table_X_invasion_by_site.csv')

invasion_by_site |>
  mutate(i_class = case_when(
    i == 1 ~ "100",
    i == 0 ~ "0",
    i > 0 & i < .5 ~ '0-50',
    i >= 0.5 & i < 1 ~ '50-100'
  )) |>
  pull(i_class) |>
  table()

# invasion by nlcd class =======================================================

dd <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
  dplyr::mutate(latitude = latitude.x,
                longitude = longitude.x,
                folded_aspect = topomicro::folded_aspect(aspect),
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect,
                rer = (nspp_exotic/nspp_total)*100) %>%
  dplyr::select(-ends_with(".y"), -ends_with(".x"))


glimpse(dd)
unique(dd$invaded)

invasion_by_nlcd <- dd |>
  mutate(invaded = ifelse(invaded == 'invaded', 1, 0)) |>
  group_by(NLCD_plot_des_main, year) |>
  summarise(n = n(),
            nspp = mean(nspp_exotic),
            i = sum(invaded),
            ir = mean(rer),
            ii = mean(rel_cover_exotic)) |>
  ungroup()  |>
  mutate(i = i/n) |>
  group_by(NLCD_plot_des_main) |>
  summarise(mean_exotic_spp = mean(nspp),
            percent_invaded_plots = mean(i) * 100,
            mean_ir = mean(ir),
            mean_ii = mean(ii) * 100) |>
  ungroup()

plot_counts <- dd |>
  group_by(NLCD_plot_des_main) |>
  summarise(n = n())

write_csv(invasion_by_nlcd |> mutate_if(is.numeric, round, 2) |> left_join(plot_counts), 'output/invasion_by_nlcd.csv')

# table of structural metrics ==================================================

dd |>
  dplyr::select()

invasion_by_plot <- dd |>
  mutate(invaded = ifelse(invaded == 'invaded', 1, 0)) |>
  group_by(plotID, year) |>
  summarise(n = n(),
            i = sum(invaded),
            ir = mean(rer),
            ii = mean(rel_cover_exotic)) |>
  ungroup()  |>
  mutate(i = i/n) |>
  group_by(plotID) |>
  summarise(i = mean(i),
            ir = mean(ir),
            ii = mean(ii)) |>
  ungroup()


invasion_by_plot |>
  mutate(i_class = case_when(
    i == 1 ~ "100",
    i == 0 ~ "0",
    i > 0 & i < .5 ~ '0-50',
    i >= 0.5 & i < 1 ~ '50-100'
  )) |>
  pull(i_class) |>
  table()

# structure metrics by site ====================================================
dd |>
  dplyr::select(ends_with('aop'), -sd.sd.aop) %>%
  pivot_longer(cols = names(.)) |>
  dplyr::group_by(name) |>
  summarise(mean = mean(value),
            sd = sd(value)) |>
  ungroup() |>
  mutate_if(is.numeric, signif, 2) |>
  mutate(name = str_remove_all(name, ".aop")|> str_remove_all('.AOP')) |>
  write_csv('output/structural_metric_table.csv')



