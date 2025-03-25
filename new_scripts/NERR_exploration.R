library(tidyverse)

dd <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
  dplyr::mutate(latitude = latitude.x,
                longitude = longitude.x,
                folded_aspect = topomicro::folded_aspect(aspect),
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect,
                rer = (nspp_exotic/nspp_total)*100) %>%
  dplyr::select(-ends_with(".y"), -ends_with(".x"))

# do it by site too

dd |>
  mutate(ent_class = case_when(entropy.aop > .6 ~ "High Entropy",
                               entropy.aop <= 0.6 ~ 'Low Entropy')) |> #lm(nspp_exotic ~ nspp_notexotic*ent_class, data = _) |> summary()
  ggplot(aes(y=nspp_exotic, x = nspp_notexotic)) +
  geom_point() +
  facet_wrap(~ent_class) +
  geom_smooth(method = 'lm')

dd |>
  mutate(ent_class = case_when(rumple.aop > 2.4 ~ "1High Rumple",
                               rumple.aop <= 2.4 ~ '0Low Rumple')) |> #lm(nspp_exotic ~ nspp_notexotic*ent_class, data = _) |> summary()
  ggplot(aes(y=nspp_exotic, x = nspp_notexotic)) +
  geom_hex() +
  facet_wrap(~ent_class) +
  geom_smooth(method = 'lm')

dd |>
  mutate(ent_class = case_when(VAI.AOP.aop > 5 ~ "1High VAI",
                               VAI.AOP.aop <= 5 ~ '0Low VAI')) |>#lm(nspp_exotic ~ nspp_notexotic*ent_class, data = _) |> summary()
  ggplot(aes(y=nspp_exotic, x = nspp_notexotic)) +
  geom_point() +
  facet_wrap(~ent_class) +
  geom_smooth(method = 'lm')
