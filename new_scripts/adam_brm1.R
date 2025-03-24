# adam's brms analysis
library(tidyverse)
library(topomicro)
library(sf)
library(terra)
library(brms)
library(performance)

# nlcd cover class
# hurdle model? investigate zeros
# mean aridity index?


d <- read_csv("data/data_with_climate_norms.csv") |>
  dplyr::rename(mean.max.canopy.ht = mean.max.canopy.ht.aop,
                max.canopy.ht = max.canopy.ht.aop,
                rumple = rumple.aop,
                deepgap.fraction = deepgap.fraction.aop,
                cover.fraction = cover.fraction.aop,
                top.rugosity = top.rugosity.aop,
                vert.sd = vert.sd.aop,
                vertCV = vertCV.aop,
                sd.sd = sd.sd.aop,
                entropy = entropy.aop,
                GFP= GFP.AOP.aop,
                VAI = VAI.AOP.aop,
                VCI = VCI.AOP.aop,
                q25 = q25.aop,
                q50 = q50.aop)


ggplot(d, aes(y=NLCD_plot_des_main, x=nspp_notexotic)) +
  geom_boxplot()
ggplot(d, aes(y=NLCD_plot_des_main, x=nspp_exotic)) +
  geom_boxplot()
# doing a princomp
pca <- d |>
  dplyr::select(VAI, VCI, vert.sd, q25, q50, GFP, entropy, rumple, vertCV, 
                top.rugosity, deepgap.fraction, mean.max.canopy.ht, max.canopy.ht,
                cover.fraction) |>
  prcomp(scale=TRUE)

pca$rotation |>
  as_tibble(rownames = "var") |>
  ggplot(aes(x=PC1, y=PC2)) +
  geom_text_repel(aes(label = var))



# rn <- brm(nspp_notexotic ~ entropy + VAI + rumple + mean.max.canopy.ht + GFP + MAP + MAT + vpdmx + (1|site/year/plotID), 
#           family = 'poisson',
#           data = d)
# summary(rn)
# 
# 
# sn <- brm(shannon_notexotic ~ s(entropy) + s(VAI) + s(rumple) + 
#             mean.max.canopy.ht + MAT + (1|site/year/plotID), 
#           family = 'lognormal',cores = 4,
#           data = d)
# summary(sn)
# marginal_effects(sn)
# performance::check_model(sn)
library(splines)
ir <- glmmTMB::glmmTMB(rer ~ #I((rer + 0.001)/100) ~ 
                         ns(rumple, 2) + VAI + vpdmx + year + (1|site/plotID), 
                 family = tweedie(),
                 data = d |> mutate_if(is.numeric,scale))
summary(ir)
performance::r2(ir)
performance::check_model(ir)
diagnose(ir)
ggeffects::ggpredict(ir, terms = c('entropy')) |> plot()
ggeffects::ggpredict(ir, terms = c('VAI')) |> plot()
ggeffects::ggpredict(ir, terms = c('mean.max.canopy.ht')) |> plot()


ii <- glmmTMB::glmmTMB(rel_cover_exotic ~ entropy + VAI + rumple + vpdmx + (1|site/plotID), 
                       family = tweedie(),
                       data = d|> mutate_if(is.numeric,scale))
summary(ii)
performance::r2(ii)
performance::check_model(ii)
ggeffects::ggpredict(ii, terms = c('entropy')) |> plot()
ggeffects::ggpredict(ii, terms = c('VAI')) |> plot()
ggeffects::ggpredict(ii, terms = c('mean.max.canopy.ht')) |> plot()

# ir <- brm(rer ~ entropy + MAT + (1|site/plotID/year), 
#           family = 'gaussian', cores = 4,
#           data = d)
# save(ir, file = 'data/ir.rda')
# summary(ir)
# conditional_effects(ir)
performance::check_model(ir)

ii <- brm(rel_cover_exotic ~ entropy + VAI + rumple + mean.max.canopy.ht + GFP + MAP + MAT + vpdmx + (1|site/plotID), 
          family = 'gaussian',
          data = d)
summary(ii)
