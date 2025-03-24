# brms analysis

# adam's analysis
library(tidyverse)
library(brms)
library(pdp)
make_caterpillar <- function(mod, title, file){
  tidy(mod) %>%
    mutate(term = str_remove_all(term, ".AOP") %>% str_remove_all(".aop")) %>%
    filter(is.na(group)) %>%
    filter(term != "(Intercept)", term != "year") %>%
    mutate(sig = case_when(
      conf.high * conf.low < 0 ~ "2",
      conf.high * conf.low > 0 & estimate >0 ~ "3",
      conf.high * conf.low > 0 & estimate <0 ~ "1"
    ),
    term = str_replace_all(term, "\\.", "\n"),
    term = fct_reorder(term,  estimate),
    ) %>%
    ggplot(aes(x=estimate, y=term,color = sig)) +
    geom_vline(xintercept = 0, lty=2, color = "grey") +
    geom_point(size=1) +
    geom_segment(aes(x=conf.low, xend = conf.high, yend = term), lwd=.5) +
    scale_color_manual(values = c(cols[1], 'grey30', cols[2])) +
    theme_linedraw() +
    theme(axis.title.y = element_blank(),
          legend.position = "none") +
    ggtitle(title)
  ggsave(filename = file, height = 3, width = 3.5, bg="white")
}
d <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
  dplyr::mutate(latitude = latitude.x,
                longitude = longitude.x,
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect)
glimpse(d)

if(!file.exists("data/rce_brm.rda")){
  rce_brm <- brm(rel_cover_exotic ~ 
              mean.max.canopy.ht.aop +
              rumple.aop +
              deepgap.fraction.aop +
              entropy.aop +
              VAI.AOP.aop +
              vertCV.aop +
              GFP.AOP.aop + 
              year +  (1|plotID) + (1|site), data = d, family = "zero_inflated_beta")
  save(rce_brm, file ="data/rce_brm.rda")
}else{load("data/rce_brm.rda")}
summary(rce_brm)
tidy(rce_brm) %>%
  mutate_if(is.numeric, signif, 3) %>%
  write_csv("output/rce_table.csv")
ce <- conditional_effects(rce_brm, spaghetti=T, ndraws = 200)
p <- plot(ce,plot=F, theme = theme(axis.title.y = element_blank()))
ggsave(plot = wrap_plots(p) +
  plot_annotation(title = "Exotic Dominance"), filename = "output/brmtest.png")


library(RColorBrewer);cols<-brewer.pal(n = 2, name="Set1")
if(!file.exists("data/nse_brm.rda")){
  nse_brm <- brm(nspp_exotic ~ 
                   mean.max.canopy.ht.aop +
                   rumple.aop +
                   deepgap.fraction.aop +
                   entropy.aop +
                   VAI.AOP.aop +
                   vertCV.aop +
                   GFP.AOP.aop + 
                   year  + (1|plotID) + (1|site), data = d, family = "poisson")
  save(nse_brm, file ="data/nse_brm.rda")
}else{load("data/nse_brm.rda")}


if(!file.exists("data/nsn_brm.rda")){
  nsn_brm <- brm(nspp_native ~ 
                   mean.max.canopy.ht.aop +
                   rumple.aop +
                   deepgap.fraction.aop +
                   entropy.aop +
                   VAI.AOP.aop +
                   vertCV.aop +
                   GFP.AOP.aop + 
                   year  + (1|plotID) + (1|site), data = d, family = "poisson")
  save(nsn_brm, file ="data/nsn_brm.rda")
}else{load("data/nsn_brm.rda")}

if(!file.exists("data/dvn_brm.rda")){
  dvn_brm <- brm(shannon_native ~ 
                   mean.max.canopy.ht.aop +
                   rumple.aop +
                   deepgap.fraction.aop +
                   entropy.aop +
                   VAI.AOP.aop +
                   vertCV.aop +
                   GFP.AOP.aop + 
                   year  + (1|plotID) + (1|site), data = d, family = "gamma")
  save(dvn_brm, file ="data/dvn_brm.rda")
}else{load("data/dvn_brm.rda")}
summary(dvn_brm)

make_caterpillar(nse_brm, title = "Exotic Richness", file = "output/caterpillar_nse.png")
make_caterpillar(nsn_brm, title = "Native Richness", file = "output/caterpillar_nsn.png")
make_caterpillar(rce_brm, title = "Exotic Dominance", file = "output/caterpillar_rce.png")
make_caterpillar(dvn_brm, title = "Native Diversity", file = "output/caterpillar_dvn.png")

library(grid)
library(patchwork)

ce <- conditional_effects(nse_brm, spaghetti=T, ndraws = 200)
p <- plot(ce,plot=F)
wrap_plots(p)

ce <- conditional_effects(rce_brm, spaghetti=T, ndraws = 200)
p <- plot(ce,plot=F)
wrap_plots(p)
as_tibble(ce)