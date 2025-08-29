# brms analysis

# adam's analysis
library(tidyverse)
library(brms)
library(pdp)
library(ggthemes)
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
d <- read_csv('data/data_with_climate_norms.csv') |>
  dplyr::select(-vertcv, -sd.sd, -vert.sd)
glimpse(d)

# bayesian zi beta tutorial: 
# https://www.andrewheiss.com/blog/2021/11/08/beta-regression-guide/#zero-inflated-beta-regression-bayesian-style
ff <-  bf(rel_cover_exotic ~ 
            rumple +
            entropy  +
            (1|plotid) + (1|site),
          phi ~ 1 ,
            # rumple +
            # entropy  +
            # (1|plotid) + (1|site),
          zi ~ 1
            # rumple +
            # nspp_native +
            # folded_aspect +
            # entropy +
            # (1|plotid) + (1|site)
            )

brms::get_prior(ff, family = zero_inflated_beta(), data = d) -> pp

priors <- c(set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
            set_prior("normal(0, 1)", class = "b"))

if(!file.exists("data/rce_brm.rda")){
  rce_brm_zi <- brm(ff, data = d, 
                    family = "zero_inflated_beta", 
                    init = 0,
                    control = list(adapt_delta = 0.97,
                                   max_treedepth = 12),
                    iter = 2000, warmup = 1000,
                    cores = 4, seed = 1234, 
                    file = 'rec_brms',
                    prior = priors)
  save(rce_brm, file ="data/rce_brm.rda")
}else{load("data/rce_brm.rda")}
summary(rce_brm)
tidy(rce_brm) %>%
  mutate_if(is.numeric, signif, 3) %>%
  write_csv("output/rce_table.csv")
ce <- conditional_effects(rce_brm, spaghetti=T, ndraws = 200)
p <- plot(ce,plot=F, theme = theme(axis.title.y = element_blank()))
ggsave(plot = wrap_plots(p) +
  plot_annotation(title = "Relative Non-Native Abundance"), filename = "output/brmtest.png")



