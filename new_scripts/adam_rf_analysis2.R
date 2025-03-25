# trying to do a good random forest analysis
library(tidyverse)
library(ranger)
# library(randomForest)
library(tidymodels)
library(pdp)
library(ggpubr)
library(ggrepel)
library(topomicro)
library(sf)
library(terra)
library(purrr)

lut_vars <- c("mean.max.canopy.ht"="Structure",
              "rumple"="Structure",
              "map"="Climate", "mat"="Climate",
              "deepgap.fraction"="Structure", 
              "latitude"="Site", "site"="Site", 
              "vpdmx"="Climate", "vpdmn"="Climate",
              "vertcv"="Structure",
              "entropy"="Structure",
              "gfp"="Structure",
              "vai"="Structure",
              "folded_aspect"="Site",
              "slope" = "Site",
              "nlcd_plot_des_main" = "Site",
              "minelev"="Site")

if(!file.exists('data/data_with_climate_norms.csv')){
  dd <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
    dplyr::mutate(latitude = latitude.x,
                  longitude = longitude.x,
                  folded_aspect = topomicro::folded_aspect(aspect),
                  cosign_aspect = cos(aspect),
                  slope_aspect = slope*cosign_aspect,
                  rer = (nspp_exotic/nspp_total)*100) %>%
    dplyr::select(-ends_with(".y"), -ends_with(".x"))
  
  
  plotids<- dd$plotID %>% unique()

  centroids <- read_csv("data/All_NEON_TOS_Plot_Centroids_V8.csv") %>%
    filter(plotID %in% plotids) %>%
    dplyr::select(latitude, longitude, plotID) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  bils <-list.files("data/prism", 
                    pattern = "annual_bil.bil$",
                    recursive = T,
                    full.names = TRUE) %>%
    terra::rast() 
  
  norms <- centroids %>%
    mutate(terra::extract(bils, vect(.), ID=F)) %>%
    rename(MAP = 3, MAT = 4, vpdmx=5, vpdmn=6) %>%
    st_set_geometry(NULL) %>%
    unique()
  
  d <- dd %>%
    left_join(norms) |>
    mutate(i_cat = ifelse(invaded =="invaded", 1, 0) |> as.factor()) |>
    rename_all(function(x) str_to_lower(x) |> str_remove_all(".aop"))
  write_csv(d, "data/data_with_climate_norms.csv")
}else{d<- read_csv('data/data_with_climate_norms.csv') |>
  mutate(i_cat = as.factor(i_cat))}

# listing out the predictors
predlist <- list()
predlist$all <- c("mean.max.canopy.ht","rumple",  "MAP", "MAT",
                  "nlcd_plot_des_main", "deepgap.fraction", "latitude", "vpdmx",
                  "vpdmn", "vertCV","entropy","GFP","VAI",
                  "folded_aspect","minElev", "slope") |> str_to_lower() 

predlist$only_structure <- c("mean.max.canopy.ht","rumple", "deepgap.fraction", 
                             "vertCV","entropy","GFP","VAI") 

predlist$not_structure <- c("MAP", "MAT", "NLCD_plot", "latitude", "vpdmx", "vpdmn",
                            "folded_aspect","minElev", "slope") 

resp <- c("rel_cover_exotic", "shannon_native","rer", "nspp_native", 'i_cat')
# higly correlated (r>80) vars: ,"q25.aop ","VCI.AOP.aop ""max.canopy.ht.aop", "cover.fraction.aop","vert.sd.aop","top.rugosity.aop",

# update methods in google doc
# add domain, site -- good random effect for LMM -  spatial regime model for lmm?
# maybe get aridity index
#"cover.fraction.aop","vert.sd.aop","top.rugosity.aop",
# super basic random forest analysis
# big random forest function ===================================================
drf <- d |>
  dplyr::filter(!nlcd_plot_des_main %in% c("Pature/hay"))
preds <- predlist$all

mod_form <- as.formula(paste("i_cat ~", paste(preds, collapse = " + ")))
mod_form_ir <- as.formula(paste("rer ~", paste(preds, collapse = " + ")))
mod_form_ii <- as.formula(paste("rel_cover_exotic ~", paste(preds, collapse = " + ")))
mod_form_sn <- as.formula(paste("shannon_notexotic ~", paste(preds, collapse = " + ")))
mod_form_nr <- as.formula(paste("nspp_notexotic ~", paste(preds, collapse = " + ")))

rs_obj <- group_vfold_cv(drf, group = plotid,v=10 )

holdout_results_pi <- function(splits, mod_form){
  rf_fit <-
    rand_forest(trees = 1500) %>%
    set_engine("ranger") %>%
    set_mode("classification") %>%
    fit(formula = mod_form, data = analysis(splits))
  
  holdout <- assessment(splits)
  res <- broom::augment(rf_fit, new_data = holdout)
  acc <- yardstick::accuracy(res, truth = as.character(mod_form[2]), .pred_class)
  acc
}

holdout_results_r <- function(splits, mod_form){
  rf_fit <-
    rand_forest(trees = 1500) %>%
    set_engine("ranger") %>%
    set_mode("regression") %>%
    fit(formula = mod_form, data = analysis(splits))
  
  holdout <- assessment(splits)
  res <- broom::augment(rf_fit, new_data = holdout)
  rsq <- yardstick::rsq(res, truth = as.character(mod_form)[2], .pred)
  
  
  # Return the assessment data set with the additional columns
  rsq
}

fit_mods_c <- function(splits, mod_form){
  rf_mod <- 
    rand_forest(trees = 1500) %>% 
    set_engine("ranger", importance = 'permutation') %>% 
    set_mode("classification")
  rf_fit <- 
    rf_mod %>% 
    fit(mod_form, data = analysis(splits))
}

fit_mods_r <- function(splits, mod_form){
  rf_mod <- 
    rand_forest(trees = 1500) %>% 
    set_engine("ranger", importance = 'permutation') %>% 
    set_mode("regression")
  rf_fit <- 
    rf_mod %>% 
    fit(mod_form, data = analysis(splits))
}

get_imps <- function(fit) {
  as_tibble(fit$fit$variable.importance, rownames = 'var')
}

rs_obj$fits_pi <- map(rs_obj$splits, fit_mods_c, mod_form)
rs_obj$fits_ii <- map(rs_obj$splits, fit_mods_r, mod_form_ii)
rs_obj$fits_ir <- map(rs_obj$splits, fit_mods_r, mod_form_ir)
rs_obj$fits_sn <- map(rs_obj$splits, fit_mods_r, mod_form_sn)
rs_obj$fits_nr <- map(rs_obj$splits, fit_mods_r, mod_form_nr)

rs_obj$imps_pi <- map(rs_obj$fits_pi, get_imps)
rs_obj$imps_ii <- map(rs_obj$fits_ii, get_imps)
rs_obj$imps_ir <- map(rs_obj$fits_ir, get_imps)
rs_obj$imps_sn <- map(rs_obj$fits_sn, get_imps)
rs_obj$imps_nr <- map(rs_obj$fits_nr, get_imps)

rs_obj$results_pi <- map(rs_obj$splits, holdout_results_pi, mod_form)
rs_obj$results_ii <- map(rs_obj$splits, holdout_results_r, mod_form_ii)
rs_obj$results_ir <- map(rs_obj$splits, holdout_results_r, mod_form_ir)
rs_obj$results_sn <- map(rs_obj$splits, holdout_results_r, mod_form_sn)
rs_obj$results_nr <- map(rs_obj$splits, holdout_results_r, mod_form_nr)

# rs_obj$accuracy_pi <- map_dbl(rs_obj$results_pi, function(x) mean(x$correct))
# summary(rs_obj$accuracy_pi)

# save(rs_obj, file = "rs_obj.rda")

p_r2 <- bind_rows(
  rs_obj$results_ii |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "c. ii", metric = 'R2'),
  rs_obj$results_ir |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "b. ir", metric = 'R2'),
  rs_obj$results_nr |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "d. nr", metric = 'R2'),
  rs_obj$results_sn |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "e. sn", metric = 'R2'),
  rs_obj$results_pi |> bind_rows() |> dplyr::select(est = .estimate) |> mutate(respo = "a. pi", metric = 'Acc.')) |>
  ggplot(aes(x=respo, y=est, fill = metric)) +
  geom_boxplot() +
  ggtitle("F. Accuracy") +
  theme_bw() + 
  ylab("Estimate") +
  theme(axis.title.x = element_blank(),
        legend.background = element_rect(fill=NA),
        legend.position = c(0,0),
        legend.justification = c(0,0))
# ggsave('output/r2_boxplot.png', width = 4, height = 4, bg = 'white')

df_imp <- rs_obj$imps_pi |>  bind_rows() |> mutate(respo = 'A. P(Invasion)') |>
  bind_rows(rs_obj$imps_ii |> bind_rows() |> mutate(respo = 'C. Invasion Impact')) |>
  bind_rows(rs_obj$imps_ir |> bind_rows() |> mutate(respo = 'B. Invasion Rate')) |>
  bind_rows(rs_obj$imps_sn |> bind_rows() |> mutate(respo = 'E. Native Diversity')) |>
  bind_rows(rs_obj$imps_nr |> bind_rows() |> mutate(respo = 'D. Native Richness')) |>
  group_by(respo, var) |>
  summarise(mean = mean(value),
            sd = sd(value)) |>
  ungroup() |>
  mutate(cat = lut_vars[var],
         var = ifelse(var == "nlcd_plot_des_main", 'NLCD', var),
         var = ifelse(var == "folded_aspect", 'aspect', var),
         var = ifelse(var == "mean.max.canopy.ht", 'max_ht', var),
         var = ifelse(var == "deepgap.fraction", 'dgf', var),
         var = ifelse(nchar(var) == 3, str_to_upper(var), var),
         var = fct_reorder2(var, mean, cat))

p <- lapply(unique(df_imp$respo),function(x){
  pp <- df_imp |>
    filter(respo == x) |>
    ggplot(aes(y=var)) +
    geom_bar(aes(x=mean, fill = cat), stat='identity', color = 'black') +
    geom_segment(aes(x=mean-sd, xend = mean+sd), linewidth = 1) +
    ggtitle(x) +
    scale_fill_brewer(palette = "Dark2") +
    # facet_wrap(~cat, scales = 'free_y', ncol = 1) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.justification = c(1,0))
  if(x == "C. Invasion Impact") pp <- pp + theme(legend.position = c(1,0),
                                              legend.title = element_blank(),
                                              legend.background = element_rect(fill = NA),
                                              legend.justification = c(1,0))
  if(x %in% c("D. Native Richness", "A. P(Invasion)")) return(pp) else pp <- pp + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
)
ggarrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p_r2, ncol = 3, nrow = 2) |>
  annotate_figure(bottom = "Importance")
ggsave("output/new_imp_fig.png", width =8, height = 7, bg = 'white')
# pdps - doing a loop

pdplist <- list()
cc <- 1
for(i in 1:nrow(rs_obj)){
  rf_fit <- rs_obj$fits[[i]]
  splits <- rs_obj$splits[[i]]
for(k in 1:length(preds)){
  if(preds[k] != "NLCD_plot_des_main"){
    pdplist[[cc]] <- 
      pdp::partial(rf_fit,
                   grid.resolution = 30,
                   train = analysis(splits),
                   pred.var = preds[k]) |>
      mutate(rep = as.character(i),
             name = preds[k]) |>
      dplyr::rename(value = 1)
  cc <- cc +1
  }else{print("skip")}
  print(preds[k])
}}

lapply(pdplist, function(x)dplyr::mutate(x, value = as.numeric(value))) |>
  bind_rows() |>
  ggplot(aes(x=value, y=yhat, group = rep)) +
  geom_line() +
  facet_wrap(~name, scales = 'free')

### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
do_rf <- function(drf){
  cell_split <- group_initial_split(drf, group = siteID)
  cell_train <- training(cell_split)
  cell_test  <- testing(cell_split)
  rf_mod <- 
    rand_forest(trees = 5000) %>% 
    set_engine("ranger") %>% 
    set_mode("regression")
  rf_fit <- 
    rf_mod %>% 
    fit(LDMC ~ ., data = cell_train)
  rf_training_pred <- 
    predict(rf_fit, cell_train) %>% 
    # Add the true outcome data back in
    bind_cols(cell_train %>% 
                dplyr::select(LDMC))
  train_rsq <- rf_training_pred %>%                # training set predictions
    yardstick::rsq(truth = LDMC, .pred) |>
    dplyr::rename(train_rsq = 3)
  
  rf_testing_pred <- 
    predict(rf_fit, cell_test) %>% 
    bind_cols(cell_test %>% dplyr::select(LDMC))
  test_rsq <- rf_testing_pred %>%                   # test set predictions
    yardstick::rsq(truth = LDMC, .pred) |>
    dplyr::rename(test_rsq = 3) |>
    left_join(train_rsq)
  return(test_rsq)
}


mods <- list()
counter <- 1
result_df <- data.frame(resp=NA, r2=NA, mse=NA, predlist=NA)

for(i in 1:length(resp)) {
  for(j in 1:length(predlist)){
    mods[[counter]]  <- formula(paste(resp[i], " ~ ",
                                      paste(predlist[[j]], 
                                            collapse = " + "))) %>%
      randomForest(data =d, ntree = 1000)
    result_df[counter,1] <- resp[i]
    result_df[counter,2] <- mods[[counter]]$rsq %>% median() %>% signif(3)
    result_df[counter,3] <- mods[[counter]]$mse %>% median() %>% signif(3)
    result_df[counter,4] <- names(predlist)[j]
    print(resp[i])
    counter <- counter + 1
  }
}
result_df %>%
write_csv("output/rf_accuracy.csv")

names(mods) <- result_df %>%
  mutate(name = str_c(resp, "_", predlist)) %>%
  pull(name)
# basically, elevation and lat are really good predictors of general diversity 
# (see plots for turnover and nestedness)
# lapply(mods, plot)


# multipanel importance plot ====================
lut_vars <- c("mean.max.canopy.ht.aop"="Structure",
              "rumple.aop"="Structure",
              "MAP"="Climate", "MAT"="Climate",
              "deepgap.fraction.aop"="Structure", 
              "latitude"="Site", "site"="Site", 
              "vpdmx"="Climate", "vpdmn"="Climate",
              "vertCV.aop"="Structure",
              "entropy.aop"="Structure",
              "GFP.AOP.aop"="Structure",
              "VAI.AOP.aop"="Structure",
              "folded_aspect"="Site",
              "slope" = "Site",
              "NLCD_plot_des_main" = "Site",
              "minElev"="Site")

dfs<-list()
for(i in 1:length(mods)){
  dfs[[i]] <- importance(mods[[i]]) %>%
    as_tibble(rownames="variable") %>%
    mutate(response = names(mods)[i])
}


# mods_all <- mods[c(1,4,7,10)]
mods_all <- mods
impvars_df<- list()
for(i in 1:length(names(mods_all))){
  impvars_df[[i]] <- importance(mods_all[[i]]) %>%
    as_tibble(rownames="variable") %>%
    arrange(desc(IncNodePurity)) %>%
    mutate(rank = rank(IncNodePurity),
           response = names(mods_all)[i])
}
impvars<- bind_rows(impvars_df) %>%
  filter(variable != "NLCD_plot") %>%
  pull(variable) %>%unique()

avg_rank <-
  impvars_df %>%
  bind_rows() %>%
  group_by(variable) %>%
  summarise(mean_rank_imp = mean(rank)) %>%
  ungroup() %>%
  mutate(type = lut_vars[variable]) %>%
  arrange(type, desc(mean_rank_imp)) %>%
  mutate(rank = letters[1:16],
         variable = str_remove_all(variable, ".aop"),
         variable = str_remove_all(variable, ".AOP"))

topranks <-
  impvars_df %>%
  bind_rows() %>%
  mutate(type = lut_vars[variable],
         toprank = ifelse(rank > 11, "top5", "bottom11")) %>%
  dplyr::select(toprank, variable, response, Importance = IncNodePurity) %>%
  filter(toprank == "top5")%>%
  mutate(response = ifelse(response == "shannon_native_all", "E. Native Alpha Diversity", response),
         response = ifelse(response == "nspp_native_all", "D. Native Richness", response),
         response = ifelse(response == "i_cat_all", "A. P(Invasion)", response),
         response = ifelse(response == "rer_all", "B. Invasion Rate", response),
         response = ifelse(response == "rel_cover_exotic_all", "C. Invasion Impact", response),
         variable = variable %>% str_remove_all(".aop") %>% 
           str_remove_all(".AOP") |> str_remove_all("_plot_des_main"))


barplot <- dfs %>%
  bind_rows() %>%
  filter(response %in% c("rel_cover_exotic_all", "rer_all", 'i_cat_all',
                         "nspp_native_all", "shannon_native_all")) |>
  mutate(response = ifelse(response == "shannon_native_all", "E. Native Alpha Diversity", response),
         response = ifelse(response == "nspp_native_all", "D. Native Richness", response),
         response = ifelse(response == "i_cat_all", "A. P(Invasion)", response),
         response = ifelse(response == "rer_all", "B. Invasion Rate", response),
         response = ifelse(response == "rel_cover_exotic_all", "C. Invasion Impact", response)) %>%
  dplyr::rename(Importance = IncNodePurity) %>%
  mutate(type = lut_vars[variable],
         variable = variable %>% str_remove_all(".aop") %>%str_remove_all(".AOP") %>%
           str_remove_all("_plot_des_main") |>
           fct_reorder2(Importance ,type),
         type = factor(type)) %>%
  ggplot() +
  geom_bar(stat = "Identity",aes(x = Importance, y=variable, fill = type)) +
  geom_bar(data = topranks, stat = "Identity", color ="black",
           fill="transparent",aes(x = Importance, y=variable)) +
  facet_wrap(~response, scales = "free_x", nrow =2) +
  scale_fill_brewer(palette = "Dark2", name = "Type") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        legend.background = element_rect(fill=NA),
        legend.position = c(.95,0.1),
        legend.justification = c(1,0),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
barplot

ggsave(plot = barplot, filename = "output/rf_fig.png", width = 8.5, height = 7, bg="white")


mods_all <- mods[c(1,4,7,10)]
if(!file.exists("data/pdp_df.csv")){
  df_partial <- list()
  for(i in 1:length(resp)){
    pdf<-list()
    for(j in unique(impvars)) {
      pdf[[j]] <- partial(mods_all[[i]], j) %>%
      mutate(variable = j) %>%
      dplyr::rename(value = j)
    }
  
  
    df_partial[[i]]<- bind_rows(pdf) %>%
      mutate(response = resp[i])
  }
  
  
  
  pdp_df <- bind_rows(df_partial) %>%
    group_by(response) %>%
    mutate(yhat = scale(yhat) %>% as.numeric()) %>%
    ungroup() %>%
    mutate(Type = lut_vars[variable],
           variable = variable %>% str_remove_all(".aop") %>% str_remove_all(".AOP")) %>%
    mutate(response = ifelse(response == "shannon_native", "Native Alpha Diversity", response),
           response = ifelse(response == "nspp_native", "Native Richness", response),
           response = ifelse(response == "rer", "Exotic Relative Richness", response),
           response = ifelse(response == "rel_cover_exotic", "Exotic Relative Cover", response))

write_csv(pdp_df, "data/pdp_df.csv")
}else(pdp_df <- read_csv("data/pdp_df.csv"))

pdp_plot <- pdp_df %>%
  left_join(avg_rank) %>%
  mutate(variable = paste0(rank, ". ", variable)) %>%
  ggplot(aes(x=value, y=yhat, color = response)) +
  geom_line(key_glyph = "timeseries", lwd=1) +
  facet_wrap(~variable, scales = "free", nrow = 4) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  ylab("Standardized Values\n                         Structure                                                       Topography                            Climate") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = c(.98, .02),
        legend.justification = c(1,0),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank());pdp_plot
  
ggsave(plot = pdp_plot,
       filename = "output/pdps_site.png", 
       bg="white", width=9.5, height=9)
 