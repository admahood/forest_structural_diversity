# adam's analysis
library(tidyverse)
library(ranger)
library(randomForest)
library(pdp)
library(topomicro)
library(sf)
library(terra)
dd <- read_csv("data_quality_control/forest_div_quality_control_v1.csv") %>%
  dplyr::mutate(latitude = latitude.x,
                longitude = longitude.x,
                folded_aspect = get_folded_aspect(aspect),
                cosign_aspect = cos(aspect),
                slope_aspect = slope*cosign_aspect) %>%
  dplyr::select(-ends_with(".y"), -ends_with(".x"))
plotids<- dd$plotID %>% unique()
# maybe separate slope and folded aspect
# maybe get rid of .aop
centroids <- read_csv("data/All_NEON_TOS_Plot_Centroids_V8.csv") %>%
  filter(plotID %in% plotids) %>%
  dplyr::select(latitude, longitude, plotID) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# 30 year normals of ppt, mean temp, and VPD min and max (resolution 800m) from https://prism.oregonstate.edu/normals/
bils <-list.files("data", recursive = T,
                        pattern = "annual_bil.bil$",
                        full.names = TRUE) %>%
  terra::rast() 

norms <- centroids %>%
  mutate(terra::extract(bils, vect(.), ID=F)) %>%
  rename(MAP = 3, MAT = 4, vpdmx=5, vpdmn=6) %>%
  st_set_geometry(NULL) %>%
  unique()
nrow(norms)
nrow(dd)
norms[15,]
dd[81,]

d <- dd %>%
  left_join(norms)
write_csv(d, "data/data_with_climate_norms.csv")

nrow(d)
dd

glimpse(d)
names(d)
summary(d)
predlist <- list()
predlist$all <- c("mean.max.canopy.ht.aop","rumple.aop",  "MAP", "MAT",
           "NLCD_plot", 
           "deepgap.fraction.aop", "latitude", "vpdmx", "vpdmn",
          "vertCV.aop","entropy.aop","GFP.AOP.aop","VAI.AOP.aop",
           "folded_aspect","minElev", "slope") 

predlist$only_structure <- c("mean.max.canopy.ht.aop","rumple.aop",  
                    "deepgap.fraction.aop", 
                    "vertCV.aop","entropy.aop","GFP.AOP.aop","VAI.AOP.aop") 

predlist$not_structure <- c("MAP", "MAT", "NLCD_plot", "latitude", "vpdmx", "vpdmn",
           "folded_aspect","minElev", "slope") 

resp <- c("rel_cover_exotic", "shannon_native","nspp_exotic", "nspp_native")
# higly correlated (r>80) vars: ,"q25.aop ","VCI.AOP.aop ""max.canopy.ht.aop", "cover.fraction.aop","vert.sd.aop","top.rugosity.aop",

# update methods in google doc
# add domain, site -- good random effect for LMM -  spatial regime model for lmm?
# maybe get aridity index
#"cover.fraction.aop","vert.sd.aop","top.rugosity.aop",
# super basic random forest analysis

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
    print(resp[i] |> paste(predlist[j]))
    counter <- counter + 1
  }
}
result_df %>%
write_csv( "output/rf_accuracy.csv")

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
              "NLCD_plot" = "Site",
              "minElev"="Site")

dfs<-list()
for(i in 1:length(mods)){
  dfs[[i]] <- importance(mods[[i]]) %>%
    as_tibble(rownames="variable") %>%
    mutate(response = names(mods)[i])
}


mods_all <- mods[c(1,4,7,10)]
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
  mutate(response = ifelse(response == "shannon_native_all", "Native Alpha Diversity", response),
         response = ifelse(response == "nspp_native_all", "Native Richness", response),
         response = ifelse(response == "nspp_exotic_all", "Exotic Richness", response),
         response = ifelse(response == "rel_cover_exotic_all", "Exotic Dominance", response),
         variable = variable %>% str_remove_all(".aop") %>%str_remove_all(".AOP"))


barplot <- dfs %>%
  bind_rows() %>%
  filter(response %in% c("rel_cover_exotic_all", "nspp_exotic_all", 
                         "nspp_native_all", "shannon_native_all")) %>%
  mutate(response = ifelse(response == "shannon_native_all", "Native Alpha Diversity", response),
         response = ifelse(response == "nspp_native_all", "Native Richness", response),
         response = ifelse(response == "nspp_exotic_all", "Exotic Richness", response),
         response = ifelse(response == "rel_cover_exotic_all", "Exotic Dominance", response)) %>%
  dplyr::rename(Importance = IncNodePurity) %>%
  mutate(type = lut_vars[variable],
         variable = variable %>% str_remove_all(".aop") %>%str_remove_all(".AOP") %>%
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
        legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
barplot

ggsave(plot = barplot, filename = "output/rf_fig.png", width = 7.5, height = 7, bg="white")


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
           response = ifelse(response == "nspp_exotic", "Exotic Richness", response),
           response = ifelse(response == "rel_cover_exotic", "Exotic Dominance", response))

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
  theme(axis.title.x = element_blank(),
        legend.position = c(.975, .02),
        legend.justification = c(1,0),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank());pdp_plot
  
ggsave(plot = pdp_plot,
       filename = "output/pdps_site.png", 
       bg="white", width=9.5, height=9)
 