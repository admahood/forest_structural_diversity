# alternate plant cover metrics extraction
library(tidyverse)
root <- paste0(getwd(), "/")

outfile_plot <- paste0(root, "data/plot_site_summaries.csv")

if(!file.exists(outfile_plot)){
  plot_level <- lapply(input_data,
                       function(x)neondiveRsity::get_diversity_info(x, scale = "plot"))%>%
    bind_rows(plot_level)
  
  site_level <- lapply(input_data,
                       function(x)neondiveRsity::get_diversity_info(x, scale = "site")) %>%
    bind_rows()
  site_plot <- bind_rows(plot_level, site_level)
  write_csv(site_plot, file = outfile_plot)
}else{site_plot <- read_csv(outfile_plot)}

sites_n_dates <- read_csv(paste0(root, "data/NEON_sites_dates_for_cover.csv"))# %>%
 # mutate_all(function(x) str_sub(x,1,4))
veg_types <- read_csv('data/field-sites.csv')
plot_data_table <- read_csv("new_scripts/output/plot_data_table.csv")
lidar_data <- read_csv("new_scripts/output/lidar_structural_metrics.csv")

# comparison of neondiveRsity vs existing script ===============================

previous_output <- read_csv("new_scripts/output/structural_metrics_by_plot.csv")

bona_plots <- previous_output %>% filter(siteID == "ABBY") %>% pull(plotID)
previous_output %>% 
  filter(siteID == "ABBY") %>%
  dplyr::select(species_richness, plotID, year, exotic_cover) %>%
  arrange(plotID)
site_plot %>% 
  filter(site == "ABBY" & year %in% c(2017,2019) & plotID %in% bona_plots) %>%
  dplyr::select(nspp_total, plotID, year, cover_exotic)
