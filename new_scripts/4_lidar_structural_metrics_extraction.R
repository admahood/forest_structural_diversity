
set.seed(123)
library(corrplot)
library(ggplot2)
library(RColorBrewer)

library(GGally)
library(caTools)
#library(randomForest)
library(caret)
library(varSelRF)
library(leaps)
library(MASS)
library(Metrics)
library(mFD)
library(lidR)
library(gstat)
library(terra)
library(raster)


#wd = "E:/NEON_Data/Structural_diversity/LiDAR_Data/plot_level_LiDAR/"

# data directory
wd = "E:/NEON_Data/Structural_diversity/LiDAR_Data/plot_level_LiDAR/"

# setting working environment
setwd(wd)


# create empty data frame
structural_diversity =  data.frame(matrix(ncol = 20, nrow = 0))




# read all file names with .las/.laz extension

file_list = dir(wd, pattern = ".laz", full.names = FALSE, ignore.case = TRUE)
names_to_cloud = substr(file_list, 1, nchar(file_list)-4) # get the file name without extension

# loop thorugh all plot scale .laz files

for (i in 1:length(names_to_cloud)) #length(names_to_cloud
  
{
  print(i)
  

  split_name = unlist(strsplit(names_to_cloud[i], "[_]"))
  
  options(digits = 11)
  #get plot center coordinates from the file name
  x = as.numeric(split_name[[6]][1])
  y = as.numeric(split_name[[7]][1])
  
  # get site name and the year and month of data collection from the file name
  site = split_name[4]
  year_mo = split_name[5]
  
  #why dropping z below zero before correcting for elevation. 
  #Probably its ok but if a site is near sea level would remove real points. However, in this study we removed all points below zero assuming no location below sea level 
  #associated in this study.
  
  point_cloud = lidR::readLAS(file.path(paste0(file_list[i])), filter = "-drop_z_below 0 ")
  

  plot_size = 20
  plot_rad = plot_size/2
  
  try (
  data.20m <- 
    clip_rectangle(point_cloud, 
                   xleft = (x - plot_rad), ybottom = (y - plot_rad),
                   xright = (x + plot_rad), ytop = (y + plot_rad))
  
  , silent = TRUE)

  #when extract lidar data at 20m x 20 x plot level, in some locations there are no points available.
  #I am not sure why. May need to investigate those plots separately. However, in this loop I skipped those plots from further processing.
  #Are the plot centroids built into the file name? Is there a chance that the centroid is off?
  
  len_rec = data.20m@header@PHB$`Number of point records`
  print(len_rec)
  
  # check point records to see if there is enough number of lidar points within the plot area for metric calculations
  
  if (len_rec<= 1 || is.nan(len_rec)== TRUE)
  {
    next
  }
  
  # In some cases, there are no ground points available to generate the DTM. So, have to skip those plot as well from further processing.
  
  #drop out finer outliers with an IVF filter
  #ivf = isolated voxels filter, which finds points with only a few other points
  #in their surrounding x * x * x = x^3 voxel neighborhood
  #res = voxel resolution, n = maximal number of other points in the 27 voxels
  #similar to lasnoise from LAStools https://rapidlasso.com/lastools/lasnoise/
  data.20m <- classify_noise(data.20m, ivf(res=3,n=0))
  #remove points with a classification (18) of noise 
  data.20m <- filter_poi(data.20m, Classification != LASNOISE)
  
  try(
  dtm <- lidR::grid_terrain(las = data.20m,
                            res = 1, # was 0.5 and 0.25 before
                            algorithm = tin())
  , silent = TRUE)
  
  
  # dtm <- grid_terrain(data.20m, 1, kriging(k = 10L)). Use above algorithm to generate dtm.
 
  data.20m <- normalize_height(data.20m, dtm)
  
  # plot(data.20m)

  ## chm data includes all heights from 0 to max chm because we are interested in all biomes including shrubs and grass
  chm <- grid_canopy(data.20m, res = 1, dsmtin())  

  #Liz turned chm smoothing off because it was messing up metric functions below (VCI, Rumple, Entropy, etc.)
  #pixels_smooth_2 <- 9 #round(((1/chm_res)-1)/2)*1 + 1 # used smoothed parameter as 1. Instead of this I used 9 as the window size
  
  #weights <- matrix(1, nrow = pixels_smooth_2, ncol = pixels_smooth_2)
  #this smoothing method seems to be taking off more than outliers > 80 m
  #smoothing method only impacts the chm raster and not the point cloud, which is used in many metrics
  #chm <- terra::focal(x = chm, weights, fun = mean,na.rm=TRUE ) # removed all values of "NA" in the smoothed products
  
  # plot(chm)
  ###########################################################################
  ########## All below indices are from https://www.neonscience.org/resources/learning-hub/tutorials/structural-diversity-discrete-return
  ###########################################################################
  
  #MEAN OUTER CANOPY HEIGHT (MOCH)
  #calculate MOCH, the mean CHM height value
  mean.max.canopy.ht <- mean(chm@data@values, na.rm = TRUE) 
  
  #MAX CANOPY HEIGHT
  #calculate HMAX, the maximum CHM height value
  max.canopy.ht <- max(chm@data@values, na.rm=TRUE) 
  
  #use Forest Service definition of trees and remove plots with vegetation < 1.4 m tall
  if(max.canopy.ht < 1.4){
    next
  }
  #RUMPLE
  #calculate rumple, a ratio of outer canopy surface area to 
  #ground surface area (1600 m^2)
  #rumple is broken with smoothed raster but works with unsmoothed
  rumple <- rumple_index(chm, y = NULL, z = NULL) 
  
  #TOP RUGOSITY
  #top rugosity, the standard deviation of pixel values in chm and 
  #is a measure of outer canopy roughness
  top.rugosity <- sd(chm@data@values, na.rm = TRUE) 
  
  #DEEP GAPS & DEEP GAP FRACTION
  #number of cells in raster (also area in m2)
  cells <- length(chm@data@values) 
  chm.0 <- chm
  chm.0[is.na(chm.0)] <- 0 #replace NAs with zeros in CHM
  #create variable for the number of deep gaps, 1 m^2 canopy gaps
  zeros <- which(chm.0@data@values == 0) 
  deepgaps <- length(zeros) #number of deep gaps
  #deep gap fraction, the number of deep gaps in the chm relative 
  #to total number of chm pixels
  deepgap.fraction <- deepgaps/cells 
  
  #COVER FRACTION
  #cover fraction, the inverse of deep gap fraction
  cover.fraction <- 1 - deepgap.fraction 
  
  #remove ground points so that only include vegetation in structural diversity metrics that use point cloud directly 
  #chm has zeros to maintain area
  data.20m <- filter_poi(data.20m, Classification != 2L)
  
  #HEIGHT SD
  #height SD, the standard deviation of height values for all points
  #in the plot point cloud
  vert.sd <- cloud_metrics(data.20m, sd(Z, na.rm = TRUE)) 
  meanH <- cloud_metrics(data.20m, mean(Z, na.rm = TRUE)) 
  vertCV <- vert.sd / meanH
  
  #quantiles 
  Q25 <- cloud_metrics(data.20m, quantile(Z, c(0.25), na.rm = TRUE))
  Q50 <- cloud_metrics(data.20m, quantile(Z, c(0.5), na.rm = TRUE))
  
  #SD of VERTICAL SD of HEIGHT
  #rasterize plot point cloud and calculate the standard deviation 
  #of height values at a resolution of 1 m^2
  sd.1m2 <- grid_metrics(data.20m, sd(Z), 1)
  #standard deviation of the calculated standard deviations 
  #from previous line 
  #This is a measure of internal and external canopy complexity
  sd.sd <- sd(sd.1m2[,3], na.rm = TRUE) 
  
  
  #some of the next few functions won't handle NAs, so we need 
  #to filter these out of a vector of Z points
  Zs <- data.20m@data$Z
  Zs <- Zs[!Zs < 0] #must remove negative numbers to make VCI/Entropy to work
  Zs <- Zs[!is.na(Zs)]
  
  #ENTROPY 
  #entropy, quantifies diversity & evenness of point cloud heights 
  #by = 1 partitions point cloud in 1 m tall horizontal slices 
  #ranges from 0-1, with 1 being more evenly distributed points 
  #across the 1 m tall slices 
  entro <- entropy(Zs, by = 1, zmax = 80) 
  
  #GAP FRACTION PROFILE 
  #gap fraction profile, assesses the distribution of gaps in the 
  #canopy volume 
  #dz = 1 partitions point cloud in 1 m horizontal slices 
  #z0 is set to a reasonable height based on the age and height of 
  #the study sites 
  gap_frac <- gap_fraction_profile(Zs, dz = 1, z0=1.4) 
  #defines gap fraction profile as the average gap fraction in each 
  #1 m horizontal slice assessed in the previous line
  GFP.AOP <- mean(gap_frac$gf) 
  
  #VAI
  #leaf area density, assesses leaf area in the canopy volume 
  #k = 0.5 is a standard extinction coefficient for foliage 
  #dz = 1 partitions point cloud in 1 m horizontal slices 
  #z0 is set to the same height as gap fraction profile above
  LADen<-LAD(Zs, dz = 1, k=0.5, z0=1.4) 
  #vegetation area index, sum of leaf area density values for 
  #all horizontal slices assessed in previous line
  VAI.AOP <- sum(LADen$lad, na.rm=TRUE) 
  
  #VCI
  #vertical complexity index, fixed normalization of entropy 
  #metric calculated above
  #set zmax comofortably above maximum canopy height
  #by = 1 assesses the metric based on 1 m horizontal slices in 
  #the canopy
  VCI.AOP <- VCI(Zs, by = 1, zmax=80) #shouldn't really have trees over 80 m
  
  
  #OUTPUT CALCULATED METRICS INTO A TABLE
  #creates a dataframe row, out.plot, containing plot descriptors 
  #and calculated metrics
  structural_diversity <- rbind(structural_diversity,data.frame(matrix(c(as.numeric(i), site, year_mo, x, y, mean.max.canopy.ht, max.canopy.ht, 
                                                                         rumple, deepgap.fraction,
                                                                         cover.fraction,top.rugosity, vert.sd, vertCV,
                                                                         sd.sd, entro,GFP.AOP, VAI.AOP, VCI.AOP, Q25, Q50),
                                                                       ncol = 20)) )
    
  
    print(as.character(site))
  

}


colnames(structural_diversity) <- 
  c("raw_id","site", "year_mo", "easting", "northing", "mean.max.canopy.ht.aop",
    "max.canopy.ht.aop", "rumple.aop",
    "deepgap.fraction.aop","cover.fraction.aop", 
    "top.rugosity.aop", "vert.sd.aop", "vertCV.aop", "sd.sd.aop", 
    "entropy.aop", "GFP.AOP.aop", "VAI.AOP.aop", "VCI.AOP.aop", "q25.aop",
    "q50.aop") 

structural_diversity$new_ID = paste(structural_diversity$site,structural_diversity$year_mo,structural_diversity$easting,structural_diversity$northing)

# reading Mukund's old data frame to merge other plot level information with the structural metrics derived from lidar data at plot level

plot_data_table = read.csv("E:/NEON_Data/Structural_diversity/plot_lidar/plot_data_table.csv")

plot_data_table$new_ID = paste(plot_data_table$siteID,plot_data_table$monthyear,plot_data_table$easting,plot_data_table$northing)

# find matching rows between two data layers based on site, data collection date and the easting and northing information
match_Id = data.frame(match(structural_diversity$new_ID,plot_data_table$new_ID)) 

pull_Col = data.frame(plot_data_table[match_Id$match.structural_diversity.new_ID..plot_data_table.new_ID.,c(1:8)]) 

structural_diversity = cbind(pull_Col,structural_diversity[,c(6:20)])
write.csv(structural_diversity, "./StructuralDiversity_Liz03152023.csv")
write.csv(structural_diversity,"E:/NEON_Data/Structural_diversity/lidar_structural_metrics.csv")
