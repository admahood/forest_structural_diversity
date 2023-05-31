
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
library(geometry)
library(vegan)
library(ggdendro)
library(FD)
library(psych)

wd = "E:/NEON_Data/Structural_diversity/"
setwd(wd)

sp_tr = read.csv("E:/NEON_Data/Structural_diversity/StructuralDiversity_Liz03152023.csv")

plot_div = read.csv("E:/forest_structural_diversity/data/plot_site_summaries.csv")
plot_div$ID2 = paste0(plot_div$year,plot_div$plotID)

slect_plots = read.csv("E:/forest_structural_diversity/output/cover_by_plot.csv")
slect_plots$ID2 = paste0(slect_plots$year,slect_plots$plotID)

match_id1 = match(slect_plots$ID2,plot_div$ID2)

slect_plots = cbind(slect_plots,plot_div[match_id1,])



slect_plots$new_ID = paste0(slect_plots$siteID," ",slect_plots$monthyear," ",slect_plots$easting," ",slect_plots$northing)
slect_plots$new_ID<-gsub(" ","",as.character(slect_plots$new_ID)) # remove spaces between new_ID attributes


# sp_tr$new_ID = paste0(sp_tr$site,sp_tr$year_mo, sp_tr$easting,sp_tr$northing)
# 
# match_id2 = match(slect_plots$new_ID,sp_tr$new_ID)
# match_id3 = match(sp_tr$new_ID,slect_plots$new_ID)
# 
# 
# slect_plots = cbind(slect_plots,sp_tr[match_id2,])

#write.csv(slect_plots, "slect_plots.csv")

sp_tr = sp_tr[,c(1:21,31:36,42:71)]


sp_tr = sp_tr[,-c(49)] # removed Lidar based VCI column as it has a lots of 'NA's
sp_tr = na.omit(sp_tr)

sp_tr2 = sp_tr[,c(21,37:48)] # extract only the lidar based diversity and  data and invasion status

# tr_cat = read.csv("test_tr_cat.csv")

#Normalized the data before PCA analysis


#normalization function

normalized<-function(y) {
  
  x<-y[!is.na(y)]
  
  x<-(x - min(x)) / (max(x) - min(x))
  
  y[!is.na(y)]<-x
  
  return(y)
}

#Normalize lidar metrics
sp_tr2 = data.frame(apply(sp_tr2[,c(2:13)],2,normalized))
sp_tr2 = cbind(sp_tr$invaded,sp_tr2)
sp_tr2 = sp_tr2[-823,] # gave a strange value in PCA. So removed for now.


# Comm_ent = sp.to.fe(sp_tr2, tr_cat, fe_nm_type = "fe_rank", check_input = TRUE)
# 
# en = Comm_ent[["fe_nb_sp"]]


#####using vegan package

# decorana (sp_tr2)
dev.new(width = 3, height = 3)
pairs.panels(sp_tr2[,-1],
             gap = 0,
             bg = c("red", "blue")[sp_tr2$`sp_tr$invaded`],
             pch=21)


##PCA analysis
PCA_lidar <- prcomp(sp_tr2[,-1],
             center = TRUE,
             scale. = TRUE)
attributes(PCA_lidar)


# PCA_lidar = rda(sp_tr2[,-1])

head(summary(PCA_lidar))

dev.new(width = 3, height = 3)
pairs.panels(PCA_lidar$x,
             gap=0,
             bg = c("red", "blue")[sp_tr2$`sp_tr$invaded`],
             pch=21)



dev.new(width = 3, height = 3)
sitetype <- as.numeric (as.factor (sp_tr$invaded))
ordiplot (PCA_lidar$x, display = 'sites', type = 'n', main = "Plot scale PCA lidar(invaded vs not invaded)")
points (PCA_lidar$x, pch = sitetype, col = sitetype)

## the most easy method (if your Group0 isn't alphabetical order, this method can't be used.)
legend(x="bottomright", legend=levels(as.factor(sp_tr$invaded)), col=1:30, pch=1:30)
#  you can use seq.int(levels(MyMeta$Group)) instead of 1:7


### Bi_plots

library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)

dev.new(width = 3, height = 3)

g <- ggbiplot(PCA_lidar,
              obs.scale = 1,
              var.scale = 1,
              groups = sp_tr$invaded[-823],
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)


##Adding landcover class per each plot

landcover = data.frame()

file_list = dir("E:/NEON_Data/Structural_diversity/drive-download-20230531T160737Z-001/", pattern = "NEON_plots*", full.names = FALSE, ignore.case = TRUE)
#names_to_cloud = substr(file_list, 1, nchar(file_list)-8)
for (i in 1:length(file_list))
{
  print(i)
  lcd = read.csv(paste0("E:/NEON_Data/Structural_diversity/drive-download-20230531T160737Z-001/",file_list[i]))
  landcover = rbind(landcover,lcd[,c(2:9)])
}




#plot_cover = read.csv("E:/forest_structural_diversity/output/cover_by_plot.csv")

landcover$new_ID = paste0(landcover$site,landcover$year_mo,landcover$easting,landcover$northing)

sp_tr$newID = paste0(sp_tr$site,sp_tr$year,sp_tr$easting,sp_tr$northing)

match_id = match(landcover$new_ID,sp_tr$newID)

cover_to_sp_tr = landcover[match_id,]

landcover = cbind(landcover,sp_tr[match_id,])

#lidar_cover = cbind(sp_tr,cover_to_sp_tr)
# lidar_cover = na.omit((lidar_cover))

match_id2 = match(landcover$new_ID,slect_plots$new_ID)
landcover = cbind(landcover,slect_plots[match_id2,])


dev.new(width = 3, height = 3)
sitetype <- as.numeric (as.factor (lidar_cover$NLCD_Classes_Abb[-823]))
ordiplot (PCA_lidar$x, display = 'sites', type = 'n', main = "PCA of plot scale cover diversity by NLCD class")
points (PCA_lidar$x, pch = sitetype, col = sitetype)

## the most easy method (if your Group0 isn't alphabetical order, this method can't be used.)
legend(x="bottomright",  cex=0.8,legend=levels(as.factor(lidar_cover$NLCD_Classes_Abb)), col=1:20, pch=1:20)



######Dendrogram for invaded vs non-invaded classification using lidar


##Dendrogram with lidar data

dev.new(width = 3, height = 3)
if(require(rpart)){
  model <- rpart(lidar_cover$invaded ~ max.canopy.ht.aop + rumple.aop+
                   deepgaps.aop+deepgap.fraction.aop+cover.fraction.aop+ top.rugosity.aop + 
                   vert.sd.aop + sd.sd.aop + entropy.aop + GFP.AOP.aop, 
                 method = "class", data = sp_tr2)
  ddata <- dendro_data(model)
  ggplot() + 
    geom_segment(data = ddata$segments, 
                 aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = ddata$labels, 
              aes(x = x, y = y, label = label), size = 4, vjust = 0) +
    geom_text(data = ddata$leaf_labels, 
              aes(x = x, y = y, label = label), size = 4, vjust = 1) +
    theme_dendro()
}

##Dendrogram with field data
sp_tr3 = sp_tr[,c(4:32)]
  
dev.new(width = 3, height = 3)  
if(require(rpart)){
  model <- rpart(invaded ~ shannon_total + evenness_total + nspp_total, 
                 method = "class", data = sp_tr3)
  ddata <- dendro_data(model)
  ggplot() + 
    geom_segment(data = ddata$segments, 
                 aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = ddata$labels, 
              aes(x = x, y = y, label = label), size = 4, vjust = 0) +
    geom_text(data = ddata$leaf_labels, 
              aes(x = x, y = y, label = label), size = 4, vjust = 1) +
    theme_dendro()
}


##Dendrogram with PCA from lidar data
XX = data.frame(PCA_lidar$x)
XX = cbind(sp_tr[-823,30],XX)

dev.new(width = 3, height = 3)  
if(require(rpart)){
  model <- rpart(sp_tr[-823, 30] ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 +
                   PC7 + PC8 + PC9 + PC10 + PC11 + PC12 , 
                 method = "class", data = XX)
  ddata <- dendro_data(model)
  ggplot() + 
    geom_segment(data = ddata$segments, 
                 aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = ddata$labels, 
              aes(x = x, y = y, label = label), size = 4, vjust = 0) +
    geom_text(data = ddata$leaf_labels, 
              aes(x = x, y = y, label = label), size = 4, vjust = 1) +
    theme_dendro()
}


X2  =sp_tr[-823, -c(1:3,30,33,34,49)]
cor_all = cor(X2,method = "spearman")

cor_greater_20 = as.data.frame(apply(cor_all, 2, function(x) ifelse (abs(x) >=0.30,x,"NA")))

plot(cor_greater_20)










##########################

# Give the chart file a name.
png(file = "max_canopy_height.png")

# Plot the chart.
b = boxplot(max.canopy.ht.aop ~ NLCD_Classes_Abb, data = sp_tr, xlab = "NLCD class",
        ylab = "Max canopy height (m)", main = "Lidar Max Height ", notch = FALSE, 
        varwidth = TRUE, 
        col = as.factor(sp_tr$NLCD_Classes_Abb)) 

# Save the file.
dev.off()



dev.new(width=2, height=1)
lab_text = element_text(face = "bold", color = "black", size = 16)
lab_text_x = element_text(face = "bold", color = "black", size = 16,angle=45)
axis_text = element_text(face = "bold", color = "black", size = 16)
ggplot() + 
  geom_smooth(data = sp_tr, aes(x=cover_native, y=easting),method = "lm",
              formula = y ~ poly(x, 2), se = TRUE ,  linetype="dashed",
              color="red") +
  geom_smooth(data = sp_tr, aes(x=cover_exotic, y=easting),method = "lm",
              formula = y ~ poly(x, 2), se = TRUE ,  linetype="dashed",
              color="green") +
  # geom_smooth(data = newdata, aes(x=chrono, y=fit2),method = "lm",
  #             formula = y ~ poly(x, 2), se = TRUE ,  linetype="dashed",
  #             color="pink") +
  # geom_hline(yintercept=0, linetype="dashed", 
  #            color = "red", size=2)+
  theme(legend.position='right') +
  labs(y = " easting ", x = "cover_native" ,title = "easting with relative cover native" )+
  theme(axis.line.x = element_line(size = 0.7, colour = "black"),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = lab_text_x, axis.text.y = lab_text, axis.title = axis_text,
        plot.title = element_text(hjust = 0.1))


dev.new(width=2, height=1)
lab_text = element_text(face = "bold", color = "black", size = 16)
lab_text_x = element_text(face = "bold", color = "black", size = 16,angle=45)
axis_text = element_text(face = "bold", color = "black", size = 16)
ggplot() + 
  geom_smooth(data = sp_tr, aes(x=cover_native, y=deepgap.fraction.aop),method = "lm",
              formula = y ~ poly(x, 2), se = TRUE ,  linetype="dashed",
              color="red") +
  geom_smooth(data = sp_tr, aes(x=cover_exotic, y=deepgap.fraction.aop),method = "lm",
              formula = y ~ poly(x, 2), se = TRUE ,  linetype="dashed",
              color="green") +
  # geom_smooth(data = newdata, aes(x=chrono, y=fit2),method = "lm",
  #             formula = y ~ poly(x, 2), se = TRUE ,  linetype="dashed",
  #             color="pink") +
  # geom_hline(yintercept=0, linetype="dashed", 
  #            color = "red", size=2)+
  theme(legend.position='right') +
  labs(y = "deepgap.fraction.aop", x = "Percent cover" ,title = "Lidar deepgap.fraction.aop with percent veg. cover" )+
  theme(axis.line.x = element_line(size = 0.7, colour = "black"),
        axis.line.y = element_line(size = 0.7, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = lab_text_x, axis.text.y = lab_text, axis.title = axis_text,
        plot.title = element_text(hjust = 0.1))

