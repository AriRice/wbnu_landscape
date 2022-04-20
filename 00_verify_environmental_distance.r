
library(rJava)
library(rgeos)
library(sp)
library(maptools)
library(raster)
library(dismo)
library(rgdal)
library(reshape2)

# raster directories at 2.5 arc minutes (worldclim 1.4 layers)
worldclim_directory <- "bioclim"

# raster directories at 2.5 arc minutes (worldclim 2.1 layers)
worldclim_directory2 <- "bioclim_new"

# read in layers (original worldclim v1 layers from sitta RAD paper)
worldclim_layers <- stack(list.files(worldclim_directory, full.names=T))
# read in layers (worldclim v2.1 layers)
worldclim_layers2 <- stack(list.files(worldclim_directory2, full.names=T))

# subset for modelling based on original sitta paper
worldclim_modelling <- stack(list.files(worldclim_directory, full.names=T)[list.files(worldclim_directory) %in% paste0(c("bio_1", "bio_2", "bio_4", "bio_5", "bio_6", "bio_9", "bio_12", "bio_15", "bio_17", "bio_18", "bio_19"), ".asc")])
# subset for modelling based on original sitta paper
worldclim_modelling2 <- stack(list.files(worldclim_directory2, full.names=T)[list.files(worldclim_directory2) %in% paste0("wc2.1_2.5m_", c("bio_1", "bio_2", "bio_4", "bio_5", "bio_6", "bio_9", "bio_12", "bio_15", "bio_17", "bio_18", "bio_19"), ".tif")])

# clip layers to North America
# max lat = -50, min lat = -165, max long = 70, min long = 15
clip_namerica <- cbind(c(-50, -50, -165, -165), c(15, 70, 70, 15))
clip_namerica <- SpatialPolygons(list(Polygons(list(Polygon(clip_namerica)), 1)))
worldclim_layers <- crop(worldclim_layers, clip_namerica)
worldclim_modelling <- crop(worldclim_modelling, clip_namerica)
worldclim_layers2 <- crop(worldclim_layers2, clip_namerica)
worldclim_modelling2 <- crop(worldclim_modelling2, clip_namerica)

# file of occurrence data
input_file <- read.csv("Model_points/Occurrence_subset_2km_train.csv")
coordinates(input_file) <- c("Long", "Lat")
proj4string(input_file) <- CRS("+init=epsg:4326")

# set up training region for each species 
# here we will use 250 km (250 *1000 m)
model_buffer <- 500 * 1000

# create a projection of the points in CRS 
points_projection <- spTransform(input_file, CRS=CRS("+init=epsg:3857"))
  
# buffer the points by the model buffer and transform back to WGS84
clipping_mask <- spTransform(gBuffer(points_projection, width=model_buffer), CRS=CRS("+init=epsg:4326"))
  
# create a training raster stack using the 500km buffer crop
worldclim_modelling <- crop(worldclim_modelling, clipping_mask)
worldclim_modelling <- mask(worldclim_modelling, clipping_mask)
worldclim_modelling2 <- crop(worldclim_modelling2, clipping_mask)
worldclim_modelling2 <- mask(worldclim_modelling2, clipping_mask)


# extract bioclim values for observed points and 10000 background points 
bioclim_values <- extract(worldclim_modelling, input_file)
random_bioclim <- sampleRandom(worldclim_modelling, 10000)
bioclim_values <- data.frame(rbind(cbind(presence=1, bioclim_values), cbind(presence=0, random_bioclim)))
bioclim_values2 <- extract(worldclim_modelling2, input_file)
random_bioclim2 <- sampleRandom(worldclim_modelling2, 10000)
bioclim_values2 <- data.frame(rbind(cbind(presence=1, bioclim_values2), cbind(presence=0, random_bioclim2)))

# perform maxent model
niche_model_R <- maxent(x=bioclim_values[ ,2:ncol(bioclim_values)], p=bioclim_values[ ,1], path="./maxent_output")
projection_contemporary <- predict(niche_model_R, worldclim_modelling)
niche_model_R2 <- maxent(x=bioclim_values2[ ,2:ncol(bioclim_values2)], p=bioclim_values2[ ,1], path="./maxent_output2")
projection_contemporary2 <- predict(niche_model_R2, worldclim_modelling2)


# model reclassify function
# model is already read in as a raster
# threshold is numeric between 0 and 1
reclass_ENM <- function(model, threshold) {
  require(raster)
  y <- model
  y[y >= threshold] <- 1
  y[y < threshold] <- NA
  return(y)
}

# reclassify inclusive of XX% of training points to get reclassified model similar
# to that shown in original paper
train_values <- extract(projection_contemporary, input_file)
threshold <- sort(train_values)[floor(length(train_values) - length(train_values) * 0.85)]
reclass_contemporary <- reclass_ENM(projection_contemporary, threshold)
train_values2 <- extract(projection_contemporary2, input_file)
threshold2 <- sort(train_values2)[floor(length(train_values2) - length(train_values2) * 0.85)]
reclass_contemporary2 <- reclass_ENM(projection_contemporary2, threshold2)

# change reclassified niche model to polygon feature
contemporary_polygons <- rasterToPolygons(reclass_contemporary, dissolve=T)
contemporary_polygons2 <- rasterToPolygons(reclass_contemporary2, dissolve=T)

# extract 5000 random points from reclassified niche model
niche_points <- spsample(contemporary_polygons, 5000, type="random")
niche_coords <- niche_points@coords
colnames(niche_coords) <- c("Long", "Lat")
niche_points2 <- spsample(contemporary_polygons2, 5000, type="random")
niche_coords2 <- niche_points2@coords
colnames(niche_coords2) <- c("Long", "Lat")
# add the genetic sampling localities (N = 9) to the niche_points object
localities <- read.csv("sitta_localities.csv", stringsAsFactors=F)
localities2 <- localities[,2:3]
niche_coords <- rbind(localities2, niche_coords)
coordinates(niche_coords) <- c("Long", "Lat")
proj4string(niche_coords) <- CRS("+init=epsg:4326")
niche_coords2 <- rbind(localities2, niche_coords2)
coordinates(niche_coords2) <- c("Long", "Lat")
proj4string(niche_coords2) <- CRS("+init=epsg:4326")

# just the real localities worldclim values
coordinates(localities2) <- c("Long", "Lat")
proj4string(localities2) <- CRS("+init=epsg:4326")
localities_wc <- extract(worldclim_layers, localities2)
localities_wc2 <- extract(worldclim_layers2, localities2)
cbind(localities, localities_wc[,c(6,15)], localities_wc2[,c(6,15)])

# extract bioclim values for all points
niche_coords_bc <- extract(worldclim_layers, niche_coords)
niche_coords_prcomp <- prcomp(niche_coords_bc, center=F, scale.=T)
niche_coords_PCs <- niche_coords_prcomp$x[,1:4]
localities_PC <- cbind(localities, niche_coords_PCs[1:9,])
niche_coords_bc2 <- extract(worldclim_layers2, niche_coords2)
niche_coords_prcomp2 <- prcomp(niche_coords_bc2, center=F, scale.=T)
niche_coords_PCs2 <- niche_coords_prcomp2$x[,1:4]
localities_PC2 <- cbind(localities, niche_coords_PCs2[1:9,])

# distances between localities
euc_distances <- dist(localities_PC[,4:ncol(localities_PC)])
euc_distances <- melt(as.matrix(euc_distances), varnames = c("Pop1", "Pop2"))
euc_distances <- euc_distances[euc_distances$Pop2 > euc_distances$Pop1,]
euc_distances2 <- dist(localities_PC2[,4:ncol(localities_PC2)])
euc_distances2 <- melt(as.matrix(euc_distances2), varnames = c("Pop1", "Pop2"))
euc_distances2 <- euc_distances2[euc_distances2$Pop2 > euc_distances2$Pop1,]

# change names of pops in euc_distances matrix
for(a in 1:nrow(localities)) {
  euc_distances$Pop1[euc_distances$Pop1 == a] <- localities$Pop1[a]
  euc_distances$Pop2[euc_distances$Pop2 == a] <- localities$Pop1[a]
  euc_distances2$Pop1[euc_distances2$Pop1 == a] <- localities$Pop1[a]
  euc_distances2$Pop2[euc_distances2$Pop2 == a] <- localities$Pop1[a]
}

# read in original distance matrix from the 2015 RAD paper
original_distances <- read.table("sitta_distances.txt", stringsAsFactors=F, header=T)

# loop through to extract the e-distances in the correct order of the 
# original_distances file
new_e_distances <- list()
new_e_distances2 <- list()
for(a in 1:nrow(original_distances)) {
  new_e_distances[[a]] <- euc_distances$value[ 
  (euc_distances$Pop1 == original_distances$Pop1[a] & 
    euc_distances$Pop2 == original_distances$Pop2[a]) | 
    (euc_distances$Pop1 == original_distances$Pop2[a] & 
       euc_distances$Pop2 == original_distances$Pop1[a])]
  new_e_distances2[[a]] <- euc_distances2$value[ 
    (euc_distances2$Pop1 == original_distances$Pop1[a] & 
       euc_distances2$Pop2 == original_distances$Pop2[a]) | 
      (euc_distances2$Pop1 == original_distances$Pop2[a] & 
         euc_distances2$Pop2 == original_distances$Pop1[a])]
}
new_e_distances <- unlist(new_e_distances)
new_e_distances2 <- unlist(new_e_distances2)
original_distances <- cbind(original_distances, new_e_distances, new_e_distances2)

write.table(original_distances, file="sitta_distances_new.txt", sep="\t", row.names=F, quote=F)
# save

save.image("e_dist_verified.RData")

load("e_dist_verified.RData")

par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
plot(original_distances$e_dist, original_distances$new_e_distances, ylim=c(0,3), xlim=c(0,3), xlab="Env. Distances (2015 Paper)", ylab="Env. Distances (Redone w/ Reproducible Code) (WC v1.4)", pch=19)
abline(lm(original_distances$new_e_distances ~ original_distances$e_dist))
summary(lm(original_distances$new_e_distances ~ original_distances$e_dist))

plot(original_distances$e_dist, original_distances$new_e_distances2, ylim=c(0,3), xlim=c(0,3), xlab="Env. Distances (2015 Paper)", ylab="Env. Distances (Redone w/ Reproducible Code) (WC v2.1)", pch=19)
abline(lm(original_distances$new_e_distances2 ~ original_distances$e_dist))
summary(lm(original_distances$new_e_distances2 ~ original_distances$e_dist))

plot(original_distances$new_e_distances, original_distances$new_e_distances2, ylim=c(0,3), xlim=c(0,3), xlab="Env. Distances (Redone w/ Reproducible Code) (WC v1.4)", ylab="Env. Distances (Redone w/ Reproducible Code) (WC v2.1)", pch=19)
abline(lm(original_distances$new_e_distances2 ~ original_distances$new_e_distances))
summary(lm(original_distances$new_e_distances2 ~ original_distances$new_e_distances))





