library(rJava)
library(rgeos)
library(sp)
library(maptools)
library(raster)
library(dismo)
library(rgdal)
library(reshape2)

countries <- readOGR("./ne_50m_admin_0_countries/ne_50m_admin_0_countries.shp")
Ethiopia <- countries[countries$NAME == "Ethiopia", ]

# raster directories at 2.5 arc minutes
worldclim_directory <- "bioclim"


# read in layers (original worldclim v1 layers from sitta RAD paper)
worldclim_layers <- stack(list.files(worldclim_directory, full.names=T))

# subset for modelling based on original sitta paper
worldclim_modelling <- stack(list.files(worldclim_directory, full.names=T)[list.files(worldclim_directory) %in% paste0(c("bio_1", "bio_2", "bio_4", "bio_5", "bio_6", "bio_9", "bio_12", "bio_15", "bio_17", "bio_18", "bio_19"), ".asc")])

# clip layers to North America
# max lat = -50, min lat = -165, max long = 70, min long = 15
clip_namerica <- cbind(c(-50, -50, -165, -165), c(15, 70, 70, 15))
clip_namerica <- SpatialPolygons(list(Polygons(list(Polygon(clip_namerica)), 1)))
worldclim_layers <- crop(worldclim_layers, clip_namerica)
worldclim_modelling <- crop(worldclim_modelling, clip_namerica)

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


# extract bioclim values for observed points and 10000 background points 
bioclim_values <- extract(worldclim_modelling, input_file)
random_bioclim <- sampleRandom(worldclim_modelling, 10000)
bioclim_values <- data.frame(rbind(cbind(presence=1, bioclim_values), cbind(presence=0, random_bioclim)))

# perform maxent model
niche_model_R <- maxent(x=bioclim_values[ ,2:ncol(bioclim_values)], p=bioclim_values[ ,1], path="./maxent_output")
projection_contemporary <- predict(niche_model_R, worldclim_modelling)


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
threshold <- sort(train_values)[floor(length(train_values) - length(train_values) * 0.9)]
reclass_contemporary <- reclass_ENM(projection_contemporary, threshold)

# change reclassified niche model to polygon feature
contemporary_polygons <- rasterToPolygons(reclass_contemporary, dissolve=T)

# extract 5000 random points from reclassified niche model
niche_points <- spsample(contemporary_polygons, 5000, type="random")
niche_coords <- niche_points@coords
colnames(niche_coords) <- c("Long", "Lat")
# add the genetic sampling localities (N = 9) to the niche_points object
localities <- read.csv("sitta_localities.csv", stringsAsFactors=F)
localities2 <- localities[,2:3]
niche_coords <- rbind(localities2, niche_coords)
coordinates(niche_coords) <- c("Long", "Lat")
proj4string(niche_coords) <- CRS("+init=epsg:4326")

# extract bioclim values for all points
niche_coords_bc <- extract(worldclim_layers, niche_coords)
niche_coords_prcomp <- prcomp(niche_coords_bc, center=F, scale.=T)
niche_coords_PCs <- niche_coords_prcomp$x[,1:4]
localities <- cbind(localities, niche_coords_PCs[1:9,])

# distances between localities
euc_distances <- dist(localities[,4:ncol(localities)])
euc_distances <- melt(as.matrix(euc_distances), varnames = c("Pop1", "Pop2"))
euc_distances <- euc_distances[euc_distances$Pop2 > euc_distances$Pop1,]

# change names of pops in euc_distances matrix
for(a in 1:nrow(localities)) {
  euc_distances$Pop1[euc_distances$Pop1 == a] <- localities$Pop1[a]
  euc_distances$Pop2[euc_distances$Pop2 == a] <- localities$Pop1[a]
}

# read in original distance matrix from the 2015 RAD paper
original_distances <- read.table("sitta_distances.txt", stringsAsFactors=F, header=T)

# loop through to extract the e-distances in the correct order of the 
# original_distances file
new_e_distances <- list()
for(a in 1:nrow(original_distances)) {
  new_e_distances[[a]] <- euc_distances$value[ 
  (euc_distances$Pop1 == original_distances$Pop1[a] & 
    euc_distances$Pop2 == original_distances$Pop2[a]) | 
    (euc_distances$Pop1 == original_distances$Pop2[a] & 
       euc_distances$Pop2 == original_distances$Pop1[a])]
}
new_e_distances <- unlist(new_e_distances)
original_distances <- cbind(original_distances, new_e_distances)
plot(original_distances$e_dist, original_distances$new_e_distances)
abline(lm(original_distances$new_e_distances ~ original_distances$e_dist))
summary(lm(original_distances$new_e_distances ~ original_distances$e_dist))

write.table(original_distances, file="sitta_distances_new.txt", sep="\t", row.names=F, quote=F)
# save

save.image("e_dist_verified.RData")
