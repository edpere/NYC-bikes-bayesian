###### STATIONS DATA (nov 2017) ######
## Loading stations data and pre-processing

trips_data = read.csv(file.path("data", "201711-citibike-tripdata.csv")) #750 stations
#tripsdata=read.csv("201306-citibike-tripdata.csv") #circa 300 stazioni

start.stations = unique(trips_data$start.station.id)

## In general, different stations can have exactly the same lat/long or there can be inconsistency in the data
## For this reason using 'unique' on lat/long is less reliable
## With 'match' we use the index of the first occurence of the station to retrieve the coordinates and avoid this problem
indices.first.occ = match(start.stations, trips_data$start.station.id)
latitudes = trips_data$start.station.latitude[indices.first.occ]
longitudes = trips_data$start.station.longitude[indices.first.occ]

## Combine stations ID with corresponding latitute and longitudes
start.stations = cbind(start.stations, latitudes, longitudes)
colnames(start.stations)[1] <- "stations.id"
rm(latitudes, longitudes, indices.first.occ)

## Order stations by ascending ID 
start.stations = start.stations[order(start.stations[,1]),]

## Count the number of trips for each station
## This count is going to be the response of the models
n.trips = rep(0, length(start.stations[,1]))
j=1
for(id in start.stations[,1]) {
  n.trips[j] = sum(trips_data$start.station.id==id)
  j=j+1
}

rm(j,id)


## Initialize a dataframe with station ID, number of trips, latitude and longitude
stations = cbind(start.stations[,1], n.trips, start.stations[,2:3])
stations = as.data.frame(stations)
colnames(stations)[1] = "station.ID"
colnames(stations)[3] = "latitude"
colnames(stations)[4] = "longitude"
rm(start.stations, n.trips)
#rm(trips_data)



##################################
### EXPLORATORY DATA ANALYSIS ###
##################################

##### SPATIAL OBJECTS + PLOTS #####
library(sf)
library(tidyr)
library(raster)
library(rgdal)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(rgeos)
library(dplyr)
library(geosphere)
library(ape)
library(spdep) 

#library(ggspatial)

#sp.stations = st_as_sf(stations, coords=c("longitude","latitude"), crs=4326)


###### ONLY MANHATTAN STATIONS (optional) #####

# Manual selection of perimeter of coordinates. 
# It should work even with a different dataset

# This first box encloses all of Manhattan (but also a lot more)
long_max=-74.02
long_min=-73.92
lat_max=40.82
lat_min=40.70

longmax1=-74.02
longmin1=-74.00
latmin1=40.70

longmax2=-74.00
longmin2=-73.97
latmin2=40.71

longmax3=-73.97
longmin3=-73.96
latmin3=40.75

longmax4=-73.96
longmin4=-73.955
latmin4=40.755

longmax5=-73.955
longmin5=-73.945
latmin5=40.76

longmax6=-73.945
longmin6=-73.94
latmin6=40.77

longmax7=-73.94
latmin7=40.785

set1=subset(stations,stations$latitude<lat_max & stations$latitude>latmin1 & stations$longitude<longmin1 & stations$longitude>longmax1)
set2=subset(stations,stations$latitude<lat_max & stations$latitude>latmin2 & stations$longitude<longmin2 & stations$longitude>longmax2)
set3=subset(stations,stations$latitude<lat_max & stations$latitude>latmin3 & stations$longitude<longmin3 & stations$longitude>longmax3)
set4=subset(stations,stations$latitude<lat_max & stations$latitude>latmin4 & stations$longitude<longmin4 & stations$longitude>longmax4)
set5=subset(stations,stations$latitude<lat_max & stations$latitude>latmin5 & stations$longitude<longmin5 & stations$longitude>longmax5)
set6=subset(stations,stations$latitude<lat_max & stations$latitude>latmin6 & stations$longitude<longmin6 & stations$longitude>longmax6)
set7=subset(stations,stations$latitude<lat_max & stations$latitude>latmin7 & stations$longitude<long_min & stations$longitude>longmax7)

stations.subset = rbind(set1,set2,set3,set4,set5,set6,set7)
sp.stations.subset = st_as_sf(stations.subset, coords=c("longitude","latitude"), crs=4326)

rm(set1,set2,set3,set4,set5,set6,set7)
rm(longmax1,longmax2,longmax3,longmax4,longmax5,longmax6,longmax7)
rm(longmin1,longmin2,longmin3,longmin4,longmin5,longmin6)
rm(latmin1,latmin2,latmin3,latmin4,latmin5,latmin6,latmin7)

# - - - - 

# NEW YORK MAP
sf.ny = readOGR(dsn = file.path(".", "maps", "borough"), layer = "borough")
#plot(sf.ny)
ny = tidy(sf.ny)
rm(sf.ny)

# BIKE LANES
lanes.ny = readOGR(dsn = file.path(".", "maps", "lanes"), "lanes")
bikelanes.ny = tidy(lanes.ny)
rm(lanes.ny)

# SUBWAY
sub.ny = readOGR(dsn = file.path(".", "maps", "subway"), layer = "subway")
subway.ny = tidy(sub.ny)
rm(sub.ny)

#----  THIS IS JUST FOR THE MAP (optional) BETTER VISUALIZATION ----
# By removing the maximum observation (in this dataset is over 15000 while
# no other observation is higher than 9050), we can have a better visualization
# of demand distribution across NYC.
# Once past the graph, however we can keep it in the dataset (it's Grand Central Station)

#stations.subset = setdiff(stations.subset, stations.subset[which.max(stations.subset$n.trips),])
stations.subset = setdiff(stations.subset, stations.subset[stations.subset$n.trips>10000,])

##### MAP OF STATIONS ##### 
x11()
sp=ggplot() +
  geom_path(data = ny, aes(x = long, y = lat, group = group)) +
  #geom_path(data = bikelanes.ny, aes(x = long, y = lat, group = group)) +
  #geom_path(data = subway.ny, aes(x = long, y = lat, group = group)) +
  #geom_point(data = stations, aes(x = longitude, y = latitude, color=n.trips))+
  geom_point(data = stations.subset, aes(x = longitude, y = latitude,color=n.trips))+
  #geom_text(data = stations.subset,aes(x=longitude, y=latitude,label=n.trips), hjust=0, vjust=1)+
  #geom_point(data = population, aes(x= longitude, y=latitude))+
  coord_sf(xlim=c(long_min, long_max), ylim=c(lat_min, lat_max))
mid = mean(stations.subset$n.trips)
sp + scale_colour_gradient2(midpoint=mid, low="green", mid="red", high="blue")

rm(mid, sp)


###### BASIC STATIONARITY CHECK ######
#---- I have considered only Manhattan in this part

# scatterplots: N.trips vs latitude and vs longitude
x11()
plot(stations.subset$n.trips, stations.subset$latitude, ylim=c(40.70,40.82), xlab="N.Trips", ylab="Latitude", main="N.trips vs Latitude")
x11()
plot(stations.subset$longitude, stations.subset$n.trips, xlim=c(-74.02,-73.92), xlab="Longitude", ylab="N.Trips", main="N.trips vs Longitude")

# Mean N.trips on smaller areas (intervals by latitude and longitude)
# THESE LINES MUST BE RUN STRICTLY IN THIS ORDER
x.lat = seq(lat_min, lat_max, by=0.005)
x.long = seq(long_max, long_min, by=0.005)

mean.lat = c()
for (i in 1:length(x.lat)){
  stsub = subset(stations.subset, stations.subset$latitude<x.lat[i+1] & stations.subset$latitude>x.lat[i])
  mean.lat[i] = mean(stsub$n.trips)
}
mean.long = c()
for (i in 1:length(x.long)){
  stsub = subset(stations.subset, stations.subset$longitude<x.long[i+1] & stations.subset$longitude>x.long[i])
  mean.long[i] = mean(stsub$n.trips)
}

c.lat = which(!is.na(mean.lat))
c.lat = x.lat[c.lat]
c.long = which(!is.na(mean.long))
c.long = x.long[c.long]
mean.lat = mean.lat[!is.na(mean.lat)]
mean.long = mean.long[!is.na(mean.long)]

# Mean values of n.trips wrt latitude
x11()
plot(mean.lat, c.lat, xlab="Mean", ylab="Latitude", main="Mean vs Latitude", col="red", pch=10)

# Mean values of n.trips wrt longitude
x11()
plot(c.long, mean.long, xlab="Longitude", ylab="Mean", main="Mean vs Longitude", col="red", pch=10)

rm(c.lat, c.long, i, mean.lat, mean.long, stsub, x.lat, x.long)

# There is no stationarity since the mean is not constant over all the map
# It is definitely higher in the center of Manhattan, and then decreases towards the
# outskirts. It seems to have less variance on smaller areas.
#- - - - - - - 

rm(bikelanes.ny, ny, subway.ny)




##### MORAN & GEARY AUTOCORRELATION TESTS #####

w = as.matrix(distm(cbind(stations.subset$longitude, stations.subset$latitude), fun = distGeo))
w.inv = 1/w
diag(w.inv) <- 0
Moran.I(stations.subset$n.trips, w.inv)

lw <- mat2listw(w.inv)
lwW <- nb2listw(lw$neighbours, glist=lw$weights, style="W")
moran.test(stations.subset$n.trips, lwW, alternative="two.sided")

stations.subset.scaled = scale(stations.subset$n.trips)[,1]
x11()
moran.plot(stations.subset.scaled, lw, lwd=1, xlab="Bike Demand", ylab="Lagged Bike Demand", labels=FALSE)

#GEARY'S test
geary.test(stations.subset$n.trips, lwW, alternative="two.sided")

rm(w,lw,w.inv,lwW,stations.subset.scaled)

##############
## DATASETS ##
##############
#################################################
## Create datasets for training and prediction ##
#################################################

## Define the area of interest
## Consider only a susbet of the stations 
### QUI SELEZIONO VALORI MAX E MIN


# South Manhattan
# long_max=-74.02
# long_min=-73.97
# lat_max=40.74
# lat_min=40.70

# Center Manhattan
long_max=-74.00
long_min=-73.95
lat_max=40.78
lat_min=40.74

# North Manhattan
#long_max=-73.99
#long_min=-73.93
#lat_max=40.82
#lat_min=40.78



stations.subset = subset(stations, stations$latitude<lat_max & stations$latitude>lat_min & stations$longitude<long_min & stations$longitude>long_max)
sp.stations.subset = st_as_sf(stations.subset, coords=c("longitude","latitude"), crs=4326)


################
## COVARIATES ##
################



##### BIKE LANES (Subset) #####
# 2) Seleziono le bike lanes nell'area che mi interessa
# estrazione coordinate delle bike lanes
# circa 9000 linestring (vedi unique(coordinates.subset$L2)) 

bike_lanes_data = read.csv(file.path("data", "nyc_bike_routes_2017.csv"))
sp.lanes = st_as_sfc(bike_lanes_data$the_geom, crs=4326) # create (multi)linestring geometry

coordinates = as.data.frame(st_coordinates(sp.lanes))

## Only lanes in the considered area
coordinates.subset = subset(coordinates, X<long_min & X>long_max & Y>lat_min & Y<lat_max)
rm(coordinates)
rownames(coordinates.subset) = c(1:length(coordinates.subset[,1]))

## Remove isolated points
tabb = as.data.frame(table(coordinates.subset$L2))  # L2 identify the original geometric object, so if only 1 it means that it is an isolated point
colnames(tabb)[1] = "id"
sum(tabb$Freq==1)  # Number of single points

ids_to_delete = as.vector(tabb$id[tabb$Freq==1])
rm(tabb)
# 43 indices to discard (with only one point, linestring is meaningless)
# Find the corresponing indices
ind = c()
k = 0
for (i in ids_to_delete){
  k=k+1
  ind[k] = which(coordinates.subset$L2==i)
}

rm(ids_to_delete)


valid_indices = setdiff(c(1:length(coordinates.subset[,1])), ind)
rm(k,i,ind)

## Points of the bike lanes
pts_lanes <- data.frame(
  x = coordinates.subset$X[valid_indices],
  y = coordinates.subset$Y[valid_indices],
  id = coordinates.subset$L2[valid_indices]) %>% 
  st_as_sf(coords = c("x","y")) %>% 
  st_set_crs(4326)
pts_lanes

rm(valid_indices, coordinates.subset)

sp.lanes.subset = pts_lanes %>%
      group_by(id) %>% 
      summarize() %>% 
      st_cast("LINESTRING")

sp.lanes.subset = st_as_sf(sp.lanes.subset, crs=4326)
sp.lanes.subset = sp.lanes.subset[,-1]

## Count how many lanes within a bound
bound = 500
lanes_count = c()
for (i in 1:length(stations.subset[,1])){
  lanes_count[i]= 0
  for(j in 1:length(as.data.frame(sp.lanes.subset)[,1])){
    dist =  st_distance(sp.stations.subset[i,3], sp.lanes.subset[1]$geometry[j], "Euclidean")
    if(as.vector(dist) <= bound) {
      lanes_count[i] = lanes_count[i] + 1
    }
  }
}

rm(i,j,bound,pts_lanes)
rm(bike_lanes_data, sp.lanes, dist)



##### POPULATION BY CENSUS BLOCKS (Subset)
library(readODS)

population = read.ods(file.path("data", "tracts.ods"))
population = as.data.frame(population)
population = population[2:2157, c(3,5,6)]
colnames(population)[1] = "Population"
colnames(population)[2] = "latitude"
colnames(population)[3] = "longitude"

# trasforma da "stringhe" a numeri (necessario per il corretto ordinamento)
population[,1] = as.numeric(population[,1])
population[,2] = as.numeric(population[,2])
population[,3] = as.numeric(population[,3])

# diseguaglianze inverse per longitude rispetto alle altre
block.subset = subset(population, longitude>long_max & longitude<long_min & latitude>lat_min & latitude<lat_max)

pts_block <- data.frame(
  x = block.subset$longitude,
  y = block.subset$latitude,
  p = block.subset$Population) %>% 
  st_as_sf(coords = c("x","y")) %>% 
  st_set_crs(4326)

block.points = pts_block %>%
  group_by(p) %>% 
  summarize() %>% 
  st_cast("POINT")

## Compute distance matrix
dist_station_block = matrix(nrow=length(stations.subset[,1]), ncol=length(block.points$geometry))
for (i in 1:length(stations.subset[,1])){
  for(j in 1:length(block.points$geometry)){
    dist_station_block[i,j] = st_distance(sp.stations.subset[i,3], block.points[1]$geometry[j], "Euclidean")
  }
}

nearest = c(ncol=length(stations.subset[,1]))
pop_block = c(ncol=length(stations.subset[,1]))
# indice del block più vicino a ciascuna stazione
for(i in 1:length(stations.subset[,1])){
  nearest[i] = which.min(dist_station_block[i,])
  pop_block[i] = block.points$p[nearest[i]]
}

rm(nearest, block.subset, dist_station_block, i, j, population, pts_block)




##### SUBWAY STATIONS POINTS (Subset) #####
subway.stat = read.csv(file.path("data", "Subway_stations.csv"))
sp.subway = st_as_sfc(subway.stat$the_geom, crs=4326) 

sub.data = as.data.frame(cbind(subway.stat$OBJECTID, st_coordinates(sp.subway)))

colnames(sub.data)[1] = "ID.subway.st"
colnames(sub.data)[2] = "longitude"
colnames(sub.data)[3] = "latitude"

subway.stat.subset = subset(sub.data, longitude<long_min & longitude>long_max & latitude>lat_min & latitude<lat_max)

pts_subway <- data.frame(
  x = subway.stat.subset$longitude,
  y = subway.stat.subset$latitude,
  id = subway.stat.subset$ID.subway.st) %>% 
  st_as_sf(coords = c("x","y")) %>% 
  st_set_crs(4326)

subway.points = pts_subway %>%
  group_by() %>% 
  summarize() %>% 
  st_cast("POINT")

## Create distance matrix
dist_station_subway = matrix(nrow=length(stations.subset[,1]), ncol=length(subway.points$geometry))
for (i in 1:length(stations.subset[,1])){
  for(j in 1:length(subway.points$geometry)){
    dist_station_subway[i,j] = st_distance(sp.stations.subset[i,3], subway.points[1]$geometry[j], "Euclidean")
  }
}

# Calcolo della distanza minima da una stazione della metro nel subset considerato
# (per ogni stazione del bike sharing)
subway_dist = rep(0, length(stations.subset[,1]))
for (i in 1:length(stations.subset[,1])){
  minim = min(dist_station_subway[i,])
  subway_dist[i] = minim
}

rm(sp.subway, dist_station_subway, i, j, minim, pts_subway)
rm(sub.data, subway.stat.subset, subway.stat)




##### PROXIMITY SCORE #####
## This covariate takes into account the effect of the stations nearby

## Create distance matrix (station to station)
dist_station_station = matrix(nrow=length(stations.subset[,1]), ncol=length(stations.subset[,1]))
for (i in 1:length(stations.subset[,1])) {
  for (j in i:length(stations.subset[,1])) {
    dist = st_distance(sp.stations.subset[i,3], sp.stations.subset[j,3])
    dist_station_station[i,j] = dist
    dist_station_station[j,i] = dist
  }
}

rm(dist, i, j)

## DIFFERENT CONSTRUCTION OF PROXIMITY MATRIX (useful for non-scaled data) ##
# With non-scaled data, values predicted by the model with proximity=sum(1/distances)
# do not make sense
# This proximity is similar the one used for W in CARBayes areal model

d = dist_station_station/1000
c = 1/d^2
c[is.infinite(c)] = 0
proximity_score = rowSums(c)

# The scores returned are between 70 and 334, instead of being very small values like 0.01
# that may drive prediction towards strange or negative values


rm(dist_station_station,c,d)


#### LANDMARKS SCORE ####
landmarks<-rbind(c(-73.966111,40.7825),#Central Park
                 c(-73.985708,40.75773),#Times Squate
                 c(-73.985428,40.748817),#Empire State
                 c(-73.979194,40.758611),#Rockefeller Center
                 c(-73.996111,40.705278),#Brooklyn Bridge
                 c(-74.005,40.7483),#High Line
                 c(-73.96367,40.77891),#Metropolitan Museum
                 c(-74.0135,40.713))#One World Trade Center
colnames(landmarks)<-c('Long','Lat')
rownames(landmarks)<-c('CP','TS','ES','RC','BB','HL','MM','OW')
library(spBayes)
#Station.data[,"Landmarks"]<-apply(iDist(Station.data[,c('Longitude','Latitude')],landmarks),1,min)
#save(Station.data, file = "Station_data_raw.RData")

# DIFFERENT LANDMARKS SCORE: use GREAT CIRCLE DISTANCE (in m)
# (useful for prediction with non-scaled data)
landmarks = as.data.frame(landmarks)
# sp.station.data = st_as_sf(Station.data, coords=c("Longitude","Latitude"), crs=4326)
sp.landmarks = st_as_sf(landmarks, coords=c("Long","Lat"), crs=4326)

dist.landmarks = matrix(nrow = length(stations.subset[,1]), ncol = length(landmarks[,1]))
for (i in 1:length(stations.subset[,1])){
  for (j in 1:length(landmarks[,1])){
    dist.landmarks[i,j] = st_distance(sp.stations.subset$geometry[i], sp.landmarks$geometry[j])
  }
}
# Score: distance to closest landmark
Landmarks = as.data.frame(apply(dist.landmarks,1,min))

rm(i,j,dist.landmarks,landmarks)


##### AGGREGAZIONE DATI: BLOCK POPULATION + LANES + SUBWAY STATIONS + PROXIMITY + LANDMARKS  #####
Station.data = data.frame(stations.subset$station.ID,
                          stations.subset$n.trips,
                          stations.subset$longitude,
                          stations.subset$latitude,
                          pop_block,
                          lanes_count,
                          subway_dist,
                          proximity_score,
                          Landmarks
)
colnames(Station.data)[1]="ID"
colnames(Station.data)[2]="N.Trips"
colnames(Station.data)[3]="Longitude"
colnames(Station.data)[4]="Latitude"
colnames(Station.data)[5]="Block.population"
colnames(Station.data)[6]="Lane.count"
colnames(Station.data)[7]="Dist.metro"
colnames(Station.data)[8]="Proximity.score"
colnames(Station.data)[9]="Landmarks"

rm(sp.station.data,sp.stations.subset)
rm(pop_block,proximity_score,subway_dist,lanes_count,Landmarks)

###### NORMALIZE DATA ######

## [0,1] normalization
normalize01 <- function(data) {
  (data - min(data))/(max(data) - min(data))
}
## zscore normalization
zscore <- function(data) {
  (data - mean(data))/sd(data)
}

norm_type = "zscore"

if (norm_type == '01') {
  for(col in setdiff(colnames(Station.data), c('ID', 'Longitude', 'Latitude'))) {
    Station.data[[col]] <- normalize01(Station.data[[col]])
  }
} else if (norm_type == 'zscore') {
  for(col in setdiff(colnames(Station.data), c('ID', 'Longitude', 'Latitude'))) {
    Station.data[[col]] <- zscore(Station.data[[col]])
  }
}

rm(i, norm_type, normalize01, zscore, col)
rm(stations.subset,stations)

#save(Station.data, file = "Station_data.RData")


## Correlation Matrix of covariates ##
cor(Station.data[,5:ncol(Station.data)])

library(corrplot)
corrplot(cor(Station.data[,5:ncol(Station.data)]), type='upper', method='circle', outline=TRUE, tl.col="black")




#### PREDICTION GRID ####
## Create a prediction grid and associated covariates
library(geoR)
pred.coords <- pred_grid(c(long_max, long_min), c(lat_min, lat_max), by=0.005)
colnames(pred.coords) <- c("longitude", "latitude")
sp.pred.coords <- st_as_sf(pred.coords, coords=c("longitude","latitude"), crs=4326)

## Population
dist_matrix = matrix(nrow=length(sp.pred.coords$geometry), ncol=length(block.points$geometry))
for (i in 1:length(sp.pred.coords$geometry)){
  for(j in 1:length(block.points$geometry)){
    dist_matrix[i,j] = st_distance(sp.pred.coords[i,1], block.points[1]$geometry[j], "Euclidean")
  }
}

nearest = c(ncol=length(pred.coords[,1]))
pop_block = c(ncol=length(pred.coords[,1]))
for(i in 1:length(pred.coords[,1])){
  nearest[i] = which.min(dist_matrix[i,])
  pop_block[i] = block.points$p[nearest[i]]
}

rm(nearest, dist_matrix, i, j)

## Bike Lanes
bound = 500
lanes_count = c()
for (i in 1:length(sp.pred.coords$geometry)){
  lanes_count[i] = 0
  for(j in 1:length(sp.lanes.subset$geometry)){
    dist =  st_distance(sp.pred.coords[i,1], sp.lanes.subset[1]$geometry[j], "Euclidean")
    if(as.vector(dist) <= bound) {
      lanes_count[i] = lanes_count[i] + 1
    }
  }
}

rm(i, j, dist, bound)

## Subway
dist_matrix = matrix(nrow=length(sp.pred.coords$geometry), ncol=length(subway.points$geometry))
for (i in 1:length(sp.pred.coords$geometry)){
  for(j in 1:length(subway.points$geometry)){
    dist_matrix[i,j] = st_distance(sp.pred.coords[i,1], subway.points[1]$geometry[j], "Euclidean")
  }
}

subway_dist = rep(0, length(pred.coords[,1]))
for (i in 1:length(pred.coords[,1])){
  minim = min(dist_matrix[i,])
  subway_dist[i] = minim
}

rm(i, j, minim, dist_matrix)

## Landmarks
dist_matrix = matrix(nrow = length(sp.pred.coords$geometry), ncol = length(sp.landmarks$geometry))
for (i in 1:length(sp.pred.coords$geometry)){
  for (j in 1:length(sp.landmarks$geometry)){
    dist_matrix[i,j] = st_distance(sp.pred.coords$geometry[i], sp.landmarks$geometry[j])
  }
}
# Score: distance to closest landmark
landmarks.score = as.data.frame(apply(dist_matrix,1,min))

rm(i,j, dist_matrix)

## Aggregation
Grid.data = data.frame(pred.coords$longitude,
                       pred.coords$latitude,
                       pop_block,
                       lanes_count,
                       subway_dist,
                       landmarks.score
                       )

colnames(Grid.data) <- c("Longitude", "Latitude", "Block.population", "Lane.count", "Dist.metro", "Landmarks")


rm(pred.coords, sp.pred.coords, pop_block, lanes_count, subway_dist, landmarks.score)
rm(block.points, sp.landmarks, sp.lanes.subset, subway.points)



ggplot() +
  geom_path(data = ny, aes(x = long, y = lat, group = group)) +
  geom_point(data = Grid.data[,1:2], aes(x = Longitude, y = Latitude), color="blue")+
  geom_text(data = Grid.data,aes(x=Longitude, y=Latitude,label=Block.population), hjust=0, vjust=1)+
  coord_sf(xlim=c(long_min, long_max), ylim=c(lat_min, lat_max))

rm(ny)

save(Grid.data, lat_max, lat_min, long_max, long_min, file = "Prediction_Grid.RData")

