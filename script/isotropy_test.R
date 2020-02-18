rm(list=ls())
library(sp)           ## Data management
library(lattice)      ## Data management
library(geoR)         ## Geostatistics
library(gstat)
library(spTest)
library(rgdal)

load("Station_data.RData")
stat.data<-list(Station.data[,c('Longitude','Latitude')],Station.data[,'N.Trips'])
names(stat.data)<-c("coords","data")
#Converting coordinates to UTM reference
cord.dec = SpatialPoints(stat.data$coords, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32748"))
stat.data$coords<-cord.UTM@coords/1e+05
rm(cord.dec,cord.UTM)
#Variogram different directions and considering mean 1stºpolynomial(trend=1st) on lat and long->DB 
x11()
par(mfrow=c(2,2))
for(i in c(0,45,90,135)){res1.v <- variog(stat.data, trend="1st",
                                          uvec=seq(0,0.04,by=0.005),direction = i*2*pi/360) 
plot(res1.v,type="b",main = expression(paste("Directional Variogram, Angle = ",i,  degree)))
}
#Variogram in the same plot
x11()
res2.v <- variog4(stat.data, trend="1st",
                  uvec=seq(0,0.04,by=0.005))#,max.dist = 0.02) 
plot(res2.v,type="b")

##MaityTest
mydata<-cbind(stat.data$coords[,"Longitude"],stat.data$coords[,"Latitude"],stat.data$data)
myA <- rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
tol <- 1e-04
my.xlims <- range(stat.data$coords[,"Longitude"])+c(-tol,tol)
my.ylims <- range(stat.data$coords[,"Latitude"])+c(-tol,tol)
xlen <- my.xlims[2] - my.xlims[1]
ylen <- my.ylims[2] - my.ylims[1]
my.grid.spacing <- c(xlen / 16, ylen / 16)
xgrid <- seq(my.xlims[1], my.xlims[2], by = my.grid.spacing[1])
ygrid <- seq(my.ylims[1], my.ylims[2], by = my.grid.spacing[2])
my.blockdims <- c(4,4)
my.nBoot <- 150#In practice, use nBoot > 50
distlag=0.015#1.5 Chilometri
mylags <- rbind(c(distlag,0), c(0, distlag), c(sqrt((distlag^2)/2),sqrt((distlag^2)/2)), c(-sqrt((distlag^2)/2),sqrt((distlag^2)/2)))
tr.maity <- MaityTest(spdata = mydata, lagmat = mylags, A = myA,
                       df = 2, xlims = my.xlims, ylims = my.ylims,
                      grid.spacing = my.grid.spacing, block.dims = my.blockdims, nBoot = my.nBoot)
tr.maity$p.value
#Ho provato per distanze più grandi a 1.5Km non c'è isotropia

