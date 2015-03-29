#setwd("data")

# ref http://www.emast.org.au/observations/climate/
#install.packages("ncdf", type = "source", configure.args="--with-netcdf-include=/usr/include")
require(ncdf)

## Loading required package: ncdf

#install.packages("raster", dependencies = T)
require(raster)

## Loading required package: raster
## Loading required package: sp

# install.packages("rgdal")
require(rgdal)

# Loading required package: rgdal
# rgdal: version: 0.9-1, (SVN revision 518)
# Geospatial Data Abstraction Library extensions to R successfully loaded
# Loaded GDAL runtime: GDAL 1.9.2, released 2012/10/08
# Path to GDAL shared files: /usr/share/gdal
# Loaded PROJ.4 runtime: Rel. 4.8.0, 6 March 2012, [PJ_VERSION: 480]
# Path to PROJ.4 shared files: (autodetected)

if extracting for points shapefile
require(ggmap)
require(rgdal)
locn <- geocode("Prospect Reservoir NSW")
# this uses google maps API, better check this

## Treat data frame as spatial points
epsg <- make_EPSG()
shp <- SpatialPointsDataFrame(cbind(locn$lon,locn$lat),locn,
                              proj4string=CRS(epsg$prj4[epsg$code %in% '4283']))
str(shp)
shp@data
# #writeOGR(shp, 'test.shp', 'test', driver  = "ESRI Shapefile")
# 
# #shp <- readOGR(dsn="test.shp", layer='test')
# #plot(shp, add = T)
# 
# # a loop through days, see comment sections that print for debugging
# strt <-'2012-01-01'
# end <- '2012-01-04'
# dates <- seq(as.Date(strt),as.Date(end),1)          
# dates
# 
# ## [1] "2012-01-01" "2012-01-02" "2012-01-03" "2012-01-04"
# 
# # if extracting to shp then set up an output dataframe to collect
# dat_out <- as.data.frame(matrix(nrow = 0, ncol = 4))
# # else just plots
# par(mfrow = c(2,2))
# for(i in 1:length(dates)){
#   #  i=1
#   date_i <- dates[i]
#   infile <- sprintf("http://dapds00.nci.org.au/thredds/dodsC/rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/day/land/tmin/e_01/2012/eMAST_ANUClimate_day_tmin_v1m0_%s.nc", gsub("-", "", date_i))
#   
#   nc <- open.ncdf(infile)
#   vals <- get.var.ncdf(nc, varid="air_temperature")
#   nc.att <- nc$var$air_temperature
#   xmin <- min(nc.att$dim[[1]]$vals)
#   xmax <- max(nc.att$dim[[1]]$vals)
#   ymin <- min(nc.att$dim[[2]]$vals)
#   ymax <- max(nc.att$dim[[2]]$vals)
#   
#   print(c(xmin,xmax))
#   print(c(ymin,ymax))
#   
#   r <- raster(t(vals),
#               xmn=xmin, xmx=xmax,
#               ymn=ymin, ymx=ymax)
#   #str(r)
#   plot(r)
#   title(date_i)
#   #image(r)
#   e <- extract(r, shp, df=T)
#   #str(e) 
#   e1 <- shp@data
#   e1$date <- date_i
#   e1$values <- e[,2]
#   dat_out <- rbind(dat_out, as.data.frame(e1))
# }
# dat_out
# 
# # monthly
# "http://dap.nci.org.au/thredds/remoteCatalogService?command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/mon/land/prec/e_01/1970_2012/catalog.xml&dataset=rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/mon/land/prec/e_01/1970_2012/eMAST_ANUClimate_mon_prec_v1m0_197001.nc"
# strt <-'1976-01-01'
# end <- '2012-12-31'
# #dates <- paste(strt:end, 1:12)
# dates <- seq(as.Date(strt),as.Date(end),31) 
yy <- as.data.frame(1970:2012)
mm <- as.data.frame(c("01","02","03","04","05","06","07","08","09","10","11","12"))
library(sqldf)
dates <- sqldf("select * from yy, mm")
head(dates)
dates <- paste(dates[,1], dates[,2], sep = "")
head(dates)
dates[1:10]
#dat_out <- as.data.frame(matrix(nrow = 0, ncol = 4))
# else just plots
#par(mfrow = c(2,2))
# setwd("data")
cfiles <-  dir(pattern="tif$")
cfiles[1:10]
tail(cfiles)
length(cfiles)
dates[426]
system("df -h")
for(i in 426:length(dates)){
  #  i=183
  date_i <- gsub("-", "", substr(dates[i],1,7))
  print(date_i)
  infile <- sprintf("http://dapds00.nci.org.au/thredds/dodsC/rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/mon/land/prec/e_01/1970_2012/eMAST_ANUClimate_mon_prec_v1m0_%s.nc", gsub("-", "", date_i))
  #infile
  nc <- open.ncdf(infile)
  #str(nc)
  vals <- get.var.ncdf(nc, varid="lwe_thickness_of_precipitation_amount")
  nc.att <- nc$var$lwe_thickness_of_precipitation_amount
  xmin <- min(nc.att$dim[[1]]$vals)
  xmax <- max(nc.att$dim[[1]]$vals)
  ymin <- min(nc.att$dim[[2]]$vals)
  ymax <- max(nc.att$dim[[2]]$vals)
  
  #  print(c(xmin,xmax))
  #  print(c(ymin,ymax))
  
  r <- raster(t(vals),
              xmn=xmin, xmx=xmax,
              ymn=ymin, ymx=ymax)
  #str(r)
  #  plot(r)
  #  title(date_i)
  writeRaster(r, paste("precipitation_",date_i,".tif",sep=''), format="GTiff")
              
  #image(r)
  #e <- extract(r, shp, df=T)
  #str(e) 
  #e1 <- shp@data
  #e1$date <- date_i
  #e1$values <- e[,2]
  #dat_out <- rbind(dat_out, as.data.frame(e1))
  #write.table(e1, file.path("kwrt_weather_drought_1888_2014_p141_output_grids_emast.csv"),
  #            col.names = i == 1, append = i>1 , sep = ",", row.names = FALSE)
  Sys.sleep(time = 2)
}
#dat_out


cfiles <-  dir(pattern="tif$")
cfiles[1:10]
tail(cfiles)
for(i in seq_len(length(cfiles))){
  #i <- 1 ## for stepping thru
  gridname <- cfiles[[i]]
  r <- raster(gridname)
  e <- extract(r, shp, df=T)
  e1 <- shp
  e1@data$values <- e[,2]
  e1@data$gridname <- gridname
  # e1@data
  # write to to target file
  write.table(e1@data, "kwrt_weather_drought_1888_2014_p141_output_grids_emast.csv",
              col.names = i == 1, append = i>1 , sep = ",", row.names = FALSE)
}

############################################################################################

dir(pattern="csv")
dat <- read.csv("kwrt_weather_drought_1888_2014_p141_output_grids_emast.csv")
head(dat)
tail(dat)

dat$raster_layer <- as.character(dat$gridname)
dat$date <- matrix(unlist(strsplit(dat$raster_layer, "_")), ncol = 2, byrow=TRUE)[,2]
head(dat)
dat$date <- gsub(".tif","",dat$date)
head(dat )
dat$date <- paste(substr(dat$date,1,4), substr(dat$date,5,6), 1, sep = "-")
head(dat )
dat$year <- substr(dat$date,1,4)
dat$month <- substr(dat$date,6,7)
dat$year <- as.numeric(dat$year)
dat$month <- as.numeric(dat$month)
dat$date <- as.Date(dat$date)
str(dat)


qc <- dat


require(HutchinsonDroughtIndex)
qc$rain <- qc$values
as.data.frame(table(qc$year))
indat <- qc[,c("date","year","month","rain")]
str(indat)
indat[(nrow(indat) - 20):nrow(indat),]


drt <- drought_index_stations(data=indat,
                              years=length(names(table(indat$year)))
)
head(drt)
tail(drt)
str(drt)
qc3 <- drt[drt$year>=1979 & drt$year < 1984,]

png(file.path(homedir,"kwrt_weather_drought_1888_2014_p141_output3.png"))
par(mfrow=c(4,1),mar=c(2.5,2,1.5,1))
plot(qc3$date,qc3$rain,type='l',main='Prospect Reservoir (67019) NSW: raw monthly rainfall')
#points(qc3$date,qc3$rain)

lines(qc3$date,qc3$sixmnthtot/6, lwd = 2) #,type='l',main='6-monthly total rainfall')
points(qc3$date,qc3$sixmnthtot/6)

plot(qc3$date,qc3$index,type='l',main='Rescaled percentiles -4 to +4, -1 is Palmer Index Mild Drought',ylim=c(-4,4))
points(qc3$date,qc3$index)
segments(min(qc3$date),-1,max(qc3$date),-1)
segments(min(qc3$date),0,max(qc3$date),0,lty=2)
plot(qc3$date,qc3$count,type='l',main='Counts below -1 threshold, count of 5 or more is a drought')
points(qc3$date,qc3$count)
segments(min(qc3$date),5,max(qc3$date),5)

plot(qc3$date,qc3$count2,type='l',main='Enhanced counts of months if already passed count of 5 and percentiles less than 50%')
points(qc3$date,qc3$count2)
segments(min(qc3$date),5,max(qc3$date),5)
dev.off()

