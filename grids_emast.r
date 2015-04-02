
#install.packages("ncdf", type = "source", configure.args="--with-netcdf-include=/usr/include")
require(ncdf)
#install.packages("raster")
require(raster)

infile <- "http://dapds00.nci.org.au/thredds/dodsC/rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/day/land/tmin/e_01/2012/eMAST_ANUClimate_day_tmin_v1m0_20120101.nc"
nc <- open.ncdf(infile)
str(nc)
names(nc$var)
tmaxscr.att <- nc$var$air_temperature
# get the maximum temperature screen values per time-step
tmaxscr.vals <- get.var.ncdf(nc, varid="air_temperature")
image(
  x=tmaxscr.att$dim[[1]]$vals,
  y=-1*tmaxscr.att$dim[[2]]$vals,
  z= tmaxscr.vals,
  xlab="Longitude",
  ylab="Latitude",
  main='test'
)
# Error in image.default(x = tmaxscr.att$dim[[1]]$vals, y = tmaxscr.att$dim[[2]]$vals,  : 
#                          increasing 'x' and 'y' values expected
#image(tmaxscr.vals)
head(tmaxscr.vals)
tmaxscr.crs <- get.var.ncdf(nc, varid="crs")
r <- raster(tmaxscr.vals,
            xmn=range(tmaxscr.att$dim[[1]]$vals)[1], xmx=range(tmaxscr.att$dim[[1]]$vals)[2],
            ymn=range(tmaxscr.att$dim[[2]]$vals)[1], ymx=range(tmaxscr.att$dim[[2]]$vals)[2])
#            crs=CRS("+proj=utm +zone=11 +datum=NAD83"))
plot(r)
#image(r)

"
Refs
http://ivanhanigan.github.io/2014/06/testing-thredds-server-with-r-code-from-tasmanian-tpac/
http://lukemiller.org/index.php/2011/03/extracting-sea-surface-temperatures-from-noaas-oisstv2/
http://ivanhanigan.github.io/2013/10/extract-weather-data-from-awap-grids/
http://www.emast.org.au/observations/climate/
http://portal.tern.org.au/daily-minimum-temperature-anuclimate-10-001-degree-australian-coverage-1970-2012
http://dap.nci.org.au/thredds/remoteCatalogService?command=subset&catalog=http://dapds00.nci.org.au/thredds/catalog/rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/day/land/tmin/e_01/2012/catalog.xml&dataset=rr9/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/day/land/tmin/e_01/2012/eMAST_ANUClimate_day_tmin_v1m0_20120101.nc
http://stackoverflow.com/questions/14513480/convert-matrix-to-raster-in-r
"
