setwd("/home/danielred/data/programs/aquaticTS/src")
library(geodata)
library(ggplot2)
library(dplyr)
library(tidyr)

# settings
point_x = 14.480255
point_y = 48.975658
width_x = 0.1
width_y = 0.1
nametoS = "/home/danielred/data/programs/aquaticTS/IN/CB"


#timewindow <- "2061-2080"
#timewindow <- "2021-2040"
timewindow <- "2041-2060"

nameto <- paste0(nametoS, timewindow, ".tsv")

# get climate data
climate <- xfun::cache_rds({
  biomin <- cmip6_tile(15, 49, "CNRM-CM6-1-HR", "585", timewindow, var="tmin", res=10, path=tempdir())
  #bioavg <- cmip6_tile(15, 49, "CNRM-CM6-1-HR", "585", "2061-2080", var="tavg", res=10, path=tempdir())
  biomax <- cmip6_tile(15, 49, "CNRM-CM6-1-HR", "585", timewindow, var="tmax", res=10, path=tempdir())
  
  list( biomin=biomin, biomax=biomax )
}, hash=list(timewindow))

# set area
myext = ext(point_x-width_x/2, point_x+width_x/2, point_y-width_y/2, point_y+width_y/2)

plot( crop(climate[[1]], myext) ) # plot tmin by month

# crop it to area
tile_tmin = crop(climate[[1]], myext)
tile_tmax = crop(climate[[2]], myext)

# get mins and maxes (by month)
mins = sapply(tile_tmin, function(x){
  sum=0
  for(col in 1:ncol(x)) for(row in 1:nrow(x)) sum = sum + x[row, col]
  return( sum/(ncol(x)*nrow(x)))
})

maxs = sapply(tile_tmax, function(x){
  sum=0
  for(col in 1:ncol(x)) for(row in 1:nrow(x)) sum = sum + x[row, col]
  return( sum/(ncol(x)*nrow(x)))
})

b2 = as.data.frame(do.call(cbind, c(mins, maxs) )) |> 
  pivot_longer(cols= everything(), names_to=c("parameter", "month"), values_to = "val", names_pattern = "wc2.1_30s_(.*)_(.*)")

b3 <- b2 |> pivot_wider(names_from = parameter, values_from = val)
#b3 <- b3[,c("month", "tmin", "tmax")]

# plotting it
times <- seq(0, 1, by = 0.01)

b3_s <- b3 |> 
  uncount(max(times), .id="period") |>
  mutate(month=as.numeric(month)/12) |>
  arrange(month) |>
  summarise(spline_x = spline(month, tmin)[["x"]], 
            spline_hi = spline(month, tmax)[["y"]], 
            spline_lo = spline(month, tmin)[["y"]])
ggplot()+
  geom_ribbon(data=b3_s[b3_s$spline_x <= 1 & b3_s$spline_x >= 0,], aes(x=spline_x*12, ymax=spline_hi, ymin=spline_lo), fill="coral")+
  geom_linerange(data= b3, aes(x=as.numeric(month), ymin=tmin, ymax=tmax), linewidth=1.5)+
  labs(x="Month", y="Ground temperature [Celsius degree]", title=paste("Ceske Budejovice", timewindow))


# saving
write.table(b3, nameto, sep="\t", col.names = F, row.names = F)


# plotting on a map
library(sf)
library("rnaturalearth")
library("rnaturalearthdata")
library(ggspatial)
library(ggplot2)

point = data.frame(
  x=c(rep(point_x-width_x/2,2), rep(point_x+width_x/2, 2), point_x-width_x/2),
  y=c(point_y-width_y/2, rep(point_y+width_y/2, 2), rep(point_y-width_y/2, 2))
)

world <- ne_countries(scale = "medium", type = 'map_units', returnclass = "sf")

eu <- list(lat_from=34, lat_to=72, long_from=-25, long_to=45)

ggplot(data = world) +
  geom_sf()+
  coord_sf(xlim=c(eu$long_from, eu$long_to ), ylim=c(eu$lat_from, eu$lat_to)) +
  theme_bw()+
  geom_polygon(data=point, aes(x=x, y=y), fill="red", color="red")+
  xlab("Longitude (°E)") + ylab("Latitude (°N)")

