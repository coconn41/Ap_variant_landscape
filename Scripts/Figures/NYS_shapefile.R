library(sf)
library(dplyr)
library(readxl)
tdir=tempdir()
stateurl = "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_500k.zip"
uaurl = "https://www.dot.ny.gov/divisions/policy-and-strategy/darb/dai-unit/ttss/repository/2010UAUC.zip"
if(file.exists(paste(tdir,"/cb_2018_us_state_500k.shp",sep=""))==F){
  download.file(stateurl, destfile = file.path(tdir, "States.zip"))
  unzip(file.path(tdir,"States.zip"),exdir=tdir)}
NYS = read_sf(paste(tdir,"/cb_2018_us_state_500k.shp",sep="")) %>%
  filter(NAME=="New York") %>%
  st_transform(.,crs=32618)

if(file.exists(paste(tdir,"/urban2010.shp",sep=""))==F){
  download.file(uaurl, destfile = file.path(tdir, "Urban_areas.zip"))
  unzip(file.path(tdir,"Urban_areas.zip"),exdir=tdir)}
Urban_areas= read_sf(paste(tdir,"/urban2010.shp",sep="")) %>%
  st_transform(.,crs=32618)
