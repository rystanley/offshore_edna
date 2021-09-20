#load libraries
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthhires)

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

# #load RV data
# get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVdata")
# rvdat <- summarize_catches()

#Basemaps
bioregion <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/MaritimesPlanningArea.shp")%>%
  st_transform(latlong)

#bounding box for the bioregion
bioregion_box <- bioregion%>%
  st_transform(utmkm)%>%
  st_buffer(50)%>% #50km buffer
  st_transform(latlong)%>%
  st_bbox()%>%st_as_sfc()
#plotlimits
plotlims <- bioregion%>%
  st_transform(utmkm)%>%
  st_buffer(20)%>% #50km buffer
  st_transform(latlong)%>%
  st_bbox()

#maritimes network
maritimes_network <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Maritimes_Draft_Network/networksites_proposed_OEM_MPA_v2.shp")%>%
  st_transform(latlong)%>%
  mutate(status=STATUS,
         status=ifelse(status=="Existing",status,
                       ifelse(status=="Proposed AOI","AOI","Proposed")),
         status=factor(status,levels=c("Existing","AOI","Proposed")))

#bioclassification polygons
bioclass <- read_sf("R:/Science/CESD/HES_MPAGroup/Data/Shapefiles/MaritimesBioclassificationPolygons.shp")%>%
  st_transform(latlong)%>%
  data.frame()%>%
  st_as_sf(crs=latlong)

#Shelfbreak contour
# GEBCO <- raster("c:/Users/stanleyr/Documents/Github/CanadaMPAs/data/Bathymetry/GEBCO/gebco_2019_Canada.tif")
# 
# bathy <- crop(GEBCO,as_Spatial(bioregion_box%>%st_transform(proj4string(GEBCO))))
# bathy[bathy>2] <- NA
# bathy[bathy<=-250] <- NA
# bathy[!is.na(bathy)] <- 1 #all one value for target depth range. UP on land (1m : -250m)
# 
# shelfbreak <- st_as_stars(bathy)%>%
#   st_as_sf(as_points = FALSE, merge = TRUE)%>%
#   st_transform(latlong)

#st_write(shelfbreak,"R:/Science/CESD/HES_MPAGroup/Data/Bathymetry/Countour_250.shp")

shelfbreak=read_sf("R:/Science/CESD/HES_MPAGroup/Data/Bathymetry/Countour_250.shp")

#Create basemap intersected with the bounding box. 
basemap_atlantic <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="Canada"),
                          ne_states(country = "United States of America",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="USA"))%>%
  st_intersection(.,bioregion_box)


#load the data (originally pulled Dec 4 2020 - stanleyr homerun7 PTRAN) --------
#get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata")
#rvdat <- summarize_catches()
#save(x=rvdat,file="R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVlinked.RData")
load("R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVlinked.RData")

#load the eDNA survey station IDs -- can't quite figure out the coordinate reference system
edna_set_info <- read.csv("R:/Science/CESD/HES_MPAGroup/Projects/eDNA GRDI 2019-2022/data/Trawl2020_eDNASamples.csv")
hydros <- edna_set_info%>%filter(!grepl("LAB",eDNA_Sample_Name),
                                 SurveyID!="",
                                 !grepl("Blank",SurveyID))%>%
  pull(SurveyID)%>%as.numeric()


sets <- read.csv("R:/Science/CESD/HES_MPAGroup/Data/eDNA/eDNA_Trawl_Catch.csv")%>%
  filter(HYDRO%in%hydros)%>%
  distinct(SETNO)%>%
  pull(SETNO)

#filter the RV data          
edna_data <- rvdat%>%
  mutate(year=year(SDATE),
         month=month(SDATE))%>%
  filter(year==2020,
         month %in% c(7,8), #just the summer survey
         SETNO %in% sets)%>%
  st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong)

#Trawls without eDNA
null_sets <- rvdat%>%
  mutate(year=year(SDATE),
         month=month(SDATE))%>%
  filter(year==2020,
         month %in% c(7,8), #just the summer survey
         !SETNO %in% sets)%>%
  st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong)
