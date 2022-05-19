### Code to make Figure 1 --- map

#load libraries ----------
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggmap)
library(ggspatial)

sf_use_s2 = FALSE

#load function to shrink polygons 

shrinkIfPossible <- function(sf, size) { # https://gis.stackexchange.com/questions/392505/can-i-use-r-to-do-a-buffer-inside-polygons-shrink-polygons-negative-buffer
  # compute inward buffer
  sg <- st_buffer(st_geometry(sf), -size)
  
  # update geometry only if polygon is not degenerate
  st_geometry(sf)[!st_is_empty(sg)] = sg[!st_is_empty(sg)]
  
  # return updated dataset
  return(sf)
}

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#Basemaps
bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%
  st_transform(latlong)

#bioclassification polygons
bioclass <- read_sf("data/Shapefiles/MaritimesBioclassificationPolygons.shp")%>%
  st_transform(latlong)%>%
  data.frame()%>%
  st_as_sf(crs=latlong)

#read in the sample coordinates
edna_samples <- read.csv("data/merged_edna_samples.csv")%>%
                st_as_sf(coords=c("longitude","latitude"),crs=latlong,remove=FALSE)

#bounding box for the bioregion
bioregion_box <- bioclass%>%
                st_combine()%>%
                st_bbox()%>%
                st_as_sfc()%>%
                st_as_sf()%>%
                st_buffer(0.5)

#plotlimits
plotlims <- bioclass%>%
            st_combine()%>%
            st_bbox()%>%
            st_as_sfc()%>%
            st_as_sf()%>%
            shrinkIfPossible(0.1)%>%
            st_bbox()

#map inset
inset_box <- bioregion%>%
            st_bbox()

inset_box[1] <- -80
inset_box[2] <- 30
inset_box[3] <- -52
inset_box[4] <- 56

inset_box <- inset_box%>%
  st_as_sfc()%>%
  st_as_sf()

#plotlimits
plotlims_inset <- shrinkIfPossible(inset_box,-1)%>%st_bbox()

#maritimes network -- note that this cannot be shared or made available online to the public
maritimes_network <- read_sf("data/Shapefiles/networksites_proposed_OEM_MPA_v2.shp")%>%
  st_transform(latlong)%>%
  mutate(status=STATUS,
         status=ifelse(status=="Existing",status,
                       ifelse(status=="Proposed AOI","AOI","Proposed")),
         status=factor(status,levels=c("Existing","AOI","Proposed")),
         status=ifelse(NAME =="Eastern Canyons","Existing",status))


#Shelfbreak contour
shelfbreak=read_sf("data/Shapefiles/Countour_250.shp")

#Create basemap intersected with the bounding box. 
basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
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
                            mutate(country="USA"))

#create the inset map first (this is the geographic context plot)
basemap_bioregion <- basemap%>%st_intersection(.,bioregion_box)

basemap_inset <- basemap%>%st_intersection(.,inset_box)

inset_map <- ggplot()+
              geom_sf(data=basemap)+
              geom_sf(data=bioregion_box,fill=NA)+
              coord_sf(xlim=plotlims_inset[c(1,3)],ylim=plotlims_inset[c(2,4)],label_graticule = "--NE")+
              theme_bw()+
              scale_x_continuous(breaks=c(-75,-70,-65,-60,-55))+
  theme(axis.text.y = element_text(margin = margin(-1,0,0,.5, unit = 'cm'),hjust=3),
            axis.text.x = element_text(vjust = 5, margin = margin(-0.5,0,0.5,0, unit = 'cm'),),
            axis.ticks.length = unit(-0.25,"cm"));inset_map

ggsave("output/inset_map.png",inset_map,height=6,width=6,units="in",dpi=600)


#map the more detailed regional view

sample_map <- ggplot()+
              geom_sf(data=basemap)+
              geom_sf(data=bioclass,aes(fill=name))+
              geom_sf(data=shelfbreak,fill=NA,col="grey40")+
              geom_sf(data=maritimes_network%>%filter(STATUS!="Proposed"),alpha=0.25)+
              geom_sf(data=edna_samples,size=2)+
              coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
              theme_bw()+
              labs(fill="")+
              annotation_scale(location="br",line_width=0.25,text_cex=1.5)+ #location is 'bottom right'
              #annotation_north_arrow(location="tr",height = unit(1,"in"),width=unit(1,"in")) +
              theme(axis.text = element_blank(),
                    axis.ticks = element_blank())

ggsave("output/sample_map.png",sample_map+theme(legend.position = "none"),height=12,width=12,units="in",dpi=600)
ggsave("output/sample_map_legend.png",sample_map+theme(legend.box.background = element_rect(colour = "black"),
                                                       legend.title = element_blank()),height=6,width=6,units="in",dpi=600)

