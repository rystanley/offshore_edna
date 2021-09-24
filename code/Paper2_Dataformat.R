#load libraries
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(rnaturalearthhires)
library(lubridate)

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
# library(Mar.datawrangling)
# library(RODBC)
# get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata/2021",
#             fn.oracle.username = "stanleyr",
#             fn.oracle.password = "homerun7",
#             fn.oracle.dsn = "PTRAN",
#             force.extract = TRUE)
# rvdat <- summarize_catches()
# save(x=rvdat,file="R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVlinked.RData")
load("R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVlinked.RData")
rvdat_all <- rvdat

rvdat <- rvdat%>%
          mutate(year=year(SDATE))%>%
          filter(year==2020, #filter to just the year of the mission
                 MISSION.GSCAT == "NED2020025")

#load the eDNA survey station IDs
edna_set_info <- read.csv("data/Trawl2020_eDNASamples.csv")%>% #this is the sampling filtering register
                  mutate(SurveyID = ifelse(is.na(SurveyID),paste0(gsub("EDNA20-","",eDNA_Sample_Name),"_Blank"),SurveyID))

hydros <- edna_set_info%>%filter(!grepl("LAB",eDNA_Sample_Name),
                                 SurveyID!="",
                                 !grepl("Blank",SurveyID))%>%
                  pull(SurveyID)%>%as.numeric()

#this is the data from the RV survey that was provided post-survey. The 'HYDRO' here is the bottom water sample that corresponds to that set
sets_fix <- read.csv("R:/Science/CESD/HES_MPAGroup/Data/eDNA/eDNA_Trawl_Catch.csv") 


#these are 'orphaned' eDNA samples that don't correspond to a trawl set based on the Hydro (SampleID) recorded on the bottle for the data provided post-survey
setdiff(hydros,unique(sets_fix$HYDRO))

#Sample 478752 also had 'set' 43 recorded on the sample, which corresponds to Hydro 478452. This means that there was a typo
#on the bottle
sets_fix%>%filter(SETNO == 43)%>%pull(HYDRO)%>%unique()

#Will change this in the trawl data so it corresponds with the sampled eDNA ## This overwrites so makes sure the code is run in sequence
edna_set_info <- edna_set_info%>%mutate(HYDRO = ifelse(SurveyID == 478752,478452,SurveyID))

#see if we can match
hydros2 <- edna_set_info%>%filter(!grepl("LAB",eDNA_Sample_Name),
                                 SurveyID!="",
                                 !grepl("Blank",SurveyID))%>%
  pull(HYDRO)%>%as.numeric()

setdiff(hydros2,unique(sets_fix$HYDRO)) #confirm two remaining as orphans

#fix orphan hydro 474254
    #based on a review of the database Mike McMahon was able to identify that 474524 was from set 201 - refer to Orphaned_Hydro_Information.docx in HES_MPAGroup/Data/eDNA folder
    edna_set_info <- edna_set_info%>%mutate(Set = ifelse(HYDRO == 474254,201,Set))
    
    sets_fix%>%filter(SETNO == 201) # for some reason 201 is not in the dataset provided post-survey
    
    #Note that 201 is not associated with any trawl so this is a dude anyway 

#fix orphan hydro 478989
    #based on a similar review there is a strong potential for the remaining orphan to be part of set 116 based on the sequence of numbers. Now this is still an assumption so this data SETNO 116
    #should be viewed with caution in the comparative analysis until this is confirmed with Jaimie Emberly who is back in October
    
    edna_set_info <- edna_set_info%>%mutate(Set = ifelse(HYDRO == 478989,116,Set))
    
    sets_fix%>%filter(SETNO == 116) #The value set here is 99999 which would seem like an error in the database corresponding to why the hydro isn't in this data. 
    
    sets_fix <- sets_fix%>%mutate(HYDRO = ifelse(SETNO == 116,478989,HYDRO))

#Now merge sets_fix into edna_set_info so that the total set list will match up
#first check that for the sets that were recorded from the samples match up with the sets and corresponding hydros from the post-survey data

  #track load the RV trawl data and match to the sets this will have two orphaned sets
  sets <- sets_fix%>%
    filter(HYDRO%in%hydros2)%>%
    distinct(SETNO)%>%
    pull(SETNO)

    #Do a quick check to see if the sets recorded from the samples match with the SampleIDs
    set_check1 <-  edna_set_info%>%
                   filter(!is.na(Set),
                          !grepl("Blank",SurveyID))%>%
                   arrange(Set)%>%
                   dplyr::select(Set,SurveyID)
    
    set_check2 <- sets_fix%>%
                  distinct(SETNO,.keep_all=TRUE)%>%
                  filter(SETNO %in% set_check1$Set)%>%
                  arrange(SETNO)%>%
                  rename(Set=SETNO)%>%
                  dplyr::select(Set,HYDRO)

    #this is the mismatch that we noted earlier
    set_check1[!set_check1$SurveyID == set_check2$HYDRO,]
    set_check2[!set_check1$SurveyID == set_check2$HYDRO,]

#Merge the datasets 
merged_master <- edna_set_info%>%
                 filter(!grepl("Blank",HYDRO))%>% #filter out the blank ids
                 mutate(HYDRO = as.numeric(HYDRO))%>%
                  left_join(sets_fix%>%
                             filter(HYDRO %in% edna_set_info$HYDRO)%>% #get rid of the sets that weren't part of the processed samples from the lab sheet
                             dplyr::select(HYDRO,SETNO)%>%
                             distinct(HYDRO,.keep_all=TRUE))%>%
                  mutate(SETNO = ifelse(HYDRO == 474254,201,SETNO)) #this is the sub we did earlier
    
#quick check to see if the sets match. Should == 0
sum(!merged_master[!is.na(merged_master$Set),"Set"] == merged_master[!is.na(merged_master$Set),"SETNO"]) #should be 0

sum(is.na(merged_master$SETNO)) #should be 0

#Now match these samples to the RV data         
trawl_edna_data <- rvdat%>%
  mutate(year=year(SDATE),
         month=month(SDATE))%>%
  filter(SETNO %in% merged_master$SETNO)%>%
  st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong)

setdiff(unique(trawl_edna_data$SETNO)%>%sort(),merged_master$SETNO%>%sort())
setdiff(merged_master$SETNO%>%sort(),unique(trawl_edna_data$SETNO)%>%sort()) #this is because 201 doesn't appear to have a trawl but an eDNA sample was still taken




#Trawls without eDNA
null_sets <- rvdat%>%
  mutate(year=year(SDATE),
         month=month(SDATE))%>%
  filter(year==2020,
         month %in% c(7,8), #just the summer survey
         !SETNO %in% sets)%>%
  st_as_sf(coords=c("LONGITUDE","LATITUDE"),crs=latlong)
