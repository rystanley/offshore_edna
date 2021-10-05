#Code for formatting the RV data for comparison to the eDNA

#load libraries ----------
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(rnaturalearth)
library(rnaturalearthhires)
library(lubridate)
library(taxize)
library(pbapply)
library(gsubfn)
library(Mar.datawrangling)

#load functions -----------
source("code/ClassifyFunction.R")

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

# #load RV data
# get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVdata")
# rvdat <- summarize_catches()

#Basemaps
bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%
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

#maritimes network -- note that this cannot be shared or made available online to the public
maritimes_network <- read_sf("data/Shapefiles/networksites_proposed_OEM_MPA_v2.shp")%>%
  st_transform(latlong)%>%
  mutate(status=STATUS,
         status=ifelse(status=="Existing",status,
                       ifelse(status=="Proposed AOI","AOI","Proposed")),
         status=factor(status,levels=c("Existing","AOI","Proposed")))

#bioclassification polygons
bioclass <- read_sf("data/Shapefiles/MaritimesBioclassificationPolygons.shp")%>%
  st_transform(latlong)%>%
  data.frame()%>%
  st_as_sf(crs=latlong)

#Shelfbreak contour
# GEBCO <- raster("c:/Users/stanleyr/Documents/Github/CanadaMPAs/data/Bathymetry/GEBCO/gebco_2019_Canada.tif") #too big for github
# 
# bathy <- crop(GEBCO,as_Spatial(bioregion_box%>%st_transform(proj4string(GEBCO))))
# bathy[bathy>2] <- NA
# bathy[bathy<=-250] <- NA
# bathy[!is.na(bathy)] <- 1 #all one value for target depth range. UP on land (1m : -250m)
# 
# shelfbreak <- st_as_stars(bathy)%>%
#   st_as_sf(as_points = FALSE, merge = TRUE)%>%
#   st_transform(latlong)

#st_write(shelfbreak,"data/Countour_250.shp")

shelfbreak=read_sf("data/Shapefiles/Countour_250.shp")

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
# un = "stanleyr"
# pw = "" #type int the password
#
# get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata/2021",
#             fn.oracle.username = un,
#             fn.oracle.password = pw,
#             fn.oracle.dsn = "PTRAN",
#             force.extract = TRUE)
# rvdat <- summarize_catches()
# #download the species taxonomic IDs that can be used to match up to itis
# Mar.datawrangling::get_data_custom(schema= "GROUNDFISH", data.dir = "R:/Science/CESD/HES_MPAGroup/Data/RVdata", 
#                                    tables = c("GSSPECIES_CODES"), usepkg = "rodbc", fn.oracle.username = un", 
#                                    fn.oracle.password = pw, fn.oracle.dsn = "PTRAN")
# rvdat <- rvdat%>%
#   mutate(year=year(SDATE))%>%
#   filter(year==2020, #filter to just the year of the mission
#          MISSION.GSCAT == "NED2020025")

# save(x=rvdat,file="data/RVlinked.RData")
load("data/RVlinked.RData")
load("data/GROUNDFISH.GSSPECIES_CODES.RData")
spec_tax <- GSSPECIES_CODES
rm(GSSPECIES_CODES)

rvdat_all <- rvdat

rvdat <- rvdat_all%>%
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
sets_fix <- read.csv("data/eDNA_Trawl_Catch.csv") 


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

merged_master[merged_master$SurveyID != merged_master$HYDRO,] #just the typo that we fixed above


#make a clean file and do an intersect with other data for things like depth and where the sites are

#Data you want from the rv_survey
rvdat <- rvdat%>% mutate(year=year(SDATE),
                         month=month(SDATE),
                         day=day(SDATE))

rv_names <- c("SETNO","year","month","day","TIME","DEPTH","DIST","LONGITUDE","LATITUDE","STRAT","STATION", #note did a check with the stations that were written down for processing and they match the RV set
              "BOTTOM_TEMPERATURE","BOTTOM_SALINITY","SURFACE_TEMPERATURE")

merged_df <- merged_master%>%
             left_join(.,rvdat%>%dplyr::select(all_of(rv_names))%>%distinct(SETNO,.keep_all=TRUE))%>%
             dplyr::select(-c(Set, SurveyID,Station,Date_Collected))%>%
             rename(edna_sample_num = Sample.)%>%
             rename_all(tolower)%>% #make them all lowercase to make it easier to code
             filter(setno != 201)%>% #there is no coordinates in the db for this one nor is there any catch anyway. Need to check with Jamie to figure this one out
             st_as_sf(coords=c("longitude","latitude"),crs=latlong,remove=FALSE)%>%
             st_join(.,bioclass%>%dplyr::select(name)%>%rename(bioclass=name),join=st_intersects)%>%
             st_join(.,maritimes_network%>%dplyr::select(NAME)%>%rename(network=NAME),join=st_intersects)%>%
             mutate(network=ifelse(is.na(network),"Outside",network),
                    nearest_site = network, #these are placeholders for the for loop that evaluates site by site what the nearest network site is, if it was not direclty in a site
                    nearest_site_distance = 0)

#key in the sites that are in very close proximity to MPAs
sites_outside <- which(merged_df$network == "Outside")

for(i in sites_outside){
  
  temp <- merged_df[i,] #just sites that are 'outside' existing sites within the proposed network
  distance_network <- as.numeric(st_distance(temp,maritimes_network)/1000)#distance in km from all sites in the network
  
  #fill in the dataframe
  merged_df[i,"nearest_site"] = maritimes_network[which.min(distance_network),]%>%pull(NAME) #name of the closest site
  merged_df[i,"nearest_site_distance"] = distance_network[which.min(distance_network)] #distance to the closest site
  
  
}

output_merged <- merged_df%>%data.frame()%>%dplyr::select(-geometry)

#Save the csv that matches the stations             
write.csv(output_merged,"data/merged_edna_samples.csv",row.names=FALSE)

##Create map of sites
p1 <- ggplot()+
  geom_sf(data=bioreg,fill=NA)+
  geom_sf(data=bioclass,aes(fill=name))+
  geom_sf(data=maritimes_network,fill="grey50",alpha=0.5)+
  geom_sf(data=merged_df)+
  theme_bw()+
  labs(fill="")+
  theme(legend.position="bottom");p1

ggsave("inst/2020_edna_survey.png",p1,width=6,height=6,units="in",dpi=300)

### Do the taxonomic classifications -------

#key out what was ID in the RV survey for which we have a paired eDNA sample ** this was when running a custom script. I was able to get access to 
#the GSSPECIES_CODES table (PTRAN) which has linked TSN numbers for each species, that can be better matched to itis

# id_names <- rvdat%>%
#               filter(SETNO %in% merged_df$setno)%>%
#               distinct(SPEC.GSSPECIES,.keep_all=TRUE)%>%
#               dplyr::select(SPEC.GSSPECIES,COMM)%>%
#               mutate(latin = SPEC.GSSPECIES, #clean up the names for entry into various online taxonomic databases via APIs and the taxise package.#basic data cleaning
#                      latin = gsub("SPP.","G.",latin,fixed=T),#this specifies that these are to the genus level
#                      latin = gsub("SP.","G.",latin,fixed=T), 
#                      latin = gsub(" SP","G.",latin,fixed=T),
#                      latin = gsub("S.P.","G.",latin,fixed=T),
#                      latin = gsub("S.C.","",latin,fixed=T), #Not sure what is meant here (superclass?) but we can ID based on the db's to class. 
#                      
#                      latin = gsub("COLUS G.","BUCCINIDAE F.",latin), #COLUS SP are a genus of Buccinidae but this doesn't return anything in itis. This (like Hayus) will have to modified after to add the family name.
#                      latin = gsub("HYAS G.","OREGONIIDAE F.",latin),#these are specified as toad crabs (genus hyas)# this will fix the gsub error caused by the first hyas sub but for some reason the online db's return plants. Have to add 'hyas' to the tax id manually
#                      
#                      latin = gsub("BRYOZOANS BRACHIOPODA P.","BRACHIOPODA P.",latin,fixed=T), #for some reason COMM 'Lampshells' are listed as under a grouped phylum id. In this case we expect them to be part of the brachiopoda and not bryozoa
#                      latin = gsub("OREGONIIDAE F. ARANEUS","HYAS ARANEUS",latin),# this will fix the gsub error caused by the first hyas sub
#                      latin = gsub("OREGONIIDAE F. COARCTATUS","HYAS COARCTATUS",latin),# this will fix the gsub error caused by the first hyas sub
#                      latin = gsub("PENNATULACEA","PENNATULACEA O.",latin), #this will specify that this is to the order level. 
#                      latin = gsub("OPHIUROIDEA","OPHIUROIDEA C.",latin), #this will specify that this is to the order level. 
#                      latin = gsub("GORGONOCEPHALIDAE,ASTERONYCHIDAE F.","PHRYNOPHIURIDA O.",latin), #there are two families grouped for basket stars, so we will go up one taxonomic level to order Phrynophiurida 
#                      latin = ifelse(grepl("GEPHYREA",latin),"SIPUNCULA P.",latin), #Gephyrea is a former taxon group which now is split into three phyla including sipuncula. This is a worm
#                      
#                      latin = trimws(latin))%>% #get rid of trailing whitespace
#   
#               filter(!grepl("eggs",tolower(latin)),# we have skate, whelk and snail/slug eggs - these are removed. 
#                       latin != "ORGANIC DEBRIS", #remove
#                       latin != "UNIDENTIFIED", #remove
#                       latin != "ASCOPHYLLUM NODOSUM", #remove - seaweed would just be captured on the upcast since the depths are much deeper than any marine algae/plants would be expected from
#                       latin !="STONES AND ROCKS", #remove
#                       latin !="PURSE LITTLE SKATE", #this is an egg case (purse)
#                      #!grepl("worm",tolower(latin)), #need to decide whether to remove these -- aslo would remove SIPUNCULA P. 
#                      )
# 

#id_names <- id_names[c(which(grepl(".",id_names$latin,fixed=T)),which(!grepl(".",id_names$latin,fixed=T))),] #for review

  id_names <- rvdat%>%
                filter(SETNO %in% merged_df$setno)%>%
                distinct(SPEC.GSSPECIES,.keep_all=TRUE)%>%
                dplyr::select(SPEC.GSSPECIES,CODE)%>%
                rename(latin=SPEC.GSSPECIES)%>%
                filter(!grepl("eggs",tolower(latin)),# we have skate, whelk and snail/slug eggs - these are removed. 
                               latin != "ORGANIC DEBRIS", #remove
                               latin != "UNIDENTIFIED", #remove
                               latin != "ASCOPHYLLUM NODOSUM", #remove - seaweed would just be captured on the upcast since the depths are much deeper than any marine algae/plants would be expected from
                               latin !="STONES AND ROCKS", #remove
                               latin !="PURSE LITTLE SKATE" #this is an egg case (purse)
                              #!grepl("worm",tolower(latin)), #need to decide whether to remove these -- aslo would remove SIPUNCULA P.
                              )%>%
                left_join(.,spec_tax)

    #species that can't be matched to ITIS using TSN (either unavailable or not confirmed)
    problem_species_codes <- c(id_names[which(is.na(id_names$TSN)),"CODE"],id_names%>%filter(!is.na(TSN_DEFINITIVE),!TSN_DEFINITIVE)%>%pull(CODE))%>%unique()
    problem_species <- id_names%>%filter(CODE %in% problem_species_codes)
  
    
  #Now extract data from itis
  PhyloNames <- c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","superclass","class","superorder","order","family","subfamily","genus","species")      
  
  id_taxonomy <- classification(id_names%>%filter(!CODE %in% problem_species_codes)%>%pull(TSN), db = 'itis',return_id = FALSE)
  
  taxInfo <- do.call("rbind", id_taxonomy)%>%
              mutate(TSN = trunc(as.numeric(rownames(.))))%>%
              group_by(TSN)%>%
              spread(rank,name)%>%
              ungroup()%>%
              dplyr::select(TSN,all_of(PhyloNames))
  
  #merge this with the taxonomic table that we got the TSN from
  species_df <- id_names%>%
                left_join(.,taxInfo)
  
  #now deal with the problem species that were not keyed out. There aren't that many of them, so it is easier to hard code
  
  #American Fourspot Flounder - itis can't ID the species, so we will pull the data for the genus and then manually add the species
      temp <- classification("HIPPOGLOSSINA",db="itis")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)%>%mutate(TSN = NA,species="Hippoglossina oblonga")
      
      colmatch <- intersect(colnames(taxInfo),colnames(temp))  
      
      species_df[problem_species%>%filter(grepl("HIPPOGLOSSINA",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]        
      
  #Lampshells - this for some reason has Byrozoans and Brachiopods. Given they are called 'lamp shells' we will classify them as phylum Brachiopoda                                               
      temp <- classification("BRACHIOPODA",db="itis")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)%>%mutate(TSN = NA)
      
      colmatch <- intersect(colnames(taxInfo),colnames(temp))  
      
      species_df[problem_species%>%filter(grepl("BRYOZOANS BRACHIOPODA",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]    
  
  #Gephyrea annelid worms - now deprecated to Sipuncula
      temp <- classification("SIPUNCULA",db="itis")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)%>%mutate(TSN = NA)
        
      colmatch <- intersect(colnames(taxInfo),colnames(temp))  
    
      species_df[problem_species%>%filter(grepl("SIPUNCULA",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]    
    
  #Neptunea decemcostata - this is available on worms with the aphiaID of 491164 ** note the change in db
        temp <- classification("491164",db="worms")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)
        colnames(temp) <- tolower(colnames(temp))
        temp$TSN <- NA
        
        colmatch <- intersect(colnames(taxInfo),tolower(colnames(temp))) 
        
        species_df[problem_species%>%filter(grepl("DECEMCOSTATA",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]    
    
   #Ceramaster granularis - this is available on worms with the aphiaID of 124020 ** note the change in db
    temp <- classification("124020",db="worms")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)
        colnames(temp) <- tolower(colnames(temp))
        temp$TSN <- NA
        
        colmatch <- intersect(colnames(taxInfo),tolower(colnames(temp))) 
        
        species_df[problem_species%>%filter(grepl("CEREMASTER GRANULARI",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]    

    #Red cushion star - as far as I can tell here the SPEC name is just misspelled. Porania pulvillus works with itis vs. PORANIA PULVILIS in the RV database
        temp <- classification("porania pulvillus",db="itis")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)%>%mutate(TSN = NA)
        
        colmatch <- intersect(colnames(taxInfo),colnames(temp))  
        
        species_df[problem_species%>%filter(grepl("PORANIA PULVILIS",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]    
    
    #Russian Hats - Vazella pourtalesii - I think the RV just misspelled. This is on WORMS with the AphiaID of 172121 ** note changes to db
        temp <- classification("172121",db="worms")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)
        colnames(temp) <- tolower(colnames(temp))
        temp$TSN <- NA
        
        colmatch <- intersect(colnames(taxInfo),tolower(colnames(temp))) 
        
        species_df[problem_species%>%filter(grepl("VAZELLA POURTALESI",SPEC))%>%pull(SPEC)==species_df$SPEC,colmatch] <- temp[,colmatch]    

     #TUNICATA S.P. - itis only keys out the subphylum and not the kingdom/phylum
        temp <- classification("146420",db="worms")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)
        colnames(temp) <- tolower(colnames(temp))
        temp$TSN <- NA
        
        colmatch <- intersect(colnames(taxInfo),tolower(colnames(temp))) 
        
        species_df["TUNICATA S.P."==species_df$SPEC,colmatch] <- temp[,colmatch] 
        
      #PYCNOGONUM LITTORALE - marine amphipod
        temp <- classification("239867",db="worms")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)
        colnames(temp) <- tolower(colnames(temp))
        temp$TSN <- NA
        
        colmatch <- intersect(colnames(taxInfo),tolower(colnames(temp))) 
        
        species_df["PYCNOGONUM LITTORALE"==species_df$SPEC,colmatch] <- temp[,colmatch] 
        
      #Paguridae (hermit crabs) aren't keyed out with the TSN number nor the APHIAID with worms. but the name search in itis works to the family
      #Paguridae that corresponds to the Paguridae F. in the RV survey data. 
      
        temp <- classification("Paguridae",db="itis")%>%.[[1]]%>%data.frame()%>%dplyr::select(rank,name)%>%spread(rank,name)
        colnames(temp) <- tolower(colnames(temp))
        temp$TSN <- NA
        
        colmatch <- intersect(colnames(taxInfo),tolower(colnames(temp))) 
        
        species_df["PAGURIDAE F."==species_df$SPEC,colmatch] <- temp[,colmatch] 
        
  #unique ID where it will get the scientific name from the RV survey if it was not key'ed to species by the itis classification
        species_df$species_id <- species_df$species
        
        for(i in species_df[is.na(species_df$species),"CODE"]){
          
          temp <- filter(species_df,CODE == i)%>%dplyr::select(all_of(PhyloNames))%>%gather()
          
          id <-  temp[rev(which(!is.na(temp$value)))[1],"value"] #first taxonomic class back from 'species' that has information 
          
          species_df[species_df$CODE == i,"species_id"] <- id #integrate into the data
          
        }
        
    #make a column indicating 'vertebrates' and 'invertebrates' using the phylum 'Chordata' - this is a logical 'are you a vertebrate' which is useful for indexing
    species_df <- species_df%>%mutate(vert=ifelse(phylum == "Chordata",TRUE,FALSE))
    
#save outputs
 write.csv(species_df%>%dplyr::select(-c(ID_SRC,COMMENTS,TSN_SPELLING,TSN_SRC,TSN_SVC,APHIAID_SPELLING,
                                         APHIAID_DEFINITIVE,APHIAID_SVC,APHIAID_MULTI_FLAG,TSN_MULTI_FLAG)),"data/rv_identified.csv",row.names=FALSE)

 save(species_df,file="data/rv_species_identification.RData")
 

## Format data in sample by species format where sample is the SETNO that be linked back in later------
 
 data_reformat <- rvdat%>%
                   filter(SETNO %in% merged_df$setno)%>%
                   mutate(std_count = TOTNO*1.75/DIST, #standardized to a set tow distance of 1.75 nm
                          std_wgt = TOTWGT*1.75/DIST)%>%
                   dplyr::select(SETNO,CODE,std_count,std_wgt)%>%
                   left_join(.,species_df%>%dplyr::select(CODE,species_id,vert))%>%
                   filter(!is.na(species_id)) %>%#pre-excluded, non-taxonomic items (eg PURSE LITTLE SKATE,SKATE UNID. EGGS, WHELK EGGS (NS), ORGANIC DEBRIS, KNOTTED WRACK,FIRST UNIDENTIFIED PER SET,STONES AND ROCKS)
                   dplyr::select(-c(CODE))
                  
  ##by counts
 count_wide_all <- data_reformat%>%
                    dplyr::select(-c(std_wgt,vert))%>%
                    spread(key=species_id,value=std_count,fill=NA)%>%
                    left_join(output_merged%>%dplyr::select(setno,edna_sample_name)%>%rename(SETNO=setno))%>% #match up to the edna sample names
                    column_to_rownames('edna_sample_name')%>%
                    dplyr::select(-SETNO)
 
 count_wide_verts <- data_reformat%>%
                     filter(vert)%>%
                     dplyr::select(-c(std_wgt,vert))%>%
                     spread(key=species_id,value=std_count,fill=NA)%>%
                     left_join(output_merged%>%dplyr::select(setno,edna_sample_name)%>%rename(SETNO=setno))%>% #match up to the edna sample names
                     column_to_rownames('edna_sample_name')%>%
                     dplyr::select(-SETNO)
   
 count_wide_inverts <- data_reformat%>%
                       filter(!vert)%>%
                       dplyr::select(-c(std_wgt,vert))%>%
                       spread(key=species_id,value=std_count,fill=NA)%>%
                       left_join(output_merged%>%dplyr::select(setno,edna_sample_name)%>%rename(SETNO=setno))%>% #match up to the edna sample names
                       column_to_rownames('edna_sample_name')%>%
                       dplyr::select(-SETNO)
 
 #by weights 
 wgt_wide_all <- data_reformat%>%
                 dplyr::select(-c(std_count,vert))%>%
                 spread(key=species_id,value=std_wgt,fill=NA)%>%
                 left_join(output_merged%>%dplyr::select(setno,edna_sample_name)%>%rename(SETNO=setno))%>% #match up to the edna sample names
                 column_to_rownames('edna_sample_name')%>%
                 dplyr::select(-SETNO)
               
 wgt_wide_verts <- data_reformat%>%
                   filter(vert)%>%
                   dplyr::select(-c(std_count,vert))%>%
                   spread(key=species_id,value=std_wgt,fill=NA)%>%
                   left_join(output_merged%>%dplyr::select(setno,edna_sample_name)%>%rename(SETNO=setno))%>% #match up to the edna sample names
                   column_to_rownames('edna_sample_name')%>%
                   dplyr::select(-SETNO)
 
 wgt_wide_inverts <- data_reformat%>%
                     filter(!vert)%>%
                     dplyr::select(-c(std_count,vert))%>%
                     spread(key=species_id,value=std_wgt,fill=NA)%>%
                     left_join(output_merged%>%dplyr::select(setno,edna_sample_name)%>%rename(SETNO=setno))%>% #match up to the edna sample names
                     column_to_rownames('edna_sample_name')%>%
                     dplyr::select(-SETNO)
  
 #save outputs
     save(count_wide_all,file="data/all_species_count_wide.RData") #standardized counts
     save(count_wide_verts,file="data/verts_count_wide.RData")
     save(count_wide_inverts,file="data/inverts_count_wide.RData")
     
     save(wgt_wide_all,file="data/all_species_wgt_wide.RData") #standardized weights
     save(wgt_wide_verts,file="data/verts_wgt_wide.RData")
     save(wgt_wide_inverts,file="data/inverts_wgt_wide.RData")
     
     write.csv(count_wide_all,file="data/all_species_count_wide.csv") #standardized counts
     write.csv(count_wide_verts,file="data/verts_count_wide.csv")
     write.csv(count_wide_inverts,file="data/inverts_count_wide.csv")
     
     write.csv(wgt_wide_all,file="data/all_species_wgt_wide.csv") #standardized weights
     write.csv(wgt_wide_verts,file="data/verts_wgt_wide.csv")
     write.csv(wgt_wide_inverts,file="data/inverts_wgt_wide.csv")
    

##now investigate which species are the most numerate -------------
 
 species_select <- rvdat%>%
                   filter(CODE %in% species_df$CODE,SETNO %in% merged_df$setno)%>%
                   mutate(std_count = TOTNO*1.75/DIST, #standardized to a set tow distance of 1.75 nm
                          std_wgt = TOTWGT*1.75/DIST)%>%
                   group_by(CODE)%>%
                   summarise(num_sets=n(), #number of stations
                             count=sum(std_count),
                             catch=sum(std_wgt))%>%
                  ungroup()%>%
                  data.frame()%>%
                  left_join(.,species_df%>%dplyr::select(-c(ID_SRC,COMMENTS,TSN_SPELLING,TSN_SRC,TSN_SVC,APHIAID_SPELLING,
                                                           APHIAID_DEFINITIVE,APHIAID_SVC,APHIAID_MULTI_FLAG,TSN_MULTI_FLAG)))%>%
                  mutate(vert_invert=ifelse(phylum == "Chordata","Vertebrate","Invertebrate"))

 #unique ID where it will get the scientific name from the RV survey if it was not key'ed to species by the itis classification
 
 species_select$species_id <- species_select$species
 
 for(i in species_select[is.na(species_select$species),"CODE"]){
   
   temp <- filter(species_select,CODE == i)%>%dplyr::select(all_of(PhyloNames))%>%gather()
   
   id <-  temp[rev(which(!is.na(temp$value)))[1],"value"] #first taxonomic class back from 'species' that has information 
   
   species_select[species_select$CODE == i,"species_id"] <- id #integrate into the data
   
  }
 
 #change the plotting labels order for ggplot
 species_select <- species_select%>%
                    mutate(species_cnt_ord = factor(species_id,levels=species_select%>%arrange(count)%>%pull(species_id)),
                           species_catch_ord = factor(species_id,levels=species_select%>%arrange(catch)%>%pull(species_id)),
                           num_sets_ord = factor(species_id,levels=species_select%>%arrange(num_sets)%>%pull(species_id)))
  
 p2 <- ggplot()+
   geom_bar(stat="identity",data=species_select,aes(y=count,x=species_cnt_ord))+
   facet_wrap(~vert_invert,scales="free_y",ncol=1)+
   coord_flip()+
   theme_bw()+
   scale_y_log10()+
   labs(y=expression(paste("Log"[10]," sum standardized count",sep=" ")),
        x="")+
   theme(strip.background = element_rect(fill="white"),
         axis.text.y = element_text(size=6));p2
 
ggsave("output/species_counts_rv.png",p2,width=5,height=8,units="in",dpi=300) 

#just the top 10 for verts and inverts

top10_count <- species_select%>%
         arrange(-count)%>%
         group_by(vert_invert)%>%
         slice(1:10)%>%
         ungroup()%>%
         data.frame()

#by overall count summed across sets
    p3 <- ggplot()+
      geom_bar(stat="identity",data=species_select%>%filter(CODE %in% top10_count$CODE),aes(y=count,x=species_cnt_ord))+
      facet_wrap(~vert_invert,scales="free_y",ncol=1)+
      coord_flip()+
      theme_bw()+
      scale_y_log10()+
      labs(y=expression(paste("Log"[10]," sum standardized count",sep=" ")),
           x="")+
      theme(strip.background = element_rect(fill="white"),
            axis.text.y = element_text(size=10));p3
    
    ggsave("output/top10_count.png",p3,width=5,height=7,units="in",dpi=300)

#by number of sets
    top10_sets <- species_select%>%
      arrange(-num_sets)%>%
      group_by(vert_invert)%>%
      slice(1:10)%>%
      ungroup()%>%
      data.frame()
     
    p4 <- ggplot()+
      geom_bar(stat="identity",data=species_select%>%filter(CODE %in% top10_sets$CODE),aes(y=num_sets,x=num_sets_ord))+
      facet_wrap(~vert_invert,scales="free_y",ncol=1)+
      coord_flip()+
      theme_bw()+
      labs(y="Total # sets",x="")+
      theme(strip.background = element_rect(fill="white"),
            axis.text.y = element_text(size=10));p4
    
    ggsave("output/top10_num_sets.png",p4,width=5,height=7,units="in",dpi=300)

#By abundance 
    top10_abund <- species_select%>%
      arrange(-catch)%>%
      group_by(vert_invert)%>%
      slice(1:10)%>%
      ungroup()%>%
      data.frame()
    
    p5 <- ggplot()+
      geom_bar(stat="identity",data=species_select%>%filter(CODE %in% top10_abund$CODE),aes(y=catch,x=species_catch_ord ))+
      facet_wrap(~vert_invert,scales="free_y",ncol=1)+
      coord_flip()+
      theme_bw()+
      labs(y=expression(paste("Log"[10]," sum standardized abundance",sep=" ")),
           x="")+
      theme(strip.background = element_rect(fill="white"),
            axis.text.y = element_text(size=10));p5
    
    ggsave("output/top10_num_sets.png",p4,width=5,height=7,units="in",dpi=300)


