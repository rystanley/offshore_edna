## Code to develop full taxonomic hierarchies for species detected in the eDNA

#load libraries
library(dplyr)
library(tidyr)
library(taxize)
library(pbapply)
library(readxl)

#load species lists per marker

markers <- c("CO1","12s","16s","18s")

edna_species <- NULL

for(i in markers){
  
  temp <- read_excel("data/Xiaoping species lists paper 2.xlsx",sheet=i,col_names=FALSE)%>%
          rename(call = 1)%>%
          mutate(marker=i)
  
  edna_species <- rbind(edna_species,temp)
  
}

#create list of unique species that will be passed to itis for taxonmic id
edna_calls <- edna_species%>%pull(call)%>%unique()

#Refer to itis for the first run through. Note that some species return more than one option so it cannot be fully autonomous
edna_taxonomy <- classification(edna_calls, db = 'itis',return_id = TRUE) #this can take a while note 'decisions' below

      #During the classification process there are several calls that require a decision because ITIS returns too many options. These
      #were made 
      # 'Cyanea' - itis 34671 is the 'valid one' - 48
      # 'Halichondria panicea' - itis 48396 is the 'valid one' - 1
      # 'Suberites ficus' - itis 48488 is the 'valid one' - 2


#flatten so we can start to build a datatable
edna_flat <- do.call("rbind", edna_taxonomy)%>%
  mutate(call = gsub("\\..*","",rownames(.))) #this creates a column with what itis recieved as input

    #some species were not key'ed out by itis.
    non_classified <- edna_flat%>%
                      filter(is.na(name))%>%
                      pull(call)%>%
                      unique()
    
    length(non_classified) #should match (28) the 'Not Found' message from the results
    
    #create a data.frame that has the itis TSN for the lowest taxonomic rank for each 'call'. Note that 'rank' is retained so you can see that it takes the lowest rank per call.
      edna_itis <- NULL
      
      for(i in edna_flat%>%filter(!is.na(name))%>%pull(call)%>%unique()){
        
        temp <- edna_flat%>%
          filter(call == i)%>%
          filter(row_number()==n())%>%
          dplyr::select(rank,call,id)%>%
          rename(TSN=id)
        
        rownames(temp) <- NULL
        
        edna_itis <- rbind(edna_itis,temp)
        
      } # change stuff

  #now reformat the eDNA data for the species that have info
 
    PhyloNames <- c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","superclass","class","superorder","order","family","subfamily","genus","species")      

    edna_format <- edna_flat%>%
                    filter(!is.na(name))%>%
                    dplyr::select(name,rank,call)%>%
                    group_by(call)%>%
                    spread(rank,name)%>%
                    ungroup()%>%
                    dplyr::select(call,all_of(PhyloNames))%>%
                    left_join(.,edna_itis%>%dplyr::select(-rank))%>%
                    mutate(db="itis") #just to flag where the data is from


    ### now fix the problem species first with gbif
    edna_taxonomy_fix <- classification(non_classified, db = 'gbif',return_id = TRUE) #this can take a while note 'decisions' below
    
    edna_flat_fix <- do.call("rbind", edna_taxonomy_fix)%>%
      mutate(call = gsub("\\..*","",rownames(.)))
    
    edna_gbif <- NULL
    for(i in edna_flat_fix%>%filter(!is.na(name))%>%pull(call)%>%unique()){
      
      temp <- edna_flat_fix%>%
        filter(call == i)%>%
        filter(row_number()==n())%>%
        dplyr::select(rank,call,id)%>%
        rename(TSN=id)
      
      rownames(temp) <- NULL
      
      edna_gbif <- rbind(edna_gbif,temp)
      
    }
    
    ##still some missing
    non_classified2 <- edna_flat_fix%>%
                      filter(is.na(name))%>%
                      pull(call)%>%
                      unique()
    
    length(non_classified2) #should be 1
    
        edna_format2 <- edna_flat_fix%>%
          filter(!is.na(name))%>%
          mutate(rank=tolower(rank))%>%
          dplyr::select(name,rank,call)%>%
          group_by(call)%>%
          spread(rank,name)%>%
          ungroup()%>%
          mutate(subkingdom=NA, #these are missing from what the worms API returns so we will set to NA
                infrakingdom=NA, # setdiff(PhyloNames,edna_flat_fix2$rank%>%unique()%>%tolower())
                subphylum=NA,
                infraphylum=NA,
                superclass=NA,
                superorder=NA,
                subfamily=NA)%>%
          dplyr::select(call,all_of(PhyloNames))%>%
          left_join(.,edna_gbif%>%dplyr::select(-rank))%>%
          mutate(db="gbif")
    
    #now with worms to get the one last remaining unclassified
    edna_taxonomy_fix2 <- classification(non_classified2, db = 'worms',return_id = TRUE) 
    
    edna_flat_fix2 <- do.call("rbind", edna_taxonomy_fix2)%>%
      mutate(call = gsub("\\..*","",rownames(.)))
    
    edna_worms <- edna_flat_fix2%>%
        filter(row_number()==n())%>%
        dplyr::select(rank,call,id)%>%
        rename(TSN=id)
    
    #everything is classified so this should be 0  
    edna_flat_fix2%>%
      filter(is.na(name))%>%
      pull(call)%>%
      unique()%>%length()
    
    edna_format3 <- edna_flat_fix2%>%
      filter(!is.na(name))%>%
      mutate(rank=tolower(rank))%>%
      dplyr::select(name,rank,call)%>%
      group_by(call)%>%
      spread(rank,name)%>%
      ungroup()%>%
      mutate(subkingdom=NA, #these are missing from what the worms API returns so we will set to NA
             infrakingdom=NA, # setdiff(PhyloNames,edna_flat_fix2$rank%>%unique()%>%tolower())
             subphylum=NA,
             infraphylum=NA,
             superclass=NA,
             superorder=NA,
             order=NA,
             genus=NA,
             species=NA)%>%
      dplyr::select(call,all_of(PhyloNames))%>%
      left_join(.,edna_worms%>%dplyr::select(-rank))%>%
      mutate(db="worms")
    
### knit it all together
    master_edna <- rbind(edna_format,edna_format2,edna_format3)%>%
                   rename(taxID = TSN)%>% #rename because TSN is only relevant to itis. 
                   mutate(CO1 = NA, #place holders for a logical identifying if this marker identified this species. 
                          m18s = NA, #m for marker because you can't start a variable name with a number
                          m12s = NA,
                          m16s = NA)
    
    #now id what markers found what
    for(i in master_edna$call){
      marker <- edna_species%>%filter(call==i)%>%pull(marker)
      
      for(j in c("CO1","m12s","m16s","m18s")){ # 
      
      master_edna[master_edna$call == i,j] <- ifelse(gsub("m","",j) %in% marker,TRUE,FALSE)
      
      } #end j 'marker' in i 'call
      
      rm(marker)
    } #end i 'call'
    
## do a quick check in worms to see if we can flag the non-marine animals -- note that because worms gives back such variable taxonomic rankings, I first completed the itis taxonomy but this can be used as a 
    ## complete surrogate to the first ones. 
    marine_check <- classification(master_edna$call,"worms") #this requires some manual choice s of the accepted names
    
    marine_check_flat <- do.call("rbind", marine_check)%>%
      mutate(call = gsub("\\..*","",rownames(.)),
             probable_marine = ifelse(!is.na(rank),TRUE,FALSE))
    
    marine_edna_format <- marine_check_flat%>%
      filter(!is.na(name))%>%
      mutate(rank=tolower(rank))%>%
      dplyr::select(name,rank,call)%>%
      group_by(call)%>%
      spread(rank,name)%>%
      ungroup()%>%
      mutate(db="worms",
             CO1 = NA, #place holders for a logical identifying if this marker identified this species. 
             m18s = NA, #m for marker because you can't start a variable name with a number
             m12s = NA,
             m16s = NA)


    #now id what markers found what
    for(i in marine_edna_format$call){
      
      marker <- edna_species%>%filter(call==i)%>%pull(marker)
      
      for(j in c("CO1","m12s","m16s","m18s")){ # 
        
        marine_edna_format[marine_edna_format$call == i,j] <- ifelse(gsub("m","",j) %in% marker,TRUE,FALSE)
        
      } #end j 'marker' in i 'call
      
      rm(marker)
    } #end i 'call'
    
## Save the outputs
  #save.image("data/step2_data_clean_edna.R") #only do this when you need to
  
  write.csv(marine_edna_format,"data/marine_edna_taxonomy.csv",row.names = FALSE) #this is just the output from worms.
  write.csv(master_edna,"data/edna_taxonomy_all.csv",row.names = FALSE) #this is all of the species id'ed
  