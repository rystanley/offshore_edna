#Function to get basic taxonomic information

Classify <- function(x,input=T,rows=1,db="itis",db_check=FALSE){
  
  #x - is the latin name of the species of interest. This can be down to the lowest taxonomic ID. If you are entering a higher 
  #    taxonomic order, you can specify with a letter. (e.g., Decapoda is an Order so you can enter Decapoda O.). The letter will
  #    elicit a check to ensure the ID is appropriate. 
  
  #input - logical variable (DEFAULT TRUE) that indicates whether the information returned will include the input info, which is useful for later merging
  
  #rows - this is a variable of the taxise::classification function. By default 1 is used. 
  
  #db - this is the database of focus. db will return different results depending on the API. 'itis' is consistent and is set as default
  
  #db_check - this is a logical (DEFAULT: TRUE) that specifies as to whether a check of other databases is required incase db returns a null response. 
  
  #code dependance
  require(dplyr)
  require(taxize)
  
  input_x <- x #before any modifications
  
  #the information we want returned - fixed for now.
  PhyloNames <- c("kingdom","phylum","subphylum","class","order","family","genus","species") #Info we want
  
  #databases we want to check - fixed for now. 
  dbs <- c("itis","worms","gbif") #note that the 'bold' database API doesn't seem to work properly with 'classification' so it is note used.
  
  #taxID for the trawl they will sometimes id things to a higher taxonomic order with a trailing letter indicating 'family (F.),'Genus (G.)' 'Order (O.)' etc
  #These need to be removed but can also be used to validate
  
  taxID = NA
  
  if(grepl(" P.",x,fixed=T)){x=gsub(" P.","",x,fixed=T)%>%trimws;taxID="phylum"}
  if(grepl(" C.",x,fixed=T)){x=gsub(" C.","",x,fixed=T)%>%trimws;taxID="class"}
  if(grepl(" O.",x,fixed=T)){x=gsub(" O.","",x,fixed=T)%>%trimws;taxID="order"}
  if(grepl(" F.",x,fixed=T)){x=gsub(" F.","",x,fixed=T)%>%trimws;taxID="family"}
  if(grepl(" G.",x,fixed=T)){x=gsub(" G.","",x,fixed=T)%>%trimws;taxID="genus"}
  
  #STEP 1 - 
  
  #Run the taxize function
  output <- suppressMessages(taxize::classification(x,db=db,return_id=FALSE,rows=rows)) #default is 1 this option will only grab the first result from the database pull. 
  output2 <- as.data.frame(output[[1]],stringsAsfactors=F) 
  
  #STEP 2 check the outputs to see if anything was returned and format
  if(!is.na(output[x])){
    output3 <- as.data.frame(t(output2$name),stringsAsFactors = F)
    names(output3) <- output2$rank
    PhyloDiff <- intersect(PhyloNames,names(output3)) #what info is missing. Some databases don't return all the same info
    
    if(length(PhyloDiff)<length(PhyloNames)){
      
      ouput3 <- output3[,PhyloDiff]
      
      #what columns are missing
      MissingCols <- setdiff(PhyloNames,PhyloDiff)
      
      #create a dummy dataframe with NAs and the corresponding missing columns
      MissingDF <- data.frame(t(rep(NA,length(MissingCols))))
      colnames(MissingDF) <- MissingCols
      
      #append dummy 'NA' dataframe to the data we do have
      output3 <- cbind(output3,MissingDF)
      
      #ensure the output is in the right order
      output3 <- output3[,PhyloNames]
      
      
    } else{output3 <- output3[,PhyloNames]} 
      
   
  } #end 
  
  #If nothing is returned either assign NAs or do a check of other databases
  if(is.na(output[x])){
    
    if(db_check){ #do you want to search the other databases for a match?
      
      
      dbs <- setdiff(dbs,db) #id the databases you haven't checked already because db specified returned nothing
      
      db_comp=NULL
      
      for(i in dbs){
        
        tempout <- suppressMessages(taxize::classification(x,db=i,return_id=FALSE,rows=rows))
        tempout <- as.data.frame(tempout[[1]],stringsAsFactors=F)%>%mutate(db=i)
        
        db_comp <- rbind(db_comp,tempout)
        
      }
      
      db_comp <- db_comp%>%
        mutate(rank=tolower(rank))%>%
        filter(rank %in% PhyloNames)
      
      db_count <- db_comp%>%
        group_by(db)%>%
        summarise(count=sum(!is.na(name)))%>%
        ungroup()%>%
        data.frame()
      
      if(sum(db_count$count == max(db_count$count)) > 1){db_rep  <- db_count[which.max(db_count$count)[1],"db"]}else{db_rep <-db_count[which.max(db_count$count),"db"]}
      
      db_new <- db_comp%>%filter(db==as.character(db_rep))
    
      output3 <- data.frame(kingdom=NA,phylum=NA,subphylum=NA,class=NA,order=NA,family=NA,genus=NA,species=NA)
    
      for(i in db_new$rank){output3[1,i]=db_new[db_new$rank == i,"name"]} #fill in the information available 
      
      db=db_rep #save for the output
      
  } else {output3 <- data.frame(kingdom=NA,phylum=NA,subphylum=NA,class=NA,order=NA,family=NA,genus=NA,species=NA)}} #no alternate check
 
  
  #if we want to add the input name for merging later (default = T)
  if(input){output3$input <- input_x} 
  
  #do taxonomic check -- this got too tangly. Not sure it is needed anyway
  #if(!is.na(taxID)){if(is.na(tolower(output3[,taxID]))){message("High taxonomic classification check failed")}else{(if(tolower(output3[,taxID]) == tolower(x)) {message("Higher taxonomic classification check passed")}else{message("High taxonomic classification check failed")}
  #output3$check <- ifelse(tolower(output3[,taxID]) == tolower(x) & !is.na(tolower(output3[,taxID])),"Passed","Failed")} else {output3$check = NA}
  
  #add the database used
  output3$db <- db
    
  #return the dataframe row
  return(output3) 
  # 
}


