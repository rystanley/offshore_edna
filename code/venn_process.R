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

##read in the data
venn_data <- read.csv("data/venn_data.csv")%>%
             gather("method","species",1:2)

master_tax_edna <- read.csv("data/edna_taxonomy_all.csv") #from step2_data_clean_edna
master_tax_trawl <- read.csv("data/rv_identified.csv") #from step1_data_clean_rv

#get the species that the trawl has captured since 1970-2020
#load the data (originally pulled Dec 4 2020 - stanleyr PTRAN) --------
    
##note this has to be run in 32bit R for RODBC to work properly
    # library(Mar.datawrangling)
    # library(RODBC)
    # un = "stanleyr"
    # pw = "" #type int the password
    # 
    # get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata/2021",
    #             fn.oracle.username = un,
    #             fn.oracle.password = pw,
    #             fn.oracle.dsn = "PTRAN",
    #             force.extract = FALSE)
    # 
    # rvdat <- summarize_catches()

#download the species taxonomic IDs that can be used to match up to itis
    # Mar.datawrangling::get_data_custom(schema= "GROUNDFISH", data.dir = "R:/Science/CESD/HES_MPAGroup/Data/RVdata",
    #                                    tables = c("GSSPECIES_CODES"), usepkg = "rodbc", fn.oracle.username = un,
    #                                    fn.oracle.password = pw, fn.oracle.dsn = "PTRAN")
    # 
    # spec_tax <- GSSPECIES_CODES
    # 
    # save(spec_tax,file="data/spec_tax.RData")

#load the outputs from the above code -------------
load("data/RVdata_all.RData")
load("data/spec_tax.RData")

#clean up the data -----------

spec_tax_clean <- spec_tax%>%
    mutate(COMM = ifelse(COMM=="PINK FLABBY WHALEFISH","Pink flabby",COMM))%>% #this is the only species 'common' name with 'whale' in it. I am going to filter 'whale' out so this is a temporary change so it is not removed
    filter(!grepl("eggs",tolower(SPEC)),#eggs can't be classified
           !grepl("egg",tolower(SPEC)),
           !grepl("obsolete",tolower(SPEC)), #the nomenclature for these has been updated
           !grepl("larvae",tolower(SPEC)), #larvae can't be id. 
           !grepl("berried",tolower(SPEC)), #lobsters key'ed in as berried
           !grepl("(short)",tolower(SPEC)), #lobster key'ed in as 'small' aka short
           !grepl("POLYMAST",SPEC),# POLYMASTIA are bugs
           !SPEC %in% c("MARINE INVERTEBRATA (NS)"),
           !SPEC %in% c("ODONTOCEIT S.O.(WHALES)","CETACEAN (NS)","ODONTOCEIT S.O.(DOLPHINS)"),#not catching whales or dolphins
           !grepl("whale",tolower(COMM)), #not getting whales in the trawl
           !grepl("unid",tolower(SPEC)), # these are 'unidentified' things
           !grepl("weed",tolower(COMM)), # not counting weeds
           !grepl("bait",tolower(SPEC)), # some things are registered as bait
           !grepl("remains",tolower(SPEC)), # some things are registered as remains
           !grepl("POLYMAST",SPEC), # POLYMASTIA are bugs
           !SPEC %in% c("RESERVED","MUD","WATER","SAND","OIL(CRUDE)","PELAGIC FISH (NS)","GROUNDFISH (NS)","SHARK (NS)","IRISH MOSS (SIMILAR)","STONES AND ROCKS",
                        "MUCUS","MIXED","ORGANIC DEBRIS","FLUID","INORGANIC DEBRIS","COPEPODA S.C.,LARGE","OPEPODA S.C.,SMALL","FOREIGN ARTICLES,GARBAGE",	
                        "PARASITES,ROUND WORMS","SHRIMP-LIKE","SAND TUBE","BLOOK ARK","OPERCULUM","SQUID BEAKS","COD WORM","SEA CORALS (NS)", "COPEPODA S.C.,SMALL",
                        "PECTINIDAE SHELLS","MOLLUSCA (EMPTY)", #these could be deposited from fishing activity
                        "CLUPEIDAE/OSMERIDAE F.", #common noted as 'herring like' would bring them back to teleostei subclass to which there are many other examples in the database keyed out further. 
                        "MORINGUIDEA F.", #cannot find any reford of this family (CODE 470). Code 471 has Morginua edwardsi which is close. Will remove
                        "ODONTOCEIT S.O.(WHALES)","CETACEAN (NS)",#we don't catch whales in a trawl
                        "CRAB"))%>%
    mutate(COMM = ifelse(COMM=="Pink flabby","PINK FLABBY WHALEFISH",COMM),#fix the temporary miss-label
           SPEC = gsub("PURSE ","",SPEC),#Some things are labeled as purses for sharks and skates. 
           format_id = 1:n())%>% #this is an identifier we can use to bring together to re-assemble the data 
    filter(SPEC !=  "DEEPWATER SKATE",#cannot find a probable common name that is associated with this common name #'crab' too unspecific ### HAVE NO IDEA WHY THIS WON'T work when nested in the !spec list above
           SPEC != "COELENTERATA P.") #Not a monophyletic group https://www.marinespecies.org/aphia.php?p=taxdetails&id=152230


#which TSN values have more than 1 entry -------------
dup_tsn <- spec_tax_clean%>%
            filter(!is.na(TSN))%>%
            group_by(TSN)%>%
            summarise(count=n())%>%
            ungroup()%>%
            filter(count>1)%>%
            arrange(-count)%>%
            data.frame()%>%
            left_join(.,spec_tax_clean)%>%
            dplyr::select(names(spec_tax_clean))

dup_tsn_fixed <- dup_tsn%>%
                  mutate(SPEC = ifelse(SPEC == "RUVETTUS PRETIOSUS","Ruvettus pretiosus",SPEC), #"RUVETTUS PRETIOSUS" also listed as 159907 should have ITIS 172364 -- aslo spelling error should be "Ruvettus pretiosus"
                       TSN = ifelse(SPEC == "Ruvettus pretiosus",172364,TSN), # "ISTIOPHORUS PLATYPTERUS" also listed as 159907 should have ITIS 172488
                       TSN = ifelse(SPEC == "ISTIOPHORUS PLATYPTERUS",172488,TSN), #wrong TSN
                       SPEC = ifelse(TSN == 48506, "Polymastia G.",SPEC),# the Polymastia is a genus according to itis. 
                       SPEC = ifelse(TSN == 51938, "Anthozoa C.",SPEC), # class Anthozoa
                       SPEC = ifelse(TSN == 52431,"Zoantharia F.",SPEC), # family Zoantharia
                       SPEC = ifelse(TSN == 57411,"NEMERTEA P.",SPEC),
                       SPEC = ifelse(TSN == 64358,"POLYCHAETA C.",SPEC), #class Polychaeta
                       SPEC = ifelse(TSN == 73297,"Nucella lapillus",SPEC),
                       SPEC = ifelse(TSN == 93294,"Amphipoda O.",SPEC), #class Amphipoda
                       SPEC = ifelse(TSN == 96797,"Lebbeus microceros",SPEC),
                       SPEC = ifelse(TSN == 157356,"Gorgonocephalidae F.",SPEC), #Family Gorgonocephalidae
                       SPEC = ifelse(TSN == 158262,"STEREODERMA UNISEMITA",SPEC),
                       TSN = ifelse(SPEC == "PRIONACE GLAUCA",160424,TSN), #wrong serial number for blue shark
                       SPEC = ifelse(TSN == 162272, "Photostomias guernei",SPEC),
                       SPEC = ifelse(TSN == 168181,"Heteropriacanthus cruentatus",SPEC), #both species names COOKEOLUS BOOPS and PRIACANTHUS CRUENTATUS are invalid on itis
                       TSN = ifelse(SPEC == "PRIACANTHUS ARENATUS",168178,TSN), #wrong TSN for the bigeye
                       SPEC = ifelse(TSN == 168684,"Selene setapinnis",SPEC), #VOMER SETAPINNIS is the invalid name
                       TSN = ifelse(SPEC == "CHIASMODON BOLANGERI",630464,TSN), #CHIASMODON BOLANGERI also called the black swallower so similar to CHIASMODON NIGER
                       TSN = ifelse(SPEC == "SCOMBEROMORUS MACULATUS",172436,TSN), #wrong tsn assigned 
                       SPEC = ifelse(TSN == 173392,"Diodon holocanthus",SPEC),
                       SPEC = ifelse(TSN == 182909,"Scopeloberyx robustus",SPEC), # SCOPELOBERYX NIGRESCENS not valid. Same species as robustus
                       SPEC = ifelse(TSN == 564219,"Rajella bigelowi",SPEC),
                       TSN = ifelse(SPEC == "MELANOSTOMIAS SPILORHYNCHUS",660867,TSN), #660867 MELANOSTOMIAS SPILORHYNCHUS is a valid species ID as is MELANOSTOMIAS BARTONBEANI
                       SPEC = ifelse(TSN == 622417, "ARISTOSTOMIAS LUNIFER",SPEC), # ARISTOSTOMIAS PHOTODACTYLUS not valid
                       SPEC = ifelse(TSN == 630981,"Lycodes terraenovae",SPEC), # LYCODES ATLANTICUS not valid
                       SPEC = ifelse(TSN == 690523,"LINOPHRYNE MACRODON",SPEC), # LINOPHRYNE BREVIBARBIS not valid
                       SPEC = ifelse(TSN == 963108,"Polysiphonia G.",SPEC),
                       TSN = ifelse(TSN == 	963108,13440,TSN)) #this had the wrong TSN - which assigned to a bacteria instead of the Polysiphonia genus. 


#fix the missing TSNs -----------
miss_TSN <- spec_tax_clean%>%
            filter(is.na(TSN),
                   is.na(APHIAID))%>%
            mutate(SPEC_org = SPEC, #this is the original 'spec' associated with each code
                   TSN = ifelse(SPEC == "APOGON SELLICAUDE",168246,TSN), # miss-spelling
                   SPEC = ifelse(SPEC == "APOGON SELLICAUDE","Apogon sellicauda",SPEC),
                   TSN = ifelse(SPEC == "MELANOSTOMIATIDAE (STOMIATIDAE) F.",622318,TSN), #misspelling Melanostomiinae 
                   SPEC = ifelse(SPEC == "MELANOSTOMIATIDAE (STOMIATIDAE) F.","Melanostomiinae",SPEC),
                   TSN = ifelse(SPEC == "STOMIATIFORMES (ORDER)",553138,TSN), #misspelling 	Stomiatiformes
                   SPEC = ifelse(SPEC == "STOMIATIFORMES (ORDER)","Stomiatiformes",SPEC),
                   TSN = ifelse(SPEC == "ASCIDIA SP. LARVAL",159174,TSN),
                   SPEC = ifelse(SPEC == "ASCIDIA SP. LARVAL","Ascidia",SPEC),
                   TSN = ifelse(SPEC == "ASCIDIA SP. ADULT",159174,TSN),
                   SPEC = ifelse(SPEC == "ASCIDIA SP. ADULT","Ascidia",SPEC),
                   TSN = ifelse(SPEC == "BRYOZOANS ECTOPROCTA P.",155469,TSN), # Ectoprocta is a bryozoan
                   SPEC = ifelse(SPEC == "BRYOZOANS ECTOPROCTA P.","Bryozoa",SPEC),
                   TSN = ifelse(SPEC == "CALAPPA MEGALOPS",98341,TSN), # there is no official species named magalops in the clappa genus according to itis and wiki https://en.wikipedia.org/wiki/Calappa_(crab)
                   SPEC = ifelse(SPEC == "CALAPPA MEGALOPS","Calappa",SPEC),
                   TSN = ifelse(SPEC == "CTENOPHORES,COELENTERATES,PORIFERA P.",118845,TSN), # Colenterata and ctenophores are both ctenophores not sure where the porifora (sponges) came from but will assign a 118845 for Ctenophora
                   SPEC = ifelse(SPEC == "CTENOPHORES,COELENTERATES,PORIFERA P.","Ctenophora",SPEC),
                   TSN = ifelse(SPEC == "TEALIA  FELINA",573747,TSN), # accepted name corrected based on worms 
                   SPEC = ifelse(SPEC == "TEALIA  FELINA","Urticina felina",SPEC),
                   TSN = ifelse(SPEC == "SAURIDA CARRIBBEA",162410,TSN), # miss-spelling
                   SPEC = ifelse(SPEC == "SAURIDA CARRIBBEA","Saurida caribbaea",SPEC),
                   TSN = ifelse(SPEC == "LAMPETRA LAMOTTEI",914061,TSN), # this classification is not longer accepted. The new accepted name is is Lethenteron appendix according to itis
                   SPEC = ifelse(SPEC == "LAMPETRA LAMOTTEI","Lethenteron appendix",SPEC),
                   TSN = ifelse(SPEC == "NOTEMIGOMUS CRYSOLEUCAS",163368,TSN), # miss-spelling
                   SPEC = ifelse(SPEC == "NOTEMIGOMUS CRYSOLEUCAS","Notemigonus crysoleucas",SPEC),
                   TSN = ifelse(SPEC == "NEONESTHES CAPERSIS",622368,TSN), # miss-spelling
                   SPEC = ifelse(SPEC == "NEONESTHES CAPERSIS","Neonesthes capensis",SPEC),
                   TSN = ifelse(SPEC == "TRACHIPTEROIDEI",638685,TSN), # TRACHIPTEROIDEI is misspelled from Trachipteroidei which is no-longer valid should be Lampriformes according to itis
                   SPEC = ifelse(SPEC == "TRACHIPTEROIDEI","Lampriformes",SPEC),
                   TSN = ifelse(SPEC == "ASTRONESTHES  SP.(RICHARDSONI?)",622316,TSN), # miss-spelling
                   SPEC = ifelse(SPEC == "ASTRONESTHES  SP.(RICHARDSONI?)","Astronesthinae",SPEC),
                   TSN = ifelse(SPEC == "SEARSIID SP.",623294,TSN),
                   SPEC = ifelse(SPEC == "SEARSIID SP.","Searsia",SPEC),
                   TSN = ifelse(SPEC == "LEPTOSTOMIAS HAPLOCAVLUS",622467,TSN),
                   SPEC = ifelse(SPEC == "LEPTOSTOMIAS HAPLOCAVLUS","Leptostomias haplocaulus",SPEC), #spelling error v instead of u
                   TSN = ifelse(SPEC == "RADIOLARIDA O.",46088,TSN),
                   SPEC = ifelse(SPEC == "RADIOLARIDA O.","Radiolaria",SPEC), # there is no 'd' in Radiolaria
                   TSN = ifelse(SPEC == "PROTOCHORDATA SP.",158853,TSN), #Protochorodates are 'lower chorodates' refering to two subphya that lack craniums - urochorodates and chephalocorodates - given the frequency of tunicates and ascidians this is going to be classified as a 'urochorodate'
                   SPEC = ifelse(SPEC == "PROTOCHORDATA SP.","Urochordata",SPEC),
                   TSN = ifelse(SPEC == "ASCIDIA",158854,TSN),
                   SPEC = ifelse(SPEC == "ASCIDIA","Ascidiacea",SPEC), #proper spelling
                   TSN = ifelse(SPEC == "THYSANOESSA ACUTIFRONS",95583,TSN),
                   SPEC = ifelse(SPEC == "THYSANOESSA ACUTIFRONS","Thysanopoda acutifrons",SPEC), #spelling error in the genus name
                   TSN = ifelse(SPEC == "EUPHAUSIACEA/MYSIDACEA O.",95496,TSN),
                   SPEC = ifelse(SPEC == "EUPHAUSIACEA/MYSIDACEA O.","Euphausiacea",SPEC), #Mysidacea is no longer valid
                   TSN = ifelse(SPEC == "ERYTHROPS ERYTHROPTHALMA",90185,TSN),
                   SPEC = ifelse(SPEC == "ERYTHROPS ERYTHROPTHALMA","Erythrops erythrophthalma",SPEC), #seems to be a spelling error
                   TSN = ifelse(SPEC == "HYPERIA OCULATA",95109,TSN),
                   SPEC = ifelse(SPEC == "HYPERIA OCULATA","Hyperia",SPEC),# there is no species in worms or itis under the genus 'Hyperia' with the name Oculata. So this is assigned to the genus
                   TSN = ifelse(SPEC == "ORCHOMONELLA SP.",94457,TSN),
                   SPEC = ifelse(SPEC == "ORCHOMONELLA SP.","Orchomenella",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ORCHOMONELLA PINGUIS",94460,TSN),
                   SPEC = ifelse(SPEC == "ORCHOMONELLA PINGUIS","Orchomenella pinguis",SPEC), #spelling error
                   TSN = ifelse(SPEC == "HIPPOMEDON DENTERULATUS",94295,TSN),
                   SPEC = ifelse(SPEC == "HIPPOMEDON DENTERULATUS","Hippomedon denticulatus",SPEC), #spelling error
                   TSN = ifelse(SPEC == "COPEPODA,NAUPLII",85257,TSN),
                   SPEC = ifelse(SPEC == "COPEPODA,NAUPLII","Copepoda",SPEC), #Naupulii are just a designation of a larval stage
                   TSN = ifelse(SPEC == "EDOTEA TRILOBA",544186,TSN),
                   SPEC = ifelse(SPEC == "EDOTEA TRILOBA","Edotia triloba",SPEC), #updated spelling
                   TSN = ifelse(SPEC == "PLEUSTES PANOPLA",94790,TSN),
                   SPEC = ifelse(SPEC == "PLEUSTES PANOPLA","Pleustes panoplus",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ORCHOMONELLA MINUTA",94458,TSN),
                   SPEC = ifelse(SPEC == "ORCHOMONELLA MINUTA","Orchomenella minuta",SPEC), #spelling error
                   TSN = ifelse(SPEC == "HELIRAGES FULVOCINCTUS",93549,TSN),
                   SPEC = ifelse(SPEC == "HELIRAGES FULVOCINCTUS","Halirages fulvocinctus",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ISCHYROCEROS ANGUIPES",94153,TSN),
                   SPEC = ifelse(SPEC == "ISCHYROCEROS ANGUIPES","Ischyrocerus anguipes",SPEC), #spelling error
                   TSN = ifelse(SPEC == "PAROEDICEROS PROPINGUIS",94562,TSN),
                   SPEC = ifelse(SPEC == "PAROEDICEROS PROPINGUIS","Paroediceros propinquus",SPEC),
                   TSN = ifelse(SPEC == "PAROEDICEROS LONGIMANUS",94581,TSN),
                   SPEC = ifelse(SPEC == "PAROEDICEROS LONGIMANUS","Perioculodes longimanus",SPEC), #note that PAROEDICEROS LONGIMANUS does not exist in ITIS or worms and only comes up in searches associatd with the survey. Assume these are miss classified to the genus paroediceros (which was itself misspelled). Distribution of Perioculodes longimanus includes the sruvey region of 
                   TSN = ifelse(SPEC == "RACHOTROPIS OCULATA",93738,TSN),
                   SPEC = ifelse(SPEC == "RACHOTROPIS OCULATA","Rhachotropis oculata",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ISCHYROCEROS SP.",94151,TSN),
                   SPEC = ifelse(SPEC == "ISCHYROCEROS SP.","Ischyrocerus",SPEC), #spelling error u for o. Can't find anything with a OS spelling though there are some terrestrial bugs with a similar name. 
                   TSN = ifelse(SPEC == "NYMPH0N LONGITARSE",83550,TSN),
                   SPEC = ifelse(SPEC == "NYMPH0N LONGITARSE","Nymphon longitarse",SPEC), #no TSN associated in original RV db
                   TSN = ifelse(SPEC == "PARACALANUS,CLAUSOCALANUS,CALOCALANUS SP.",85333,TSN), 
                   SPEC = ifelse(SPEC == "PARACALANUS,CLAUSOCALANUS,CALOCALANUS SP.","Copepoda",SPEC), #all different taxonomies of the subclass copepoda. Since the'sp' was associated with the genus calocalanus we will use that for this term
                   TSN = ifelse(SPEC == "BATHYPTYPHLOPS MARIONAE",162457,TSN),
                   SPEC = ifelse(SPEC == "BATHYPTYPHLOPS MARIONAE","Bathytyphlops marionae",SPEC), #spelling error
                   TSN = ifelse(SPEC == "NUCULA  TENUIS",567536,TSN),
                   SPEC = ifelse(SPEC == "NUCULA  TENUIS","Ennucula tenuis",SPEC), # accepted name is now Ennucula tenuis or the 'smooth nutclam' according to itis. This taxonomy was formally referred to as Nucula tenuis
                   TSN = ifelse(SPEC == "TEUTHOIDEA O.",82367,TSN),
                   SPEC = ifelse(SPEC == "TEUTHOIDEA O.","Teuthida",SPEC), # spelling error
                   TSN = ifelse(SPEC == "LOLIGINIDAE,OMMASTREPHIDAE F.",555706,TSN),
                   SPEC = ifelse(SPEC == "LOLIGINIDAE,OMMASTREPHIDAE F.","Decabrachia",SPEC), # note that these are two families under the Dechabrancia superorder. Note that WORMS and ITIS have different names and even trees here (different subordinate classes). On WORMS it they are referred to Decapodiformes
                   APHIAID = ifelse(SPEC == "Decabrachia",325342,APHIAID), #updating this aphia to match the slight disconnect between ITIS and worms
                   TSN = ifelse(SPEC == "OMMASTREPHES PTEROPSUS",82530,TSN),
                   SPEC = ifelse(SPEC == "OMMASTREPHES PTEROPSUS","Sthenoteuthis pteropus",SPEC), #spelled incorrectly and no-longer valid. updated now. 
                   TSN = ifelse(SPEC == "ENOPLOTEUTHINAE S.F.",82397,TSN),
                   SPEC = ifelse(SPEC == "ENOPLOTEUTHINAE S.F.","Enoploteuthidae",SPEC), #spelling error
                   
                   TSN = ifelse(SPEC == "BATHYNECTES MARAVIGNA (SUPERBUS)",NA,TSN), #note the 'maravigna' species is not in the itis database so we will leave as NA the genus bathynectes is in itis with a TSN  - 98692
                   SPEC = ifelse(SPEC == "BATHYNECTES MARAVIGNA (SUPERBUS)","Bathynectes maravigna",SPEC),
                   APHIAID = ifelse(SPEC == "Bathynectes maravigna",107377,APHIAID), 
                  
                    TSN = ifelse(SPEC == "PENTAGONASTER TOSIA",157035,TSN),
                   SPEC = ifelse(SPEC == "PENTAGONASTER TOSIA","Pentagonaster",SPEC), #Tosia is not a species, it is a former family name (as far as I can tell). Assigned to the family level now
                   APHIAID = ifelse(SPEC == "Pentagonaster",123302,APHIAID), 
                   TSN = ifelse(SPEC == "CALATHURA BRANCHIATA",92212,TSN),
                   SPEC = ifelse(SPEC == "CALATHURA BRANCHIATA","Calathura brachiata",SPEC), #spelling error
                   TSN = ifelse(SPEC == "POLYCHAETA C.,SMALL",64358,TSN),
                   SPEC = ifelse(SPEC == "POLYCHAETA C.,SMALL","Polychaeta",SPEC), #can't id based on size
                   TSN = ifelse(SPEC == "LUMBRINEREIDAE F.",66335,TSN),
                   SPEC = ifelse(SPEC == "LUMBRINEREIDAE F.","Lumbrineridae",SPEC),# spelling error
                   TSN = ifelse(SPEC == "MELINNIA ELIZABETHAE",67765,TSN),
                   SPEC = ifelse(SPEC == "MELINNIA ELIZABETHAE","Melinna elisabethae",SPEC), # spelling error
                   TSN = ifelse(SPEC == "OPHELIA ACUMINUTA",67391,TSN),
                   SPEC = ifelse(SPEC == "OPHELIA ACUMINUTA","Ophelina acuminata",SPEC), #spelling error
                   TSN = ifelse(SPEC == "PHERUSA PHERUSA",67241,TSN),
                   SPEC = ifelse(SPEC == "PHERUSA PHERUSA","Pherusa",SPEC), # there is a family called Pherusa but there is no species called pherusa. the common name is noted as flabelligerid worm which can be associted with all worms in the family pherusa. keeping this a the family level 
                   TSN = ifelse(SPEC == "MARGARITES GROENLANDICA",69897,TSN),
                   SPEC = ifelse(SPEC == "MARGARITES GROENLANDICA","Margarites groenlandicus",SPEC),#spelling error on the unaccepted name. This is now the accepted nomenclature
                   TSN = ifelse(SPEC == "NASSARIIDAE OR THAISIDAE F.",74102,TSN),
                   SPEC = ifelse(SPEC == "NASSARIIDAE OR THAISIDAE F.","Nassariidae",SPEC), # there is not Thaisidae family on worms or itis
                   APHIAID = ifelse(SPEC == "Nassariidae",151,APHIAID), #low aphia!
                   TSN = ifelse(SPEC == "MARGARITES HELICINA",69879,TSN),
                   SPEC = ifelse(SPEC == "MARGARITES HELICINA","Margarites (Margarites) helicinus",SPEC),
                   TSN = ifelse(SPEC == "PITAR MORRHUANA",81501,TSN),
                   SPEC = ifelse(SPEC == "PITAR MORRHUANA","Pitar morrhuanus",SPEC), #new spelling - updated
                   
                   TSN = ifelse(SPEC == "STROMBUS AND BUSYCON SP.",NA,TSN), #see aphia
                   SPEC = ifelse(SPEC == "STROMBUS AND BUSYCON SP.","Caenogastropoda",SPEC), # both are part of the Caenogastropoda  subclass (which is not key'ed out on ITIS). However their genus is correct. Thes can be added in after the fact for the comparison, though the code is only commonly associated to this point
                   APHIAID = ifelse(SPEC == "Caenogastropoda",224570,APHIAID),
                   
                   TSN = ifelse(SPEC == "SCOMBROIDEI (SUBORDER)",172353,TSN),
                   SPEC = ifelse(SPEC == "SCOMBROIDEI (SUBORDER)","Scombroidei",SPEC),
                   TSN = ifelse(SPEC == "GADOIDEI S.O.",164665,TSN),
                   SPEC = ifelse(SPEC == "GADOIDEI S.O.","Gadiformes",SPEC), #gadoids under the order Gadiformes
                   TSN = ifelse(SPEC == "CENTROBRANCHUS NIGRO-OCELLATUS",162732,TSN),
                   SPEC = ifelse(SPEC == "CENTROBRANCHUS NIGRO-OCELLATUS","Centrobranchus nigroocellatus",SPEC), #hyphen between the 'o's' removed
                   TSN = ifelse(SPEC == "LOBIANCHIA GEMELLERII",623876,TSN),
                   SPEC = ifelse(SPEC == "LOBIANCHIA GEMELLERII","Lobianchia gemellarii",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ENCHELYOPUS/UROPHYCIS SP.","164748",TSN),
                   SPEC = ifelse(SPEC == "ENCHELYOPUS/UROPHYCIS SP.","Enchelyopus cimbrius",SPEC), #Urophicis is not a rockling hake - the common name here was given as 'rockling hake' so this was assigned the cimbrius species 164748 http://gma.org/fogm/Enchelyopus_cimbrius.htm
                   TSN = ifelse(SPEC == "ASTROTECTEN DUPLICATUS","156903",TSN),
                   SPEC = ifelse(SPEC == "ASTROTECTEN DUPLICATUS","Astropecten duplicatus",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ARABACIA SP.",157905,TSN),
                   SPEC = ifelse(SPEC == "ARABACIA SP.","Arbacia",SPEC), #spelling error
                   TSN = ifelse(SPEC == "PSOLUSES,THYONES,ETC. (NS)",158140,TSN),
                   SPEC = ifelse(SPEC == "PSOLUSES,THYONES,ETC. (NS)","Holothuroidea",SPEC), #class that covers seacucumbers which this is ultimately referring too. 
                   
                   TSN = ifelse(SPEC == "PORROCAECUM DECIPIENS",NA,TSN), # not ITIS id for this species. 
                   SPEC = ifelse(SPEC == "PORROCAECUM DECIPIENS","Pseudoterranova decipiens",SPEC), #this is the new valid taxonomy according to WORMS
                   APHIAID = ifelse(SPEC == "Pseudoterranova decipiens",123078,APHIAID),
                   
                   TSN = ifelse(SPEC == "BONAPARTIA PEDIOLOTA",162197,TSN),
                   SPEC = ifelse(SPEC == "BONAPARTIA PEDIOLOTA","Bonapartia pedaliota",SPEC), #spelling error
                   
                   TSN = ifelse(SPEC == "MELAMPHAEIDAE",166083,TSN),
                   SPEC = ifelse(SPEC == "MELAMPHAEIDAE","Beryciformes",SPEC), # these are a family of Beryciformes but both worms and itis don't have a record of that family name, though there are several references online https://researcharchive.calacademy.org/research/ichthyology/catalog/SpeciesByFamily.asp#Melamphaidae
                   
                   TSN = ifelse(SPEC == "ASCONEMA FOLIATA",NA,TSN), #no species listed under the genus Asconema in ITIS (#659654)
                   SPEC = ifelse(SPEC == "ASCONEMA FOLIATA","Asconema foliata",SPEC),
                   APHIAID = ifelse(SPEC == "Asconema foliata",1496787,APHIAID),
                   
                   TSN = ifelse(SPEC == "MONOMITOPIS AGASSIZII",622859,TSN),
                   SPEC = ifelse(SPEC == "MONOMITOPIS AGASSIZII","Monomitopus agassizii",SPEC), #spelling error
                   
                   TSN = ifelse(SPEC == "SELACHII (CHONDRICHTHYES) C.",NA,TSN),
                   SPEC = ifelse(SPEC == "SELACHII (CHONDRICHTHYES) C.","Chondrichthyes",SPEC), #is an infraclass of the class Elasmobranchii
                   APHIAID = ifelse(SPEC == "Chondrichthyes",368408,APHIAID),
                   
                   TSN = ifelse(SPEC == "SCORPAENIFORMES (ORDER)",166702,TSN),
                   SPEC = ifelse(SPEC == "SCORPAENIFORMES (ORDER)","Scorpaeniformes",SPEC),
                   TSN = ifelse(SPEC == "EVERMANELLA BALBOA",162542,TSN),
                   SPEC = ifelse(SPEC == "EVERMANELLA BALBOA","Evermannella balbo",SPEC), #spelling error
                   TSN = ifelse(SPEC == "PSEUDOPHICHTHUS SPLENDENS",635810,TSN),
                   SPEC = ifelse(SPEC == "PSEUDOPHICHTHUS SPLENDENS","Pseudophichthys splendens",SPEC), #spelling error
                   TSN = ifelse(SPEC == "DYSOMMINAE S.F.",161285,TSN),
                   SPEC = ifelse(SPEC == "DYSOMMINAE S.F.","Dysommina",SPEC), #spelling error this should be the genus Dysommina. The subfamily was previously known as Dysommatinae (so spelling difference but close) now referred to as Ilyophinae. 
                   TSN = ifelse(SPEC == "MORGINUA EDWARDSI",161138,TSN),
                   SPEC = ifelse(SPEC == "MORGINUA EDWARDSI","Moringua edwardsi ",SPEC), #spelling error
                   TSN = ifelse(SPEC == "EVERMANELLA INDICA",162544,TSN),
                   SPEC = ifelse(SPEC == "EVERMANELLA INDICA","Evermannella indica",SPEC), #spelling error
                   TSN = ifelse(SPEC == "GONATUS STEENSTRUPII",205721,TSN),
                   SPEC = ifelse(SPEC == "GONATUS STEENSTRUPII","Gonatus steenstrupi",SPEC), #spelling error
                   TSN = ifelse(SPEC == "MONOSTRTOMA SPP.",6481,TSN),
                   SPEC = ifelse(SPEC == "MONOSTRTOMA SPP.","Monostroma",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ANTIPATHARIAN O.",51940,TSN),
                   SPEC = ifelse(SPEC == "ANTIPATHARIAN O.","Antipatharia",SPEC),
                   TSN = ifelse(SPEC == "HETROPODIDAE",155757,TSN),
                   SPEC = ifelse(SPEC == "HETROPODIDAE","Heteroporidae",SPEC), #spelling error - staghorn corals (northern) - genus heteropora associated with the family Heteroporidae. 
                   TSN = ifelse(SPEC == "OCTOPODA, CIRRATA",555711,TSN),
                   SPEC = ifelse(SPEC == "OCTOPODA, CIRRATA","Cirrata",SPEC), #suborder cirrata this was not spelled correctly on itis this is called 'cirrina' with the same child families 
                   APHIAID = ifelse(SPEC == "Cirrata",11732,APHIAID),
                   TSN = ifelse(SPEC == "LEPODOCHELYS KEMPII",551770,TSN),
                   SPEC = ifelse(SPEC == "LEPODOCHELYS KEMPII","Lepidochelys kempii",SPEC), #spelling error
                   TSN = ifelse(SPEC == "LEPODOCHELYS OLIVACEA",173840,TSN),
                   SPEC = ifelse(SPEC == "LEPODOCHELYS OLIVACEA","Lepidochelys olivacea",SPEC), #spelling error
                   TSN = ifelse(SPEC == "CRAB(ANOMURA)",97698,TSN),
                   SPEC = ifelse(SPEC == "CRAB(ANOMURA)","Anomura",SPEC),
                   TSN = ifelse(SPEC == "CRAB (OVALIPES SP)",98710,TSN),
                   SPEC = ifelse(SPEC == "CRAB (OVALIPES SP)","Ovalipes",SPEC),
                   TSN = ifelse(SPEC == "STYPOCAULON SCOPARIA",11125,TSN),
                   SPEC = ifelse(SPEC == "STYPOCAULON SCOPARIA","Halopteris scoparia",SPEC), #was a spelling error from Stypocaulon scoparium which is no longer valid. Valid taxonomy replaced 
                   TSN = ifelse(SPEC == "HALOSAUROPSIS MACHAROCHIR",161672,TSN),
                   SPEC = ifelse(SPEC == "HALOSAUROPSIS MACHAROCHIR","Halosauropsis macrochir ",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ONERIDES SP.",164629,TSN),
                   SPEC = ifelse(SPEC == "ONERIDES SP.","Oneirodes",SPEC), #assume this is a spelling error
                   TSN = ifelse(SPEC == "NOTOCANTHUS BONAPARTE",635871,TSN),
                   SPEC = ifelse(SPEC == "NOTOCANTHUS BONAPARTE","Notacanthus bonaparte",SPEC), #spelling error
                   TSN = ifelse(SPEC == "PSEUDOSCOPELUS ALTIPINNUS",171094,TSN),
                   SPEC = ifelse(SPEC == "PSEUDOSCOPELUS ALTIPINNUS","Pseudoscopelus altipinnis",SPEC), #spelling error
                   TSN = ifelse(SPEC == "ACANTHEPHYRA EXEMIA",96120,TSN),
                   SPEC = ifelse(SPEC == "ACANTHEPHYRA EXEMIA","Acanthephyra eximia",SPEC), #spelling error
                   TSN = ifelse(SPEC == "NOTOMASTUS ROBUSTUS",67423,TSN),
                   SPEC = ifelse(SPEC == "NOTOMASTUS ROBUSTUS","Notomastus",SPEC), #Both Notomastus Robustus and Elegans were not found in ITIS or Worms. There is a genus Notomastus which I assume they belong. They will be assigned to that level
                   TSN = ifelse(SPEC == "NOTOMASTUS ELEGANS",67423,TSN),
                   SPEC = ifelse(SPEC == "NOTOMASTUS ELEGANS","Notomastus",SPEC), 
                   TSN = ifelse(SPEC == "ARISTEIDAE F",96049,TSN),
                   SPEC = ifelse(SPEC == "ARISTEIDAE F","Aristeidae",SPEC),
                   TSN = ifelse(SPEC == "HALIPTERUS (BALTICINA) SP.",52406,TSN),
                   SPEC = ifelse(SPEC == "HALIPTERUS (BALTICINA) SP.","Balticina",SPEC),
                   
                   #purses
                   TSN = ifelse(SPEC == "BARNDOOR SKATE",564139,TSN),
                   SPEC = ifelse(SPEC == "BARNDOOR SKATE","Dipturus laevis",SPEC),
                   TSN = ifelse(SPEC == "THORNY SKATE",564149,TSN),
                   SPEC = ifelse(SPEC == "THORNY SKATE","Amblyraja radiata",SPEC),
                   TSN = ifelse(SPEC == "SMOOTH SKATE",564151,TSN),
                   SPEC = ifelse(SPEC == "SMOOTH SKATE","Malacoraja senta",SPEC),
                   TSN = ifelse(SPEC == "LITTLE SKATE",564130,TSN),
                   SPEC = ifelse(SPEC == "LITTLE SKATE","Leucoraja erinacea",SPEC),
                   TSN = ifelse(SPEC == "WINTER SKATE",564145,TSN),
                   SPEC = ifelse(SPEC == "WINTER SKATE","Leucoraja ocellata",SPEC),
                   TSN = ifelse(SPEC == "SPINYTAIL SKATE",160932,TSN),
                   SPEC = ifelse(SPEC == "SPINYTAIL SKATE","Bathyraja spinicauda",SPEC),
                   
                   TSN = ifelse(SPEC == "BRIER SKATE",NA,TSN),
                   SPEC = ifelse(SPEC == "BRIER SKATE","Rostroraja eglanteria",SPEC), # Rostroraja eglanteria or 'clearnose skate' is not listed in itis.
                   APHIAID = ifelse(SPEC == "Rostroraja eglanteria",1460141,APHIAID),
                   
                   TSN = ifelse(SPEC == "ROUND SKATE",564135,TSN),
                   SPEC = ifelse(SPEC == "ROUND SKATE","Rajella fyllae",SPEC),
                   TSN = ifelse(SPEC == "SOFT SKATE",564152,TSN),
                   SPEC = ifelse(SPEC == "SOFT SKATE","Malacoraja spinacidermis",SPEC),
                   TSN = ifelse(SPEC == "ARCTIC SKATE",564138,TSN),
                   SPEC = ifelse(SPEC == "ARCTIC SKATE","Amblyraja hyperborea",SPEC),
                   TSN = ifelse(SPEC == "ABYSSAL SKATE",564125,TSN),
                   SPEC = ifelse(SPEC == "ABYSSAL SKATE","Rajella bathyphila ",SPEC),
                   
                   TSN = ifelse(SPEC == "WHITE SKATE",NA,TSN), #not on itis
                   SPEC = ifelse(SPEC == "WHITE SKATE","Rajella lintea",SPEC), #also the 'linen skate'
                   APHIAID = ifelse(SPEC == "Rajella lintea",1019159,APHIAID),
                   
                   TSN = ifelse(SPEC == "RAJA BIGELOW SKATE",564219,TSN),
                   SPEC = ifelse(SPEC == "RAJA BIGELOW SKATE","Rajella bigelowi",SPEC),
                   TSN = ifelse(SPEC == "ROSETTE SKATE",564136,TSN),
                   SPEC = ifelse(SPEC == "ROSETTE SKATE","Leucoraja garmani",SPEC),
                   TSN = ifelse(SPEC == "FRECKLED SKATE",564140,TSN),
                   SPEC = ifelse(SPEC == "FRECKLED SKATE","Leucoraja lentiginosa",SPEC),
                   TSN = ifelse(SPEC == "JENSENS SKATE",564207,TSN),
                   SPEC = ifelse(SPEC == "JENSENS SKATE","Amblyraja jenseni",SPEC), #aka the shorttail skate
                   
                   TSN = ifelse(SPEC == "HENRICA SP.",157152,TSN),
                   SPEC = ifelse(SPEC == "HENRICA SP.","Henricia",SPEC) #assume spelling error for a genus of sea stars henricia
                  
                   )  

#get the data with some sort of taxonomic identifier  ----------------

dup_ids <- dup_tsn_fixed%>%pull(format_id)
miss_TSN_ids <- miss_TSN%>%pull(format_id)

clean_data_itis <- rbind(
                         spec_tax_clean%>%
                             filter(!format_id %in% c(dup_ids,miss_TSN_ids))%>%
                             filter(!is.na(TSN)), #the ones which weren't cleaned up
                         dup_tsn_fixed, #the ones that were cleaned up
                         miss_TSN%>%filter(!is.na(TSN))%>%dplyr::select(names(spec_tax_clean))
                         )

clean_data_aphia <- rbind(
                            spec_tax_clean%>%
                                filter(!format_id %in% c(dup_ids,miss_TSN_ids))%>%
                                filter(is.na(TSN)),
                            miss_TSN%>%
                                filter(is.na(TSN))%>%
                                dplyr::select(names(spec_tax_clean))
                          )

#Run the classification from itis and those without taxinomic serial numbers run from aphia. ----------------
itis_taxonomy <- pblapply(clean_data_itis%>%pull(TSN),Classify)%>%
                 do.call("rbind",.)%>%
                 left_join(.,clean_data_itis%>%select(SPEC,COMM,CODE,format_id,TSN)%>%rename(id=TSN)%>%distinct(id,.keep_all=T))


aphia_taxonomy <- pblapply(clean_data_aphia%>%pull(APHIAID),FUN=function(x) Classify(x,db="worms"))%>%
                  do.call("rbind",.)%>%
                  left_join(.,clean_data_aphia%>%select(SPEC,COMM,CODE,format_id,APHIAID)%>%rename(id=APHIAID)%>%distinct(id,.keep_all=T))

# ##now use the cleaned itis species names to do an aphia search 
# 
        PhyloNames <- c("kingdom","subkingdom","infrakingdom","phylum","subphylum","infraphylum","superclass",
                        "class","superorder","order","family","subfamily","genus","species")

        # itis_names <- itis_taxonomy%>%
        #               select(PhyloNames,format_id)%>%
        #               gather(key=PhyloNames,value="value",-format_id)%>%
        #               group_by(format_id)%>%
        #               summarise(lowest_id = name_extract(value))%>%
        #               ungroup()%>%
        #               left_join(.,itis_taxonomy)
        # 

#save intermediate outputs
save(itis_taxonomy,file="output/itis_taxonomy.csv")
save(aphia_taxonomy,file="output/worms_taxonomy.csv")

### Now run the comparison -----------------------

    rv_names_all <- rbind(itis_taxonomy,aphia_taxonomy)%>%
                    select(PhyloNames,format_id)%>%
                    gather(key=PhyloNames,value="value",-format_id)%>%
                    group_by(format_id)%>%
                    summarise(lowest_id = name_extract(value))%>%
                    ungroup()%>%
                    left_join(.,rbind(itis_taxonomy,aphia_taxonomy))%>%
                    distinct(lowest_id,.keep_all = TRUE)

    rv_names <- master_tax_trawl%>%
                mutate(format_id = 1:n())%>%
                select(PhyloNames,format_id)%>%
                gather(key=PhyloNames,value="value",-format_id)%>%
                group_by(format_id)%>%
                summarise(lowest_id = name_extract(value))%>%
                ungroup()%>%
                left_join(.,master_tax_trawl%>%mutate(format_id = 1:n()))%>%
                distinct(lowest_id,.keep_all = TRUE)
    
    #check for synonyms to make sure we have the right name
    
    
    eDNA_names <- master_tax_edna %>%
                  mutate(format_id = 1:n())%>%
                  select(PhyloNames,format_id)%>%
                  gather(key=PhyloNames,value="value",-format_id)%>%
                  group_by(format_id)%>%
                  summarise(lowest_id = name_extract(value))%>%
                  ungroup()%>%
                  left_join(.,master_tax_edna%>%mutate(format_id = 1:n()))%>%
                  distinct(lowest_id,.keep_all = TRUE)
                  
    #mutate(new_id = ifelse(is.na(taxID),taxID_extract(x=lowest_id,db=db),taxID)) #fill in the missing IDs
                  
    #Fill in missing IDs this doesn't work in a mutate step for some reason 
    eDNA_names$newID = eDNA_names$taxID
    for(i in 1:nrow(eDNA_names)){
        
        if(is.na(eDNA_names[i,"newID"])){
            
            eDNA_names[i,"newID"] <- taxID_extract(x=eDNA_names[i,]%>%pull(lowest_id),
                                                   db=eDNA_names[i,]%>%pull(db))    
            
        }
        
    }
    
    #check the itis names to make sure they are valid. I noticed that raja radiata was there which is not the right name. 
    # this might be a hold over from step2 where the first row was selected from classification(), in some cases that is not the valid name
    #
    
    eDNA_names <- eDNA_names%>%
                  mutate(valid = lowest_id,
                         valid_id = newID,
                         valid = ifelse(valid == "Hippoglossina oblonga","Paralichthys oblongus",valid),#for some reason this doesn't work with get_tsn in synonym_check()
                         valid = ifelse(valid == "Phyllodoce groenlandica","Anaitides groenlandica",valid),
                         valid = ifelse(valid == "Pseudopotamilla reniformis","Potamilla reniformis",valid))
    
    message("Note that Cyanea returns multiple itis responses so you have to review the list and enter the right response, usually '32'")
    message("This takes a while but it is in a for loop so you can intervene if an species id error crops up.")
    
    for(i in 215:nrow(eDNA_names)){ #this check takes quite a while (~ 40 minutes)
        
        message(paste0("working on ",i," of ",eDNA_names%>%filter(db=="itis")%>%nrow()))
        
        if(eDNA_names[i,"db"]=="itis"){
            eDNA_names[i,"valid"] = synonym_check(eDNA_names[i,]%>%pull(valid),db="itis",return_tsn = FALSE)
            eDNA_names[i,"valid_id"] = synonym_check(eDNA_names[i,]%>%pull(valid),db="itis",return_tsn = TRUE)
        }
        
    } #not all are itis so the last will print really quick
    

    
mismatch <- eDNA_names%>%
            filter(lowest_id != valid)%>%
            data.frame()

mismatch_tax <- pblapply(mismatch%>%pull(valid_id),Classify)%>%
                do.call("rbind",.)%>%
                dplyr::select(-input,-db)%>%
                rename(valid_id = id)%>%
                mutate(valid_id = as.numeric(valid_id))%>%
                left_join(.,mismatch%>%select(-PhyloNames))%>%
                dplyr::select(eDNA_names%>%names)

eDNA_cleaned <- rbind(eDNA_names%>%filter(lowest_id == valid),
                      mismatch_tax)%>%
                mutate(lowest_id = valid,
                       taxID = valid_id)%>%
                dplyr::select(-newID,-valid,-valid_id)%>%
                data.frame()
    
save.image("June3.RData")   
save(eDNA_cleaned,file="data/eDNA_cleaned.RData")
                
