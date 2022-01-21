library(dplyr)
setwd("~/Documents/GitHub/offshore_edna")

####################                    16S Verts Count                     ####################
GRDI_data1 <- read.table("data/comparison/merged_16S_verts_count_for_formatting.txt", header=TRUE)
View(GRDI_data1)
dim(GRDI_data1)
#30, 109

# change column names
names(GRDI_data1) <- c("Taxa","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S",
                       "verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count","eDNA_16S","verts_count")

# formatting
Col1 <- select(GRDI_data1, 1)
View(Col1)

Sam1 <- select(GRDI_data1, 2:3); Sample1 <-cbind(Col1, Sam1)
View(Sample1)

Sam2<-select(GRDI_data1, 4:5); Sam3<-select(GRDI_data1, 6:7); Sam4<-select(GRDI_data1, 8:9); Sam5 <- select(GRDI_data1, 10:11); Sam6 <- select(GRDI_data1, 12:13)
Sam7<-select(GRDI_data1, 14:15); Sam8<-select(GRDI_data1, 16:17); Sam9<-select(GRDI_data1, 18:19); Sam10 <- select(GRDI_data1, 20:21); Sam11 <- select(GRDI_data1, 22:23)
Sam12<-select(GRDI_data1, 24:25); Sam13<-select(GRDI_data1, 26:27); Sam14<-select(GRDI_data1, 28:29); Sam15 <- select(GRDI_data1, 30:31); Sam16 <- select(GRDI_data1, 32:33)
Sam17<-select(GRDI_data1, 34:35); Sam18<-select(GRDI_data1, 36:37); Sam19<-select(GRDI_data1, 38:39); Sam20 <- select(GRDI_data1, 40:41); Sam21 <- select(GRDI_data1, 42:43)
Sam22<-select(GRDI_data1, 44:45); Sam23<-select(GRDI_data1, 46:47); Sam24<-select(GRDI_data1, 48:49); Sam25 <- select(GRDI_data1, 50:51); Sam26 <- select(GRDI_data1, 52:53)
Sam27<-select(GRDI_data1, 54:55); Sam28<-select(GRDI_data1, 56:57); Sam29<-select(GRDI_data1, 58:59); Sam30 <- select(GRDI_data1, 60:61); Sam31 <- select(GRDI_data1, 62:63)
Sam32<-select(GRDI_data1, 64:65); Sam33<-select(GRDI_data1, 66:67); Sam34<-select(GRDI_data1, 68:69); Sam35 <- select(GRDI_data1, 70:71); Sam36 <- select(GRDI_data1, 72:73)
Sam37<-select(GRDI_data1, 74:75); Sam38<-select(GRDI_data1, 76:77); Sam39<-select(GRDI_data1, 78:79); Sam40 <- select(GRDI_data1, 80:81); Sam41 <- select(GRDI_data1, 82:83)
Sam42<-select(GRDI_data1, 84:85); Sam43<-select(GRDI_data1, 86:87); Sam44<-select(GRDI_data1, 88:89); Sam45 <- select(GRDI_data1, 90:91); Sam46 <- select(GRDI_data1, 92:93)
Sam47<-select(GRDI_data1, 94:95); Sam48<-select(GRDI_data1, 96:97); Sam49<-select(GRDI_data1, 98:99); Sam50 <- select(GRDI_data1, 100:101); Sam51 <- select(GRDI_data1, 102:103)
Sam52<-select(GRDI_data1, 104:105); Sam53<-select(GRDI_data1, 106:107); Sam54<-select(GRDI_data1, 108:109)

Sample2 <-cbind(Col1, Sam2); Sample3 <-cbind(Col1, Sam3); Sample4 <-cbind(Col1, Sam4); Sample5 <-cbind(Col1, Sam5); Sample6 <-cbind(Col1, Sam6)
Sample7 <-cbind(Col1, Sam7); Sample8 <-cbind(Col1, Sam8); Sample9 <-cbind(Col1, Sam9); Sample10 <-cbind(Col1, Sam10); Sample11 <-cbind(Col1, Sam11)
Sample12 <-cbind(Col1, Sam12); Sample13 <-cbind(Col1, Sam13); Sample14 <-cbind(Col1, Sam14); Sample15 <-cbind(Col1, Sam15); Sample16 <-cbind(Col1, Sam16)
Sample17 <-cbind(Col1, Sam17); Sample18 <-cbind(Col1, Sam18); Sample19 <-cbind(Col1, Sam19); Sample20 <-cbind(Col1, Sam20); Sample21 <-cbind(Col1, Sam21)
Sample22 <-cbind(Col1, Sam22); Sample23 <-cbind(Col1, Sam23); Sample24 <-cbind(Col1, Sam24); Sample25 <-cbind(Col1, Sam25); Sample26 <-cbind(Col1, Sam26)
Sample27 <-cbind(Col1, Sam27); Sample28 <-cbind(Col1, Sam28); Sample29 <-cbind(Col1, Sam29); Sample30 <-cbind(Col1, Sam30); Sample31 <-cbind(Col1, Sam31)
Sample32 <-cbind(Col1, Sam32); Sample33 <-cbind(Col1, Sam33); Sample34 <-cbind(Col1, Sam34); Sample35 <-cbind(Col1, Sam35); Sample36 <-cbind(Col1, Sam36)
Sample37 <-cbind(Col1, Sam37); Sample38 <-cbind(Col1, Sam38); Sample39 <-cbind(Col1, Sam39); Sample40 <-cbind(Col1, Sam40); Sample41 <-cbind(Col1, Sam41)
Sample42 <-cbind(Col1, Sam42); Sample43 <-cbind(Col1, Sam43); Sample44 <-cbind(Col1, Sam44); Sample45 <-cbind(Col1, Sam45); Sample46 <-cbind(Col1, Sam46)
Sample47 <-cbind(Col1, Sam47); Sample48 <-cbind(Col1, Sam48); Sample49 <-cbind(Col1, Sam49); Sample50 <-cbind(Col1, Sam50); Sample51 <-cbind(Col1, Sam51)
Sample52 <-cbind(Col1, Sam52); Sample53 <-cbind(Col1, Sam53); Sample54 <-cbind(Col1, Sam54)

verts_count_16S_formatted<-rbind(Sample1,Sample2,Sample3,Sample4,Sample5,Sample6,Sample7,Sample8,Sample9,Sample10,Sample11,Sample12,
         Sample13,Sample14,Sample15,Sample16,Sample17,Sample18,Sample19,Sample20,Sample21,Sample22,Sample23,
         Sample24,Sample25,Sample26,Sample27,Sample28,Sample29,Sample30,Sample31,Sample32,Sample33,Sample34,
         Sample35,Sample36,Sample37,Sample38,Sample39,Sample40,Sample41,Sample42,Sample43,Sample44,Sample45,
         Sample46,Sample47,Sample48,Sample49,Sample50,Sample51,Sample52,Sample53,Sample54)

write.csv(verts_count_16S_formatted,'data/comparison/verts_count_16S_formatted.csv')
