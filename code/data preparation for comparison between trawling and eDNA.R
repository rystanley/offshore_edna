setwd("~/Documents/GitHub/offshore_edna/data/comparison")

#read data
data_16S <- read.csv("16S_final3_2021Oct13.csv", header=TRUE)
data_verts_wgt <- read.csv("verts_wgt_wide_2021Oct13.csv", header=TRUE)
data_verts_count <- read.csv("verts_count_wide_2021Oct13.csv", header=TRUE)

merged_16S_verts_wgt <- merge(data_16S,data_verts_wgt,by="Taxa")
View(merged_16S_verts_wgt)
dim(merged_16S_verts_wgt)
write.csv(merged_16S_verts_wgt,'merged_16S_verts_wgt.csv')

merged_16S_verts_count <-merge(data_16S, data_verts_count, by="Taxa")
write.csv(merged_16S_verts_count,"merged_16S_verts_count.csv")
