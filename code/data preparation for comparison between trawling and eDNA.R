setwd("~/Documents/GitHub/offshore_edna/data/comparison")

############ 16S reads number
data_16S <- read.csv("16S_final3_2021Oct13.csv", header=TRUE)
data_verts_wgt <- read.csv("verts_wgt_wide_2021Oct13.csv", header=TRUE)
data_verts_count <- read.csv("verts_count_wide_2021Oct13.csv", header=TRUE)

merged_16S_verts_wgt <- merge(data_16S,data_verts_wgt,by="Taxa")
View(merged_16S_verts_wgt)
dim(merged_16S_verts_wgt)
write.csv(merged_16S_verts_wgt,'merged_16S_verts_wgt.csv')

merged_16S_verts_count <-merge(data_16S, data_verts_count, by="Taxa")
write.csv(merged_16S_verts_count,"merged_16S_verts_count.csv")

############ 16S relative abundance
data_16S_rel <- read.csv("2020_16S_fish_rel_table_2022Jan23.csv", header=TRUE)

merged_16S_rel_verts_wgt <- merge(data_16S_rel,data_verts_wgt,by="Taxa")
View(merged_16S_rel_verts_wgt)
dim(merged_16S_rel_verts_wgt)
write.csv(merged_16S_rel_verts_wgt,'merged_16S_rel_verts_wgt.csv')

merged_16S_rel_verts_count <-merge(data_16S_rel, data_verts_count, by="Taxa")
write.csv(merged_16S_rel_verts_count,"merged_16S_rel_verts_count.csv")
