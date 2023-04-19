#PC02 plot- match up to Figure 1

##load libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(scales)
library(tidyr)

#load PCoA ordination data outputs from CHIME II
trawling <- read.table("data/trawling_all_species_jaccard_PCoA.txt", header=T)%>%
            mutate(type = "Ecosystem Trawl")%>%
            select(PC1,PC2,bioclass,coarse_class,strata_class,type)

eDNA <- read.table("data/12S16SCO1_PCoA_results.txt", header=T)%>%
          mutate(type = "eDNA")%>%
          select(PC1,PC2,bioclass,coarse_class,strata_class,type)

plot_data <- rbind(trawling,eDNA)%>%
             mutate(bioclass = case_when(bioclass == "ESS_Banks" ~ "ESS: Banks",
                                         bioclass == "LaurentianChannel_ShelfBreak" ~ "Laurentian Channel/Shelf Break",
                                         bioclass == "WSS_Banks_InnerBoF" ~ "WSS/Inner BoF",
                                         bioclass == "WSS_OuterBoF" ~ "WSS: Banks/Inner BoF",
                                         TRUE ~ "ESS"))

#fill cols to match the default 6 fill colours for ggplot that were used in Figure 1
fill_cols <- hue_pal()(6)[c(1:3,5,6)] # the 4'th fill colour was for the 'slope' to which we don't have samples

#shapes to match the code by XH 
pco2_shapes<-c(19,15,17,9,7)


#assemble each plot -----------
trawl_plot <- ggplot(plot_data%>%filter(type=="Ecosystem Trawl"),aes(x=PC1,y=PC2,fill=bioclass,shape=bioclass))+
              geom_point(alpha=.9, size=4,col="black") +
              scale_fill_manual(values=fill_cols) + 
              scale_shape_manual(values=c(21:25)) +
              labs(x= "PCoA 1: 13.1%", y="PCoA 2: 10.4%",fill="",shape="",title = "A) Ecosystem trawl" ) + 
              theme(legend.position = c(0.82, 0.8),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    text=element_text(family="sans")) +
              theme_bw()

eDNA_plot <- ggplot(plot_data%>%filter(type=="eDNA"),aes(x=PC1,y=PC2*-1,fill=bioclass,shape=bioclass))+
            geom_point(alpha=.9, size=4,col="black") +
            scale_fill_manual(values=fill_cols) + 
            scale_shape_manual(values=c(21:25)) +
            labs(x= "PCoA 1: 9.2%", y="PCoA 2: 8.3%",fill="",shape="",title = "B) eDNA" ) + 
            theme(legend.position = c(0.82, 0.8),
                  axis.line=element_line(color="black"),
                  axis.ticks=element_line(color="black"),
                  axis.text=element_text(colour="black", size=10),
                  text=element_text(family="sans")) +
            theme_bw()+
            scale_y_continuous(position="right")

#combine the plots together. 
combo_plot <- trawl_plot + eDNA_plot + plot_layout(ncol=2,guides = "collect") & theme(legend.position = "bottom")

ggsave("output/PCo2_plot.png",combo_plot,width=8,height=6,units="in",dpi=300)

