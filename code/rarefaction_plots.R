#load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(patchwork)

#load data and format it for plotting
taxdat_12s <- read.csv("data/12S_final5_taxa_table_richness-rarefaction-10000.csv")%>%
          gather(key = "var",value = "count",depth.1_iter.1:depth.10000_iter.10)%>%
          separate(var,into=c("var1","var2"),sep="_")%>%
          mutate(depth=as.numeric(gsub("depth.","",var1)),
                 iteration=gsub("iter.","",var2),
                 marker="12s")%>%
          select(sample.id,depth,iteration,count,marker)

taxdat_16s <- read.csv("data/16S_final5_taxa_table_richness_rarefaction-10000.csv")%>%
              gather(key = "var",value = "count",depth.1_iter.1:depth.10000_iter.10)%>%
              separate(var,into=c("var1","var2"),sep="_")%>%
              mutate(depth=as.numeric(gsub("depth.","",var1)),
                     iteration=gsub("iter.","",var2),
                     marker="16s")%>%
              select(sample.id,depth,iteration,count,marker)

taxdat_CO1 <- read.csv("data/CO1_final2_taxa_table_richness-rarefaction-10000.csv")%>%
              gather(key = "var",value = "count",depth.1_iter.1:depth.10000_iter.10)%>%
              separate(var,into=c("var1","var2"),sep="_")%>%
              mutate(depth=as.numeric(gsub("depth.","",var1)),
                     iteration=gsub("iter.","",var2),
                     marker="CO1")%>%
              select(sample.id,depth,iteration,count,marker)



taxdat <- rbind(taxdat_12s,taxdat_16s,taxdat_CO1)

taxdat_processed <- taxdat%>%
                    group_by(marker,sample.id,depth)%>%
                    summarise(mn=mean(count,na.rm=T),
                              sd=sd(count,na.rm=T),
                              se=sd/sqrt(10))%>%
                    ungroup()%>%
                    data.frame()

groups <- taxdat_processed%>%
          group_by(marker,sample.id)%>%
          summarise(max=max(mn,na.rm=T))%>%
          ungroup()%>%
          data.frame()%>%
          arrange(marker,max)

taxdat_processed2 <- taxdat_processed%>%
                     left_join(.,groups)

group_num <- taxdat_processed2%>%
             group_by(marker,max)%>%
             summarise(sample_n = length(unique(sample.id)))%>%
             ungroup()%>%
             data.frame()

group_depth <- taxdat_processed2%>%
               filter(!is.na(mn))%>%
               group_by(marker,max)%>%
               summarise(max_depth = max(depth,na.rm=T))%>%
               ungroup()%>%
               data.frame()
                     
taxdat_processed3 <- taxdat_processed2%>%
                     left_join(.,group_num)%>%
                     mutate(alpha = 1/sample_n) #this will get the alpha right

count_df <- taxdat_processed3%>%
            group_by(marker,max,depth,mn)%>%
            summarise(count=length(unique(sample.id)))%>%
            ungroup()%>%
            left_join(.,group_num)%>%
            mutate(alpha = 1/sample_n*count)%>%
            data.frame()

#now construct the plotting elements for each marker, with slight adjustments to the plots for assembly using patchwork

#12s primer
  plot_12s <- ggplot()+
    geom_point(data=count_df%>%filter(depth>1,marker=="12s"),
               aes(x=depth,y=mn,size=count),pch=21,col="black",show.legend=FALSE)+
    geom_line(data=taxdat_processed3%>%filter(!is.na(mn),marker=="12s"),
              aes(x=depth,y=mn,group=sample.id,col=factor(max),alpha=alpha),lwd=1.5,show.legend=FALSE)+
    geom_point(data=count_df%>%filter(depth>1,marker=="12s"),
               aes(x=depth,y=mn,fill=factor(max),alpha=alpha,size=count),pch=21,col="black",show.legend=FALSE)+
    geom_point(aes(x=0,y=1),col="black",fill="white",pch=21,size=3)+
    labs(x="Sequencing depth",y="",size="# obs",title="A) 12s")+
    theme_bw()+
    scale_y_continuous(breaks=seq(2,15,2))+
    scale_color_viridis(discrete=T,option="C")+
    scale_size_continuous(breaks=c(1,5,9,14))+
    scale_fill_viridis(discrete=T,option="C")+
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
#16s primer
  plot_16s <- ggplot()+
    geom_point(data=count_df%>%filter(depth>1,marker=="16s"),
               aes(x=depth,y=mn,size=count),pch=21,col="black")+
    geom_line(data=taxdat_processed3%>%filter(!is.na(mn),marker=="16s"),
              aes(x=depth,y=mn,group=sample.id,col=factor(max),alpha=alpha),lwd=1.5,show.legend=FALSE)+
    geom_point(data=count_df%>%filter(depth>1,marker=="16s"),
               aes(x=depth,y=mn,fill=factor(max),alpha=alpha,size=count),pch=21,col="black",show.legend=FALSE)+
    geom_point(aes(x=0,y=1),col="black",fill="white",pch=21,size=3)+
    labs(x="Sequencing depth",y="",size="# obs",title="B) 16s")+
    theme_bw()+
    scale_y_continuous(breaks=seq(2,15,2))+
    scale_color_viridis(discrete=T,option="C")+
    scale_size_continuous(breaks=c(1,5,9,14))+
    scale_fill_viridis(discrete=T,option="C")+
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

#CO1 primer
  plot_co1 <- ggplot()+
    geom_point(data=count_df%>%filter(depth>1,marker=="CO1"),
               aes(x=depth,y=mn,size=count),pch=21,col="black",show.legend=FALSE)+
    geom_line(data=taxdat_processed3%>%filter(!is.na(mn),marker=="CO1"),
              aes(x=depth,y=mn,group=sample.id,col=factor(max),alpha=alpha),lwd=1.5,show.legend=FALSE)+
    geom_point(data=count_df%>%filter(depth>1,marker=="CO1"),
               aes(x=depth,y=mn,fill=factor(max),alpha=alpha,size=count),pch=21,col="black",show.legend=FALSE)+
    geom_point(aes(x=0,y=1),col="black",fill="white",pch=21,size=3)+
    labs(x="Sequencing depth",y="",size="# obs",title="C) CO1")+
    theme_bw()+
    scale_y_continuous(breaks=seq(2,25,2))+
    scale_color_viridis(discrete=T,option="C")+
    scale_size_continuous(breaks=c(1,5,9,14),)+
    scale_fill_viridis(discrete=T,option="C")+
    theme(legend.position = "none",
          axis.title.y = element_blank()) 
  
  p_lab <- ggplot() + 
            annotate(geom = "text", x = 1, y = 1, label = "Mean taxon richness", angle = 90) +
            coord_cartesian(clip = "off")+
            theme_void()

#combined into one plot 
  combo_plot = p_lab + {(plot_12s+plot_16s+plot_co1+guide_area()) + plot_layout(nrow=2,guides="collect") & theme(legend.position="right")} + plot_layout(ncol=2,widths=c(0.1,1))

#save plot
  ggsave("output/taxonrichness_seqdepth.png",combo_plot,width=7,height=6,units="in",dpi=300)
