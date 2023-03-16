#load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)

#load data and format it for plotting
taxdat <- read.csv("data/12S_final5_taxa_table_richness-rarefaction-10000.csv")%>%
          gather(key = "var",value = "count",depth.1_iter.1:depth.10000_iter.10)%>%
          separate(var,into=c("var1","var2"),sep="_")%>%
          mutate(depth=as.numeric(gsub("depth.","",var1)),
                 iteration=gsub("iter.","",var2))%>%
          select(sample.id,depth,iteration,count)

taxdat_processed <- taxdat%>%
                    group_by(sample.id,depth)%>%
                    summarise(mn=mean(count,na.rm=T),
                              sd=sd(count,na.rm=T),
                              se=sd/sqrt(10))%>%
                    ungroup()%>%
                    data.frame()

groups <- taxdat_processed%>%
          group_by(sample.id)%>%
          summarise(max=max(mn,na.rm=T))%>%
          ungroup()%>%
          data.frame()%>%
          arrange(max)

taxdat_processed2 <- taxdat_processed%>%
                     left_join(.,groups)

group_num <- taxdat_processed2%>%
             group_by(max)%>%
             summarise(sample_n = length(unique(sample.id)))%>%
             ungroup()%>%
             data.frame()

group_depth <- taxdat_processed2%>%
               filter(!is.na(mn))%>%
               group_by(max)%>%
               summarise(max_depth = max(depth,na.rm=T))%>%
               ungroup()%>%
               data.frame()
                     
taxdat_processed3 <- taxdat_processed2%>%
                     left_join(.,group_num)%>%
                     mutate(alpha = 1/sample_n) #this will get the alpha right

count_df <- taxdat_processed3%>%
            group_by(max,depth,mn)%>%
            summarise(count=length(unique(sample.id)))%>%
            ungroup()%>%
            left_join(.,group_num)%>%
            mutate(alpha = 1/sample_n*count)%>%
            data.frame()

p1 <- ggplot()+
  geom_point(data=count_df%>%filter(depth>1),
             aes(x=depth,y=mn,size=count),pch=21,col="black")+
  geom_line(data=taxdat_processed3%>%filter(!is.na(mn)),
            aes(x=depth,y=mn,group=sample.id,col=factor(max),alpha=alpha),lwd=1.5,show.legend=FALSE)+
  geom_point(data=count_df%>%filter(depth>1),
             aes(x=depth,y=mn,fill=factor(max),alpha=alpha,size=count),pch=21,col="black",show.legend=FALSE)+
  geom_point(aes(x=0,y=1),col="black",fill="white",pch=21,size=3)+
  labs(x="Sequencing depth",y="Mean observed features",size="# obs")+
  theme_bw()+
  theme()+
  scale_y_continuous(breaks=seq(2,15,2))+
  scale_color_viridis(discrete=T,option="C")+
  scale_fill_viridis(discrete=T,option="C");p1 

ggsave("output/obsfeaturs_depth_plot.png",p1,width=6,height=6,units="in",dpi=300)
