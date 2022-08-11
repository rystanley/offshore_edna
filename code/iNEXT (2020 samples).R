#load libaries -----
    library(iNEXT)
    library(ggplot2)
    library(gridExtra)
    library(patchwork)
    library(dplyr)
    library(viridis)

#Richness analysis

  #read in data
    east_fish <-read.table("data/fish_iNEXT.txt",header=T)
    east_invert <-read.table("data/data_iNEXT_invertebrates.txt",header=T)
   
  #set number of samples 
    east.num <- 108
  
  #Richness analysis
    vert_next <- iNEXT(east_fish, q=0, datatype="incidence_freq", endpoint=east.num, knots = 40, se = TRUE, conf = 0.95,
                     nboot = 1000)
    
    invert_next <- iNEXT(east_invert, q=0, datatype="incidence_freq", endpoint=east.num, knots = 40, se = TRUE, conf = 0.95,
                             nboot = 1000)
    
   #Extract plotting data
    plotdata <- rbind(fortify.iNEXT(vert_next,type=1)%>%mutate(type="Vertebrates"),
                      fortify.iNEXT(vert_next,type=3)%>%mutate(type="Vertebrates"),
                      fortify.iNEXT(invert_next,type=1)%>%mutate(type="Invertebrates"),
                      fortify.iNEXT(invert_next,type=3)%>%mutate(type="Invertebrates"))%>%
      mutate(var=gsub("fish_","",site),
             var=gsub("_invertebrates","",var),
             var=ifelse(var == "trawling","Trawling",var),
             var=factor(var,levels=c("Trawling","12S","16S","CO1","eDNA")),
             denominator = ifelse(plottype==1,"Number of sampling sites","Sample coverage"))%>%
      arrange(var,type,x)%>%
      dplyr::select(site,var,type,method,x,y,y.lwr,y.upr,denominator)
    
    #assemble data for the asymptotic diversity estimates
    assym_div <- rbind(vert_next$AsyEst,invert_next$AsyEst)%>%
                 mutate(type=ifelse(grepl("fish",Site,),"Fish","Invertebrates"),
                       var=gsub("fish_","",Site),
                       var=gsub("_invertebrates","",var),
                       var=ifelse(var == "trawling","Trawling",var),
                       var=factor(var,levels=c("Trawling","12S","16S","CO1","eDNA")))
    
    
    #Rarefaction plots
    plotcols <- RColorBrewer::brewer.pal(5,"Dark2")
    plotcols[3] <- "grey30" #increase the contrast of 16S
    
    p1 <- ggplot()+
      geom_ribbon(data=plotdata,aes(x=x,ymin=y.lwr,ymax=y.upr,fill=var),alpha=0.5)+
      geom_line(data=plotdata,aes(x=x,y=y,col=var),lty=2,lwd=1.25)+
      geom_line(data=filter(plotdata,method!="extrapolated"),aes(x=x,y=y,col=var),lwd=1.25)+
      geom_point(data=filter(plotdata,method=="observed"),aes(x=x,y=y,fill=var),shape=21,size=4)+ #change in type to make it clearer.
      theme_bw()+
      facet_grid(type~denominator,scales="free")+
      theme(strip.background = element_rect(, colour = "black", fill = "white"),
            strip.text.x = element_text(colour = "black",size=14), 
            strip.text.y = element_text(colour = "black",size=14),
            axis.text = element_text(colour = "black",size=12),
            axis.title = element_text(colour = "black",size=12),
            legend.position="bottom")+
      scale_y_continuous(expand=c(0,0.02))+
      scale_fill_manual(values=plotcols)+ # to increase contrast
      scale_colour_manual(values=plotcols)+
      labs(x="",y="Species Richness",fill="",col="");p1
    
    ggsave("output/Figure3.png",p1,width=9,height=9,units="in",dpi=600)
    
    #break the plots out
    p2 <- ggplot()+
      geom_ribbon(data=plotdata%>%filter(denominator == "Number of sampling sites"),aes(x=x,ymin=y.lwr,ymax=y.upr,fill=var),alpha=0.5)+
      geom_line(data=plotdata%>%filter(denominator == "Number of sampling sites"),aes(x=x,y=y,col=var),lty=2,lwd=1.25)+
      geom_line(data=filter(plotdata,method!="extrapolated",denominator == "Number of sampling sites"),aes(x=x,y=y,col=var),lwd=1.25)+
      geom_point(data=filter(plotdata,method=="observed",denominator == "Number of sampling sites"),aes(x=x,y=y,fill=var),shape=21,size=4)+ #change in type to make it clearer.
      theme_bw()+
      facet_grid(type~.,scales="free")+
      theme(strip.background = element_blank(),
            strip.text = element_blank(),
            axis.text = element_text(colour = "black",size=12),
            axis.title = element_text(colour = "black",size=12),
            legend.position="bottom")+
      scale_y_continuous(expand=c(0,0.02))+
      scale_fill_manual(values=plotcols)+ # to increase contrast
      scale_colour_manual(values=plotcols)+
      labs(x="Number of sampling sites",y="Species Richness",fill="",col="")
    
    #just the sample size rarefaction plot
    plotdata <- plotdata%>%
                mutate(type2 = ifelse(type=="Vertebrates","Fish",type))
    
    p2_solo <- ggplot()+
              geom_ribbon(data=plotdata%>%filter(denominator == "Number of sampling sites"),aes(x=x,ymin=y.lwr,ymax=y.upr,fill=var),alpha=0.5)+
              geom_line(data=plotdata%>%filter(denominator == "Number of sampling sites"),aes(x=x,y=y,col=var),lty=2,lwd=1.25)+
              geom_line(data=filter(plotdata,method!="extrapolated",denominator == "Number of sampling sites"),aes(x=x,y=y,col=var),lwd=1.25)+
              geom_point(data=filter(plotdata,method=="observed",denominator == "Number of sampling sites"),aes(x=x,y=y,fill=var),shape=21,size=4)+ #change in type to make it clearer.
              theme_bw()+
              facet_grid(type2~.,scales="free")+
              theme(strip.background = element_rect(, colour = "black", fill = "white"),
                    strip.text.x = element_text(colour = "black",size=14), 
                    strip.text.y = element_text(colour = "black",size=14),
                    axis.text = element_text(colour = "black",size=12),
                    axis.title = element_text(colour = "black",size=12),
                    legend.position="bottom")+
              scale_y_continuous(expand=c(0,0.02))+
              scale_fill_manual(values=plotcols)+ # to increase contrast
              scale_colour_manual(values=plotcols)+
              labs(x="Number of sampling sites",y="Species Richness",fill="",col="")
    
    ggsave("output/richness_sites.png",p2_solo,width=9,height=9,units="in",dpi=600)
    
    p3 <- ggplot()+
      geom_ribbon(data=plotdata%>%filter(denominator == "Sample coverage"),aes(x=x,ymin=y.lwr,ymax=y.upr,fill=var),alpha=0.5)+
      geom_line(data=plotdata%>%filter(denominator == "Sample coverage"),aes(x=x,y=y,col=var),lty=2,lwd=1.25)+
      geom_line(data=filter(plotdata,method!="extrapolated",denominator == "Sample coverage"),aes(x=x,y=y,col=var),lwd=1.25)+
      geom_point(data=filter(plotdata,method=="observed",denominator == "Sample coverage"),aes(x=x,y=y,fill=var),shape=21,size=4)+ #change in type to make it clearer.
      theme_bw()+
      facet_grid(type~.,scales="free")+
      theme(strip.background = element_rect(, colour = "black", fill = "white"),
            strip.text.x = element_text(colour = "black",size=14), 
            strip.text.y = element_text(colour = "black",size=14),
            axis.text.y = element_blank(),
            axis.text.x = element_text(colour = "black",size=12),
            axis.title.y = element_blank(),
            axis.title.x = element_text(colour = "black",size=12),
            plot.margin = unit(c(5.5, 5.5, 5.5, 0), "pt"),# have to adjust this to make it look like a facet wrap
            #axis.ticks.y = element_blank(), #if you don't want the ticks
            legend.position="bottom")+
      scale_y_continuous(expand=c(0,0.02))+
      scale_fill_manual(values=plotcols)+ # to increase contrast
      scale_colour_manual(values=plotcols)+
      labs(x="Sample coverage",y="Species Richness",fill="",col="")
    
    output <- p2 + p3 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    #save the plot
    ggsave("output/RichnessPlots.png",output,width=9,height=9,units="in",dpi=600)
    
    
    