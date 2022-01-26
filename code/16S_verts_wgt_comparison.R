library(ggplot2)
setwd("~/Documents/GitHub/offshore_edna")

require(lme4)
require(DHARMa)
require(sjPlot)
require(sjmisc)
require(sjlabelled)
require(ggplot2)
require(ggeffects)
require(RColorBrewer)

main_theme <- theme(panel.background=element_blank(),
                    panel.grid=element_blank(),
                    axis.line=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    text=element_text(family="sans"))

############ detected by two methods ###########
verts_wgt_16S <- read.table("data/comparison/verts_wgt_16S_formatted_for_abundance_comparison.txt", header=TRUE)

lme_verts_wgt_16S <- lmer(log_eDNA_16S ~ log_verts_wgt + (1|Taxa), data = verts_wgt_16S)
summary(lme_verts_wgt_16S)

lme_verts_wgt_16S_pd <- ggpredict(lme_verts_wgt_16S,terms="log_verts_wgt",full.data=FALSE)
lme_verts_wgt_16S_pd

p_verts_wgt_16S <- ggplot(lme_verts_wgt_16S_pd, aes(x, predicted)) + 
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) + main_theme +
  labs(x=paste("trawling | lg(verts_wgt)"),
       y=paste("16S | lg(reads + 1)")) + annotate(geom="text", x=0.42, y=6, label="16S reads number vs trawling_wgt") +
  geom_point(data=verts_wgt_16S, aes(x=log_verts_wgt, y=log_eDNA_16S), color="blue")

p_verts_wgt_16S

cor(verts_wgt_16S$log_eDNA_16S, verts_wgt_16S$log_verts_wgt, method="spearman")
cor(verts_wgt_16S$log_eDNA_16S, verts_wgt_16S$log_verts_wgt, method="pearson")
