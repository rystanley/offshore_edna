setwd("~/GRDI_202002/AbundanceComparison")
library(ggplot2)
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

verts_count_16S <- read.csv("merged_16S_verts_count_sorted_formatted2.csv", header=TRUE)

lme_verts_count_16S <- lmer(log_eDNA_16S ~ log_verts_count + (1|Taxa), data = verts_count_16S)
summary(lme_verts_count_16S)

lme_verts_count_16S_pd <- ggpredict(lme_verts_count_16S,terms="log_verts_count",full.data=FALSE)
lme_verts_count_16S_pd

p_verts_count_16S <- ggplot(lme_verts_count_16S_pd, aes(x, predicted)) + 
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) + main_theme +
  labs(x=paste("trawling | lg(verts_count + 1)"),
       y=paste("16S | lg(reads + 1)")) + annotate(geom="text", x=0.42, y=6, label="(a) 16S") +
  geom_point(data=verts_count_16S, aes(x=log_verts_count, y=log_eDNA_16S), color="blue")

p_verts_count_16S

cor(verts_count_16S$log_eDNA_16S, verts_count_16S$log_verts_count, method="spearman")
cor(verts_count_16S$log_eDNA_16S, verts_count_16S$log_verts_count, method="pearson")

