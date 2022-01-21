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
############ detected by at least one method ###########
verts_count_16S <- read.table("data/comparison/verts_count_16S_formatted_for_abundance_comparison.txt", header=TRUE)

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

############ detected by two methods ###########
verts_count_16S_2m <- read.table("data/comparison/verts_count_16S_formatted_for_abundance_comparison_2methods.txt", header=TRUE)

lme_verts_count_16S_2m <- lmer(log_eDNA_16S ~ log_verts_count + (1|Taxa), data = verts_count_16S_2m)
summary(lme_verts_count_16S_2m)

lme_verts_count_16S_pd_2m <- ggpredict(lme_verts_count_16S_2m,terms="log_verts_count",full.data=FALSE)
lme_verts_count_16S_pd_2m

p_verts_count_16S_2m <- ggplot(lme_verts_count_16S_pd_2m, aes(x, predicted)) + 
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) + main_theme +
  labs(x=paste("trawling | lg(verts_count)"),
       y=paste("16S | lg(reads + 1)")) + annotate(geom="text", x=0.42, y=6, label="(a) 16S") +
  geom_point(data=verts_count_16S_2m, aes(x=log_verts_count, y=log_eDNA_16S), color="blue")

p_verts_count_16S_2m


