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
verts_count_16S_rel <- read.table("data/comparison/verts_count_16S_rel_formatted_for_abundance_comparison.txt", header=TRUE)

lme_verts_count_16S_rel <- lmer(log_eDNA_16S ~ log_verts_count + (1|Taxa), data = verts_count_16S_rel)
summary(lme_verts_count_16S_rel)

lme_verts_count_16S_rel_pd <- ggpredict(lme_verts_count_16S_rel,terms="log_verts_count",full.data=FALSE)
lme_verts_count_16S_rel_pd

p_verts_count_16S_rel <- ggplot(lme_verts_count_16S_rel_pd, aes(x, predicted)) + 
  geom_line() + geom_point() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) + main_theme +
  labs(x=paste("trawling | lg(verts_count + 1)"),
       y=paste("16S | lg(relative abundance + 1)")) + annotate(geom="text", x=0.42, y=2.2, label="(a) 16S relative abundance vs vert_count") +
  geom_point(data=verts_count_16S_rel, aes(x=log_verts_count, y=log_eDNA_16S), color="blue")

p_verts_count_16S_rel


cor(verts_count_16S_rel$log_eDNA_16S, verts_count_16S_rel$log_verts_count, method="spearman")
cor(verts_count_16S_rel$log_eDNA_16S, verts_count_16S_rel$log_verts_count, method="pearson")

