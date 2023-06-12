library(ggplot2)
library(tidyverse)

dS_table <- read.table("dS_by_group.tsv", col.names = c("dS", "group")) %>% filter(dS <= 0.8) 

p <- ggplot(dS_table, aes(x=group, y=dS, fill=group)) + 
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width=0.2, color="black") +
#  geom_point(position="jitter", alpha=0.2) + 
  theme_classic() + 
  scale_fill_manual(values=c("#5B84B1FF", "#FC766AFF")) +
  scale_x_discrete(labels=c("Single-copy orthologues\nin non-divergent haplotypes\n(N=9773)","Single-copy orthologues\nin shared haplotypes\n(N=189)")) + 
  coord_cartesian(ylim=c(0, 0.2)) + 
  theme(legend.position = "none") + 
  xlab("Haplotype status") + ylab(expression(d[S]))

ggsave("dS_boxplot.png", plot = p, width=4, height=4, units="in")
ggsave("dS_boxplot.pdf", plot = p, width=4, height=4, units="in")

wilcox.test(dS_table$dS ~ dS_table$group, alternative = "two.sided")  

