library(ggplot2)
library(tidyverse)
library(patchwork)

benchmark <- read.table("benchmark_values.combined.reformatted.tsv", col.names=c("flank_identity", "flank_length", "length", "snp_dens", "bins", "cov", "spanLS", "intersect", "spanDL", "sens", "spec"))

p1 <- ggplot(data=benchmark, aes(x=sens, y=spec, color=factor(flank_identity))) + 
  geom_point(size=2, alpha=0.75) + 
  theme_bw() + 
  scale_colour_discrete(name="Flanking alignment\nidentity (%)") + 
  theme(legend.position = c(0.13, 0.33), legend.background = element_rect(linewidth=0.2, linetype="solid", color="black")) + 
  xlab("Sensitivity (%)") + ylab("Specificity (%)") + 
  coord_cartesian(xlim=c(40,90), ylim=c(40,90)) +
  geom_abline(intercept = 0, slope = 1, linetype=2)  + 
  geom_point(aes(x=81.3482, y=79.0856), colour="black", pch=23, size=5)

p2 <- ggplot(data=benchmark, aes(x=sens, y=spec, color=factor(flank_length))) + 
  geom_point(size=2, alpha=0.75) + 
  theme_bw() + 
  scale_colour_discrete(name="Flanking alignment\nlength (kb)") + 
  theme(legend.position = c(0.13, 0.27), legend.background = element_rect(linewidth=0.2, linetype="solid", color="black"))  + 
  xlab("Sensitivity (%)") + ylab("Specificity (%)") + 
  coord_cartesian(xlim=c(40,90), ylim=c(40,90)) +
  geom_abline(intercept = 0, slope = 1, linetype=2)  + 
  geom_point(aes(x=81.3482, y=79.0856), colour="black", pch=23, size=5)


p3 <- ggplot(data=benchmark, aes(x=sens, y=spec, color=factor(length))) + 
  geom_point(size=2, alpha=0.75) + 
  theme_bw() + 
  scale_colour_discrete(name="Length (kb)") + 
  theme(legend.position = c(0.1, 0.18), legend.background = element_rect(linewidth=0.2, linetype="solid", color="black"))  + 
  xlab("Sensitivity (%)") + ylab("Specificity (%)") + 
  coord_cartesian(xlim=c(40,90), ylim=c(40,90)) +
  geom_abline(intercept = 0, slope = 1, linetype=2)  + 
  geom_point(aes(x=81.3482, y=79.0856), colour="black", pch=23, size=5)

p4 <- ggplot(data=benchmark, aes(x=sens, y=spec, color=factor(snp_dens))) + 
  geom_point(size=2, alpha=0.75) + 
  theme_bw() + 
  scale_colour_discrete(name="SNP density") + 
  theme(legend.position = c(0.1, 0.25), legend.background = element_rect(linewidth=0.2, linetype="solid", color="black")) + 
  xlab("Sensitivity (%)") + ylab("Specificity (%)") + 
  coord_cartesian(xlim=c(40,90), ylim=c(40,90)) +
  geom_abline(intercept = 0, slope = 1, linetype=2)  + 
  geom_point(aes(x=81.3482, y=79.0856), colour="black", pch=23, size=5)
 
p <- p1 + p2 + p3 + p4 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18, face="bold"))

ggsave("all_benchmarking.png", plot = p, width = 40, height = 24, units = "cm")
ggsave("all_benchmarking.pdf", plot = p, width = 40, height = 24, units = "cm")
