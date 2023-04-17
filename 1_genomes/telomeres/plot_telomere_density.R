library(ggplot2)
library(tidyverse)
library(patchwork)

HB1 <- read.table("nxHelBake1.1.primary.fa_telomere_density_1kb.bed", col.names=c("chr", "start", "end", "count", "span", "window", "dens"))

p1 <- ggplot(data=HB1, aes(x=start/1e6, y=count)) + 
  geom_line() + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + xlab("Position (Mb)") + ylab("Counts of TTAGGC\nper 1kb window") + 
  ggtitle(expression(paste(italic("H. bakeri"), " nxHelBake1.1")))

HP1 <- read.table("ngHelPoly1.1.primary.fa_telomere_density_1kb.bed", col.names=c("chr", "start", "end", "count", "span", "window", "dens"))

p2 <- ggplot(data=HP1, aes(x=start/1e6, y=count)) + 
  geom_line() + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + xlab("Position (Mb)") + ylab("Counts of TTAGGC\nper 1kb window") + 
  ggtitle(expression(paste(italic("H. polygyrus"), " ngHelPoly1.1")))


p <- p1/p2 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18, face="bold"))

ggsave("telomere_distributions.png", plot = p, width=14, height=7, units="in")
ggsave("telomere_distributions.pdf", plot = p, width=14, height=7, units="in")
