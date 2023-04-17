library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

bakeri_polygyrus <- read.table("nxHelBake1_vs_ngHelPoly1_protein_identity_location.tsv", col.names=c("seq", "chr", "start", "stop", "id")) %>% 
  filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X") %>% mutate(pos = (start+stop)/2) %>% filter(id >= 85)
briggsae_nigoni <- read.table("cbriggsae_vs_cnigoni_protein_identity_location.tsv", col.names=c("seq", "chr", "start", "stop", "id")) %>% 
  filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X")  %>% mutate(pos = (start+stop)/2) %>% filter(id >= 85)
remanei_latens <- read.table("cremanei_vs_clatens_protein_identity_location.tsv", col.names=c("seq", "chr", "start", "stop", "id")) %>% 
  filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X")  %>% mutate(pos = (start+stop)/2) %>% filter(id >= 85)

p1 <- ggplot(data=bakeri_polygyrus, aes(x=pos/1e6, y=id)) + 
  geom_point(alpha=0.2) + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Amino acid identity (%)") + 
  geom_smooth(span=1, colour="#49abaa") + 
  ggtitle(label=expression(paste(italic("H. bakeri"), " vs ", italic("H. polygyrus"), " (98.8%)"))) + 
  coord_cartesian(ylim=c(95,100)) 

p2 <- ggplot(data=briggsae_nigoni, aes(x=pos/1e6, y=id)) + 
  geom_point(alpha=0.2) + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Amino acid identity (%)") + 
  geom_smooth(span=1, colour="#49abaa") + 
  ggtitle(label=expression(paste(italic("Caenorhabditis briggsae"), " vs ", italic("Caenorhabditis nigoni"), " (97.9%)"))) + 
  coord_cartesian(ylim=c(95,100)) 

p3 <- ggplot(data=remanei_latens, aes(x=pos/1e6, y=id)) + 
  geom_point(alpha=0.2) + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Amino acid identity (%)") + 
  geom_smooth(span=1, colour="#49abaa") + 
  ggtitle(label=expression(paste(italic("Caenorhabditis remanei"), " vs ", italic("Caenorhabditis latens"), " (98.1%)"))) + 
  coord_cartesian(ylim=c(95,100)) 

p <- p1 / p2 / p3 + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18, face="bold"))

ggsave("compare_protein_divergence.png", plot = p, width=12, height=10, units="in")
ggsave("compare_protein_divergence.pdf", plot = p, width=12, height=10, units="in")





