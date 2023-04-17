library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

bed <- read.table("nxHelBake1_vs_ngHelPoly1_minimap2.bed", col.names=c("chr", "start", "stop", "length", "div")) %>% filter(length > 5e3) %>% 
  filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X")

p1 <- ggplot(data=bed) + 
  geom_rect(aes(xmin=start/1e6, xmax=stop/1e6, ymax=(100-div)+0.5, ymin=(100-div)-0.5), fill="black", color="black", linewidth=0.1) + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Nucleotide identity (%)") + 
  coord_cartesian(ylim=c(85,100)) 

protein_identity_location <- read.table("nxHelBake1_vs_ngHelPoly1_protein_identity_location.tsv", col.names=c("seq", "chr", "start", "stop", "id")) %>% 
  filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X") %>% mutate(pos = (start+stop)/2) %>% filter(id >= 85)

p2 <- ggplot(data=protein_identity_location, aes(x=pos/1e6, y=id)) + 
  geom_point(alpha=0.2) + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Amino acid identity (%)") + 
  geom_smooth(span=1, colour="#49abaa") + 
  coord_cartesian(ylim=c(95,100)) 

p <- p1 / p2

ggsave("nxHelBake1_vs_ngHelPoly1_nucleotide_protein_divergence.png", plot = p, width=10, height=5, units="in")
ggsave("nxHelBake1_vs_ngHelPoly1_nucleotide_protein_divergence.pdf", plot = p, width=10, height=5, units="in")