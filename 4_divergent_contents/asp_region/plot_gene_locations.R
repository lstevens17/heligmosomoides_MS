library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

nxHelBake1_primary <- read.table("nxHelBake1_primary_coordinates.tsv", header=T, sep="\t")
nxHelBake2_primary <- read.table("nxHelBake2_primary_coordinates.tsv", header=T, sep="\t")
ngHelPoly1_primary <- read.table("ngHelPoly1_primary_coordinates.tsv", header=T, sep="\t")
ngHelPoly1_alternate <- read.table("ngHelPoly1_alternate_coordinates.tsv", header=T, sep="\t")
ngHelPoly2_primary <- read.table("ngHelPoly2_primary_coordinates.tsv", header=T, sep="\t")

p <- ggplot() + 
  geom_hline(yintercept=26) + 
  geom_hline(yintercept=20) +
  geom_hline(yintercept=14) + 
  geom_hline(yintercept=8) + 
  geom_hline(yintercept=2) + 
  geom_rect(data=nxHelBake1_primary, aes(xmin=start+40e3, xmax=end+40e3, ymin=28, ymax=24, fill=type), color="black") + 
  geom_rect(data=nxHelBake2_primary, aes(xmin=start-50e3, xmax=end-50e3, ymin=22, ymax=18, fill=type), color="black") +
  geom_rect(data=ngHelPoly1_primary, aes(xmin=start+15e3, xmax=end+15e3, ymin=16, ymax=12, fill=type), color="black") +
  geom_rect(data=ngHelPoly1_alternate, aes(xmin=start+0e3, xmax=end+0e3, ymin=10, ymax=6, fill=type), color="black") +
  geom_rect(data=ngHelPoly2_primary, aes(xmin=start+25e3, xmax=end+25e3, ymin=4, ymax=0, fill=type), color="black") +
  theme_void()

p

ggsave("gene_locations.pdf", plot = p, width=12, height=12, units="in")