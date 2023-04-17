library(ggplot2)
library(tidyverse)
library(ggrepel)

HB1 <- read.table("nxHelBake1.primary.no_mt.contigs.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "nxHelBake1")
HB2 <- read.table("nxHelBake2.primary.no_mt.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "nxHelBake2")
HB3 <- read.table("nxHelBake3.primary.no_mt.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "nxHelBake3")
HP1 <- read.table("ngHelPoly1.primary.no_mt.contigs.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "ngHelPoly1")
HP2 <- read.table("ngHelPoly2.primary.no_mt.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "ngHelPoly2")
PRJEB1203 <- read.table("heligmosomoides_polygyrus.PRJEB1203.WBPS17.genomic.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "PRJEB1203")
PRJEB15396 <- read.table("heligmosomoides_polygyrus.PRJEB15396.WBPS17.genomic.fa_cumulative_length.tsv", col.names=c("count", "total")) %>% mutate(dataset = "PRJEB15396")

df <- rbind(HB1, HB2, HB3, HP1, HP2, PRJEB15396, PRJEB1203)

p <- ggplot(data=df, aes(x=count, y=total, colour=dataset)) + 
  geom_line(size=1.5) +
  xlab("Number of contigs") + ylab("Span (Mb)") + 
  scale_y_continuous(labels=function(x)x/1e6) +
  xlim(0,10000) + 
  theme_bw() + labs(color="Assembly") + 
  scale_colour_manual(values=c("#d98900", "#f5c067", "#1a7574", "#49abaa", "#94ebea", "lightgrey", "darkgrey"))

ggsave("cumulative_length_curves.pdf", plot = p, width=10, height=6, units="in")
ggsave("cumulative_length_curves.pdf", plot = p, width=10, height=6, units="in")