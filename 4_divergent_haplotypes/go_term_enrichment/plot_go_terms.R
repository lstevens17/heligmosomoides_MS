library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

csv <- read.csv("tableS5_goterm.csv", header=T)

p <- ggplot(csv, aes(x=reorder(description, -p_weight), y=-log2(p_weight), size=div_count/expected)) + 
  geom_point(fill="#0092c1", pch=21) + coord_flip() + 
  geom_hline(yintercept=-log2(0.05), linetype=2) + 
  ylab("-log2(p-value)") + xlab("GO term") + 
  theme_bw() + scale_size(breaks = c(2, 4, 6, 8), name="Fold enrichment") + 
  theme(legend.position = c(0.85, 0.2), legend.box.background = element_rect(color="black", size=0.5))

p
ggsave("goterm_enrichment.png", plot = p, width=10, height=5, units="in")
  