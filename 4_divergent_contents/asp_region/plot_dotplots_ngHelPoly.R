library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

ngHelPoly1_alt <- read.table("ngHelPoly1_alternate_ragtag_vs_ngHelPoly1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p1 <- ngHelPoly1_alt %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+8624310)/1e6, xend=(E1+8624310)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), size=2) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in ngHelPoly1 chromosome I (Mb)") + ylab("Position in ngHelPoly1 alternate (kb)") + 
  coord_cartesian(xlim=c(8624310/1e6,9007409/1e6), ylim=c(0/1e3,373797/1e3))

ngHelPoly2_pri <- read.table("ngHelPoly2_primary_ragtag_vs_ngHelPoly1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p2 <- ngHelPoly2_pri %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+8624310)/1e6, xend=(E1+8624310)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), size=2) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in ngHelPoly1 chromosome I (Mb)") + ylab("Position in ngHelPoly2 primary (kb)") + 
  coord_cartesian(xlim=c(8624310/1e6,9007409/1e6), ylim=c(23833/1e3,302680/1e3))

p <- p1 / p2

ggsave("ngHelPoly_ragtag_dotplots.png", plot = p, width=8, height=8, units="in")
ggsave("ngHelPoly_ragtag_dotplots.pdf", plot = p, width=8, height=8, units="in")
