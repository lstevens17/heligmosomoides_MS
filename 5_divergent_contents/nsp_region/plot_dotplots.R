library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

nxHelBake2 <- read.table("nxHelBake2_ragtag_vs_nxHelBake1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p1 <- nxHelBake2 %>% filter(LEN1 > 500) %>%
  ggplot(.) +
  geom_segment(aes(x=(S1+580000)/1e6, xend=(E1+580000)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) +
  theme_bw() +
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Position in nxHelBake2 primary (kb)") +
  coord_cartesian(xlim=c(580000/1e6,1150000/1e6), ylim=c(132635/1e3,712589/1e3), expand=c(0,0)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

ngHelPoly1 <- read.table("ngHelPoly1_vs_nxHelBake1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p2 <- ngHelPoly1 %>% filter(LEN1 > 500) %>%
  ggplot(.) +
  geom_segment(aes(x=(S1+580000)/1e6, xend=(E1+580000)/1e6, y=(S2+1066773)/1e6, yend=(E2+1066773)/1e6, col=IDY), linewidth=4) +
  theme_bw() +
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ylab("Position in ngHelPoly1 primary (Mb)") +
  coord_cartesian(xlim=c(580000/1e6,1150000/1e6), expand=c(0,0)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

p <- p1 / p2 

ggsave("nxHelBake2_ngHelPoly1_dotplots.png", plot = p, width=8, height=8, units="in")
ggsave("nxHelBake2_ngHelPoly1_dotplots.pdf", plot = p, width=8, height=8, units="in")
