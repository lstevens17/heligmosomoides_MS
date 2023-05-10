library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

nxHelBake2 <- read.table("nxHelBake2_ragtag_vs_nxHelBake1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p1 <- nxHelBake2 %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+8923671)/1e6, xend=(E1+8923671)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), size=2) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 chromosome I (Mb)") + ylab("Position in nxHelBake2 primary (kb)") + 
  coord_cartesian(xlim=c(8923671/1e6,9199833/1e6), ylim=c(54100/1e3,429344/1e3))

nxHelBake3 <- read.table("nxHelBake3_ragtag_vs_nxHelBake1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p2 <- nxHelBake3 %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+8923671)/1e6, xend=(E1+8923671)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), size=2) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 chromosome I (Mb)") + ylab("Position in nxHelBake3 primary (kb)") + 
  coord_cartesian(xlim=c(8923671/1e6,9199833/1e6), ylim=c(44872/1e3,411789/1e3))

p <- p1 / p2

ggsave("nxHelBake2_nxHelBake3_ragtag_dotplots.png", plot = p, width=8, height=8, units="in")
ggsave("nxHelBake2_nxHelBake3_ragtag_dotplots.pdf", plot = p, width=8, height=8, units="in")
