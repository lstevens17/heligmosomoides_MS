library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

left = 8923671
right = 9199833

t0 <- read.table("compare_to_nxHelBake1/nxHelBake2_asp_contigs.ragtag_vs_nxHelBake1_I_8923671-9199833.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p0 <- t0 %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+left)/1e6, xend=(E1+left)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 chromosome I (Mb)") + ylab("Position in\nnxHelBake2 primary (kb)") + 
  coord_cartesian(xlim=c(left/1e6,right/1e6), ylim=c(54100/1e3,429344/1e3)) +
  geom_vline(xintercept=9027639/1e6, linetype=2) + geom_vline(xintercept=9040340/1e6, linetype=2) + 
  geom_vline(xintercept=9113978/1e6, linetype=2) + geom_vline(xintercept=9119297/1e6, linetype=2) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))

t1 <- read.table("compare_to_nxHelBake1/ngHelPoly1_I_8624310-9007409_vs_nxHelBake1_I_8923671-9199833.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p1 <- t1 %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+left)/1e6, xend=(E1+left)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 chromosome I (Mb)") + ylab("Position in\nngHelPoly1 primary (kb)") + 
  coord_cartesian(xlim=c(left/1e6,right/1e6)) +
  geom_vline(xintercept=9027639/1e6, linetype=2) + geom_vline(xintercept=9040340/1e6, linetype=2) + 
  geom_vline(xintercept=9113978/1e6, linetype=2) + geom_vline(xintercept=9119297/1e6, linetype=2) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))

t2 <- read.table("compare_to_nxHelBake1/ngHelPoly1_alternate_asp_contigs.ragtag_vs_nxHelBake1_I_8923671-9199833.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p2 <- t2 %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+left)/1e6, xend=(E1+left)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 chromosome I (Mb)") + ylab("Position in\nngHelPoly1 alternate (kb)") + 
  coord_cartesian(xlim=c(left/1e6,right/1e6)) +
  geom_vline(xintercept=9027639/1e6, linetype=2) + geom_vline(xintercept=9040340/1e6, linetype=2) + 
  geom_vline(xintercept=9113978/1e6, linetype=2) + geom_vline(xintercept=9119297/1e6, linetype=2) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin = unit(c(0,0,0,0), "cm"))

t3 <- read.table("compare_to_nxHelBake1/ngHelPoly2_primary_asp_contigs.ragtag_vs_nxHelBake1_I_8923671-9199833.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p3 <- t3 %>% filter(LEN1 > 500) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+left)/1e6, xend=(E1+left)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 chromosome I (Mb)") + ylab("Position in\nngHelPoly2 primary (kb)") + 
  coord_cartesian(xlim=c(left/1e6,right/1e6)) +
  geom_vline(xintercept=9027639/1e6, linetype=2) + geom_vline(xintercept=9040340/1e6, linetype=2) + 
  geom_vline(xintercept=9113978/1e6, linetype=2) + geom_vline(xintercept=9119297/1e6, linetype=2) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p <- p0 / p1 / p2 / p3 

ggsave("all_haplotypes_vs_nxHelBake1_dotplots.png", plot = p, width=8, height=12, units="in")
ggsave("all_haplotypes_vs_nxHelBake1_dotplots.pdf", plot = p, width=8, height=12, units="in")
