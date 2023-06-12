library(ggplot2)
library(tidyverse)
library(scales)
library(patchwork)

# p1
table1 <- read.table("hap_ptg002319l_1_1_vs_nxHelBake1_V.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p1 <- table1 %>% filter(LEN1 > 0) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+71784000)/1e6, xend=(E1+71784000)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 primary chromosome V (Mb)") + ylab("Position in nxHelBake1 alternate\nhap_ptg002319l_1_1 (kb)") + 
  geom_vline(xintercept=71815752/1e6, linetype=2) + geom_vline(xintercept=71846478/1e6, linetype=2) + 
  coord_cartesian(ylim=c(125000/1e3,232000/1e3), xlim=c(71784000/1e6,71878000/1e6))  +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# p2
table2 <- read.table("ngHelPoly1_V_vs_nxHelBake1_V.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p2 <- table2 %>% filter(LEN1 > 0) %>%
  ggplot(.) + 
  geom_segment(aes(x=(S1+71784000)/1e6, xend=(E1+71784000)/1e6, y=(S2+68770407)/1e6, yend=(E2+68770407)/1e6, col=IDY), linewidth=4) + 
  theme_bw() + 
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  xlab("Position in nxHelBake1 primary chromosome V (Mb)") + ylab("Position in ngHelPoly1 primary\nchromosome V (Mb)") + 
  geom_vline(xintercept=71815752/1e6, linetype=2) + geom_vline(xintercept=71846478/1e6, linetype=2) + 
  coord_cartesian(ylim=c((26000+68770407)/1e6,(130000+68770407)/1e6), xlim=c(71784000/1e6,71878000/1e6)) 

p <- p1 / p2  

ggsave("ngHelPoly1_primary_nxHelBake1_alternate_vs_nxHelBake1_primary_dotplots.png", plot = p, width=8, height=6, units="in")
ggsave("ngHelPoly1_primary_nxHelBake1_alternate_vs_nxHelBake1_primary_dotplots.pdf", plot = p, width=8, height=6, units="in")