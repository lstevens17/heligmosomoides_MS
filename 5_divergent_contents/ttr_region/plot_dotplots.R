library(ggplot2)
library(tidyverse)
library(scales)
library(patchwork)

# p1
table1 <- read.table("nxHelBake1.alternate_ragtag_vs_nxHelBake1.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p1 <- table1 %>% filter(LEN1 > 0) %>%
  ggplot(.) +
  geom_segment(aes(x=(S1+20891115)/1e6, xend=(E1+20891115)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) +
  theme_bw() +
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Position in nxHelBake1 primary chromosome V (Mb)") + ylab("Position in nxHelBake1 alternate\nhap_ptg002319l_1_1 (kb)") +
  #geom_vline(xintercept=20925054/1e6, linetype=2) + geom_vline(xintercept=20925705/1e6, linetype=2) + # g651
  #geom_vline(xintercept=20929090/1e6, linetype=2) + geom_vline(xintercept=20932648/1e6, linetype=2) + # g652
  geom_vline(xintercept=20937142/1e6, linetype=2) + geom_vline(xintercept=20938184/1e6, linetype=2) + # g653
  geom_vline(xintercept=20948347/1e6, linetype=2) + geom_vline(xintercept=20949981/1e6, linetype=2) + # g654
  geom_vline(xintercept=20966215/1e6, linetype=2) + geom_vline(xintercept=20974311/1e6, linetype=2) + # g655
  #geom_vline(xintercept=20991583/1e6, linetype=2) + geom_vline(xintercept=20994085/1e6, linetype=2) + # g657
  #coord_cartesian(xlim=c(20900000/1e6,21050000/1e6), ylim=c(40000/1e3, 190000/1e3)) + 
  coord_cartesian(xlim=c(20900000/1e6,21008000/1e6), ylim=c(40000/1e3, 159000/1e3)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# p2
table2 <- read.table("ngHelPoly1_II_vs_nxHelBake1_II.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p2 <- table2 %>% filter(LEN1 > 0) %>%
  ggplot(.) +
  geom_segment(aes(x=(S1+20891115)/1e6, xend=(E1+20891115)/1e6, y=(S2+20228318)/1e6, yend=(E2+20228318)/1e6, col=IDY), linewidth=4) +
  theme_bw() +
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Position in nxHelBake1 primary chromosome V (Mb)") + ylab("Position in\nngHelPoly1 primary (Mb)") +
  geom_vline(xintercept=20937142/1e6, linetype=2) + geom_vline(xintercept=20938184/1e6, linetype=2) + # g653
  geom_vline(xintercept=20948347/1e6, linetype=2) + geom_vline(xintercept=20949981/1e6, linetype=2) + # g654
  geom_vline(xintercept=20966215/1e6, linetype=2) + geom_vline(xintercept=20974311/1e6, linetype=2) + # g655
  coord_cartesian(xlim=c(20900000/1e6,21008000/1e6), ylim=c((20000+20228318)/1e6, (140000+20228318)/1e6)) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# p3
table3 <- read.table("ngHelPoly2_primary_vs_nxHelBake1_II.filtered.coords", col.names = c("S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "CHROM1", "CHROM2"))

p3 <- table3 %>% filter(LEN1 > 0) %>%
  ggplot(.) +
  geom_segment(aes(x=(S1+20891115)/1e6, xend=(E1+20891115)/1e6, y=S2/1e3, yend=E2/1e3, col=IDY), linewidth=4) +
  theme_bw() +
  scale_colour_gradient(low = "#eb8034", high = "#1b365f", name = "Nucleotide\nidentity (%)", limits = c(85,100), oob=squish) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Position in nxHelBake1 primary chromosome V (Mb)") + ylab("Position in\nngHelPoly2 primary (kb)") +
  geom_vline(xintercept=20937142/1e6, linetype=2) + geom_vline(xintercept=20938184/1e6, linetype=2) + # g653
  geom_vline(xintercept=20948347/1e6, linetype=2) + geom_vline(xintercept=20949981/1e6, linetype=2) + # g654
  geom_vline(xintercept=20966215/1e6, linetype=2) + geom_vline(xintercept=20974311/1e6, linetype=2) + # g655
  coord_cartesian(xlim=c(20900000/1e6,21008000/1e6), ylim=c(157980/1e3, 280000/1e3))

p <- p1 / p2 / p3

ggsave("transthyretin_dotplots.png", plot = p, width=8, height=8, units="in")
ggsave("transthyretin_dotplots.pdf", plot = p, width=8, height=8, units="in")
