library(ggplot2)
library(tidyverse)
library(patchwork)
library(forcats)

compare <- read.table("nxHelBake1_ngHelPoly1_nigon_position_table.tsv", col.names=c("buscoID", "nigon", "gene1", "chr1", "start1", "stop1", "num1", "gene2", "chr2", "start2", "stop2", "num2"))

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "unassigned" = "grey")

p1 <- compare %>% filter(chr1 == "X") %>% filter(chr2 == "X") %>%
  ggplot(., aes(x=start1/1e6, y=start2/1e6, color=nigon)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(fct_rev(chr1)~chr2, scales="free") + theme_bw() + 
  scale_colour_manual(values=cols) + xlab(expression(paste("Position in ", italic("H. bakeri"), " (Mb)"))) + 
  ylab(expression(paste("Position in ", italic("H. polygyrus"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), 
        panel.grid.minor = element_blank(), 
        legend.position = c(0.9, 0.15),
        legend.background = element_rect(linetype="solid", color="black", linewidth=0.2),
        strip.background = element_blank(),
        strip.text = element_blank()) + 
  xlim(0, 92046022/1e6) + ylim(0,97362115/1e6) + 
  geom_vline(xintercept=40078085/1e6, linetype=2, alpha=0.5) + 
  geom_vline(xintercept=51748017/1e6, linetype=2, alpha=0.5) + 
  geom_hline(yintercept=56814194/1e6, linetype=2, alpha=0.5) + 
  geom_hline(yintercept=43767944/1e6, linetype=2, alpha=0.5) 

nxHelBake1 <- read.table("nxHelBake1.1.primary.fa.out.gff_100kb.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens"))
ngHelPoly1 <- read.table("ngHelPoly1.1.primary.fa.out.gff_100kb.bed", col.names=c("chr", "start", "stop", "count", "span", "window", "dens"))

p2 <- nxHelBake1 %>% filter(chr == "X") %>%
  ggplot(aes(x=start/1e6, y=dens*100)) + 
  geom_point(alpha=0.25) + 
  geom_smooth(span=0.5) +
  theme_bw() + 
  ylab("Repeat density (%)") + 
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt")) + ylim(0,100) + 
  geom_vline(xintercept=40078085/1e6, linetype=2, alpha=0.5) + 
  geom_vline(xintercept=51748017/1e6, linetype=2, alpha=0.5) 

p3 <- ngHelPoly1 %>% filter(chr == "X") %>%
  ggplot(aes(x=start/1e6, y=dens*100)) + 
  geom_point(alpha=0.25) + 
  geom_smooth(span=0.5) +
  theme_bw() + 
  ylab("Repeat density (%)") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0,0,0,0), "pt")) + 
  geom_vline(xintercept=56814194/1e6, linetype=2, alpha=0.5) + 
  geom_vline(xintercept=43767944/1e6, linetype=2, alpha=0.5) +
  coord_flip() + ylim(0,100)

design <- "AAAAA##
           BBBBBCC
           BBBBBCC
           BBBBBCC"

p <- p2 + p1 + p3 + plot_layout(design = design)

ggsave("x_inversion_repeat_density.png", plot = p, width=8, height=7, units="in")
ggsave("x_inversion_repeat_density.pdf", plot = p, width=8, height=7, units="in")
