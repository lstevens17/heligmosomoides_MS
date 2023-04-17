library(ggplot2)
library(tidyverse)
library(patchwork)
library(forcats)

compare <- read.table("nxHelBake1_ngHelPoly1_nigon_position_table.tsv", col.names=c("buscoID", "nigon", "gene1", "chr1", "start1", "stop1", "num1", "gene2", "chr2", "start2", "stop2", "num2"))

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "unassigned" = "grey")

p <- compare %>% 
  ggplot(., aes(x=start1/1e6, y=start2/1e6, color=nigon)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(fct_rev(chr1)~chr2, scales="free") + theme_bw() + 
  scale_colour_manual(values=cols) + xlab(expression(paste("Position in ", italic("H. bakeri"), " (Mb)"))) + 
  ylab(expression(paste("Position in ", italic("H. polygyrus"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank())

ggsave("nxHelBake1_ngHelPoly1_Nigon_oxford.legend.pdf", plot = p, width = 20, height = 18, units = "cm")
ggsave("nxHelBake1_ngHelPoly1_Nigon_oxford.legend.png", plot = p, width = 20, height = 18, units = "cm")

x <- compare %>% filter(chr1 == "X") %>% filter(chr2 == "X") %>%
  ggplot(., aes(x=start1/1e6, y=start2/1e6, color=nigon)) + 
  geom_point(alpha=0.5, size=2) + facet_grid(fct_rev(chr1)~chr2, scales="free") + theme_bw() + 
  scale_colour_manual(values=cols) + xlab(expression(paste("Position in ", italic("H. bakeri"), " (Mb)"))) + 
  ylab(expression(paste("Position in ", italic("H. polygyrus"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank())

ggsave("nxHelBake1_ngHelPoly1_Nigon_oxford.legend.X.pdf", plot = x, width = 20, height = 18, units = "cm")
ggsave("nxHelBake1_ngHelPoly1_Nigon_oxford.legend.X.png", plot = x, width = 20, height = 18, units = "cm")
