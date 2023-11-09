library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

compare <- read.table("nxHelBake1_ngHelPoly1_nigon_position_table.tsv", col.names=c("buscoID", "nigon", "gene1", "chr1", "start1", "stop1", "num1", "gene2", "chr2", "start2", "stop2", "num2"))
minimap2_bed <- read.table("nxHelBake1_vs_ngHelPoly1_minimap2.bed", col.names=c("chr", "start", "stop", "length", "div")) %>% filter(length > 5e3) %>%
  filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X")
dS_locations <- read.table(file = "nxHelBake1_codeml_dS.with_location.updated.tsv", col.names = c("id", "chr", "start", "stop", "dS")) %>% filter(chr == "I" | chr == "II" | chr == "III" | chr == "IV" | chr == "V" | chr == "X") %>% filter(dS < 0.8)
write.table(dS_locations, file="nxHelBake1_codeml_dS.with_location.updated.filtered.tsv", row.names=F, sep="\t", quote=F)

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b", "unassigned" = "grey")

p1 <- compare %>%
  ggplot(., aes(x=start1/1e6, y=start2/1e6, color=nigon)) +
  geom_point(alpha=0.5, size=2) + facet_grid(fct_rev(chr1)~chr2, scales="free") + theme_bw() +
  scale_colour_manual(values=cols, name="Nigon") + xlab(expression(paste("Position in ", italic("H. bakeri"), " (Mb)"))) +
  ylab(expression(paste("Position in ", italic("H. polygyrus"), " (Mb)")))  +
  theme(panel.spacing=unit(0, "lines"), panel.grid.minor = element_blank(), 
        legend.position=c(0.91,0.29), 
        legend.background = element_rect(size=0.5, linetype="solid", colour="black"))

p2 <- ggplot(data=minimap2_bed) +
  geom_rect(aes(xmin=start/1e6, xmax=stop/1e6, ymax=(100-div)+0.5, ymin=(100-div)-0.5), fill="black", color="black", linewidth=0.15) +
  facet_grid(~chr, scales="free_x") +
  theme_bw() +
  ylab("Nucleotide identity (%)") +
  coord_cartesian(ylim=c(85,100)) + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

p3 <- ggplot(data=dS_locations, aes(x=((start+stop)/2)/1e6, y=dS*100)) + 
  geom_point(alpha=0.2) + 
  geom_smooth(alpha=0.2, colour="#49abaa") + 
  geom_hline(yintercept=mean(dS_locations$dS*100), linetype=2) +
  coord_cartesian(ylim=c(0,17)) + 
  facet_grid(~chr, scales="free_x") + 
  theme_bw() + 
  xlab("Position (Mb)") + ylab("Synonymous site divergence (%)") + 
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))

p <- (p1 + p1) / p2 / p3 +
  plot_layout(heights = c(2, 1, 1)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18, face="bold")) 

ggsave("Figure2_draft.pdf", plot = p, width=12, height=11, units="in")
