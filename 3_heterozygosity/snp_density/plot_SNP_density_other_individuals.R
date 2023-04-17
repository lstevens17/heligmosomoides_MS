library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

# read in files and transform SNP density
HB2 <- read.table("HB2_vs_nxHelBake1.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

HB3 <- read.table("HB3_vs_nxHelBake1.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

HP2 <- read.table("HP2_vs_ngHelPoly1.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

# write transformed dataframes to file
write.table(HB2, "HB2_vs_nxHelBake1.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)
write.table(HB3, "HB3_vs_nxHelBake1.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)
write.table(HP2, "HP2_vs_ngHelPoly1.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)

# before plotting, remove any tiny bins and remove non-chromosomal scaffolds
HB2 <- HB2 %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>%
  filter(nonrepspan>500)

HB3 <- HB3 %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>% 
  filter(nonrepspan>500)

HP2 <- HP2 %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>% 
  filter(nonrepspan>500)

# plot
p1 <- ggplot(data=HP2, aes(x=start/1e6, y=transformed_snpdens)) +
  geom_point(alpha=0.3, size=0.4) +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1", size=1) + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.05)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle(expression(paste(italic("H. polygyrus"), " ngHelPoly2"))) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


# plot
p2 <- ggplot(data=HB2, aes(x=start/1e6, y=transformed_snpdens)) +
  geom_point(alpha=0.3, size=0.4) +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1") + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.05)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle(expression(paste(italic("H. bakeri"), " nxHelBake2"))) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


# plot
p3 <- ggplot(data=HB3, aes(x=start/1e6, y=transformed_snpdens)) +
  geom_point(alpha=0.3, size=0.4) +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1") + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.05)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle(expression(paste(italic("H. bakeri"), " nxHelBake2"))) + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


p <-p1 / p2 / p3 + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 18, face="bold"))

# save plots
ggsave("nxHelBake2_nxHelBake3_ngHelPoly2_SNP_density_10kb.png", plot = p, width=10, height=9, units="in")
ggsave("nxHelBake2_nxHelBake3_ngHelPoly2_SNP_density_10kb.pdf", plot = p, width=10, height=9, units="in")
