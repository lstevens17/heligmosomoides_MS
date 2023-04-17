library(ggplot2)
library(tidyverse)
library(dplyr)
library(patchwork)

# read in SNP files and transform SNP density
nxHelBake1.alternate <- read.table("nxHelBake1.alternate.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

nxHelBake2.primary <- read.table("nxHelBake2.primary.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

nxHelBake2.alternate <- read.table("nxHelBake2.alternate.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

nxHelBake3.primary <- read.table("nxHelBake3.primary.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

nxHelBake3.alternate <- read.table("nxHelBake3.alternate.variant_repeat_density_10kb.bed", col.names=c("chr", "start", "stop", "snp_count", "snp_span", "window_size", "snp_dens", "rep_count", "rep_span", "rep_dens")) %>%
  mutate(nonrepspan=window_size-rep_span) %>%
  mutate(transformed_snpdens=snp_span/nonrepspan)

# read in divergent haplotype files
nxHelBake1.alternate.div <- read.table("nxHelBake1.alternate_divergent.bed", col.names=c("chr", "start", "stop")) %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X")

nxHelBake2.primary.div <- read.table("nxHelBake2.primary_divergent.bed", col.names=c("chr", "start", "stop")) %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X")

nxHelBake2.alternate.div <- read.table("nxHelBake2.alternate_divergent.bed", col.names=c("chr", "start", "stop")) %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X")

nxHelBake3.primary.div <- read.table("nxHelBake3.primary_divergent.bed", col.names=c("chr", "start", "stop")) %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X")

nxHelBake3.alternate.div <- read.table("nxHelBake3.alternate_divergent.bed", col.names=c("chr", "start", "stop")) %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X")


# write transformed dataframes to file
write.table(nxHelBake1.alternate, "nxHelBake1.alternate.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)
write.table(nxHelBake2.primary, "nxHelBake2.primary.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)
write.table(nxHelBake2.alternate, "nxHelBake2.alternate.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)
write.table(nxHelBake3.primary, "nxHelBake3.primary.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)
write.table(nxHelBake3.alternate, "nxHelBake3.alternate.variant_repeat_density_10kb.transformed.bed", sep="\t", row.names=FALSE)

# before plotting, remove any tiny bins and remove non-chromosomal scaffolds
nxHelBake1.alternate <- nxHelBake1.alternate %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>%
  filter(nonrepspan>500)

nxHelBake2.primary <- nxHelBake2.primary %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>% 
  filter(nonrepspan>500)

nxHelBake2.alternate <- nxHelBake2.alternate %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>% 
  filter(nonrepspan>500)

nxHelBake3.primary <- nxHelBake3.primary %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>% 
  filter(nonrepspan>500)

nxHelBake3.alternate <- nxHelBake3.alternate %>%
  filter(chr == "I" | chr == "II" | chr == "III" |  chr == "IV" |  chr == "V" | chr == "X") %>% 
  filter(nonrepspan>500)

# plot
p1 <- ggplot() +
  geom_point(data=nxHelBake1.alternate, aes(x=start/1e6, y=transformed_snpdens), alpha=0.3, size=0.4) +
  geom_rect(data=nxHelBake1.alternate.div, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0.1, ymax=0.11), colour="red") +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1", size=1) + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.11)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle("nxHelBake1 alternate") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p2 <- ggplot() +
  geom_point(data=nxHelBake2.primary, aes(x=start/1e6, y=transformed_snpdens), alpha=0.3, size=0.4) +
  geom_rect(data=nxHelBake2.primary.div, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0.1, ymax=0.11), colour="red") +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1", size=1) + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.11)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle("nxHelBake2 primary") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p3 <- ggplot() +
  geom_point(data=nxHelBake2.alternate, aes(x=start/1e6, y=transformed_snpdens), alpha=0.3, size=0.4) +
  geom_rect(data=nxHelBake2.alternate.div, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0.1, ymax=0.11), colour="red") +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1", size=1) + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.11)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle("nxHelBake2 alternate") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p4 <- ggplot() +
  geom_point(data=nxHelBake3.primary, aes(x=start/1e6, y=transformed_snpdens), alpha=0.3, size=0.4) +
  geom_rect(data=nxHelBake3.primary.div, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0.1, ymax=0.11), colour="red") +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1", size=1) + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.11)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle("nxHelBake3 primary") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p5 <- ggplot() +
  geom_point(data=nxHelBake3.alternate, aes(x=start/1e6, y=transformed_snpdens), alpha=0.3, size=0.4) +
  geom_rect(data=nxHelBake3.alternate.div, aes(xmin=start/1e6, xmax=stop/1e6, ymin=0.1, ymax=0.11), colour="red") +
  geom_smooth(method="loess", span=0.05, se=F, color="#0092c1", size=1) + 
  facet_wrap(~chr, scales="free_x", nrow=1) +
  xlab("Position (Mb)") + ylab("SNP density") + theme_bw() + 
  coord_cartesian(ylim=c(0,0.11)) + scale_y_continuous(expand=c(0,0)) + 
  ggtitle("nxHelBake3 alternate") + 
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

p <-p1 / p2 / p3 + p4 + p5 + plot_annotation(tag_levels = 'A')  &
  theme(plot.tag = element_text(size = 18, face="bold"))

# save plots
ggsave("SNP_density_divergent_haplotype_locations.png", plot = p, width=10, height=12, units="in")
ggsave("SNP_density_divergent_haplotype_locations.pdf", plot = p, width=10, height=12, units="in")
