library(tidyverse)

# resd in codeML table
bmalayi <- read.table("bmalayi_codeml_dS.updated.tsv", col.names=c("gene", "dS")) %>% filter(dS <= 0.8) %>% mutate(sp = "bmalayi_vs_bpahangi")
cbriggsae <- read.table("cbriggsae_codeml_dS.updated.tsv", col.names=c("gene", "dS")) %>% filter(dS <= 0.8) %>% mutate(sp = "cbriggsae_vs_cnigoni")
cremanei <- read.table("cremanei_codeml_dS.updated.tsv", col.names=c("gene", "dS")) %>% filter(dS <= 0.8) %>% mutate(sp = "cremanei_vs_clatens")
hbakeri <- read.table("nxHelBake1_codeml_dS.updated.tsv", col.names=c("gene", "dS")) %>% filter(dS <= 0.8) %>% mutate(sp = "hbakeri_vs_hpolygyrus")
ovolvulus <- read.table("ovolvulus_codeml_dS.updated.tsv", col.names=c("gene", "dS")) %>% filter(dS <= 0.8) %>% mutate(sp = "ovolvulus_vs_oochengi")

# merge tables and write out a TSV
merged <- rbind(bmalayi, cbriggsae, cremanei, hbakeri, ovolvulus)
write.table(merged, file="comparative_dS.tsv", sep="\t", row.names=F, quote=F)

# make a table out of mean dS for each species
mean_dS <- matrix(c("Brugia malayi\nvs\nBrugia pahangi", "Caenorhabditis briggsae\nvs\nCaenorhabditis nigoni", "Caenorhabditis remanei\nvs\nCaenorhabditis latens", "Heligmosomoides bakeri\nvs\nHeligmosomoides polygyrus", "Onchocerca volvulus\nvs\nOnchocerca ochengi", 
         mean(bmalayi$dS), mean(cbriggsae$dS), mean(cremanei$dS), mean(hbakeri$dS), mean(ovolvulus$dS)), ncol=2)
colnames(mean_dS) <- c("species","dS")
mean_dS <- data.frame(mean_dS, stringsAsFactors = FALSE)
cols.num <- c("dS")
mean_dS[cols.num] <- sapply(mean_dS[cols.num],as.numeric)

# plot
p <- ggplot() + 
  geom_col(data=mean_dS, aes(x=reorder(species, dS), y=dS*100), fill="#798cbd", width = 0.8) + 
  theme_bw() + 
  xlab("Species pair") + ylab("Synonymous site divergence (%)")

# save plots
ggsave("comparative_dS.png", plot = p, width=8, height=6, units="in")
ggsave("comparative_dS.pdf", plot = p, width=8, height=6, units="in")
