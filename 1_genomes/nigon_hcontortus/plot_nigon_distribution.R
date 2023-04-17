library(ggplot2)
library(tidyverse)
library(patchwork)

nxHelBake1_nigon_location <- read_tsv("nxHelBake1_nigon_location.tsv")
window_size = 3000000

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

p1 <- filter(nxHelBake1_nigon_location) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(query_chr ~ ., scales = "free") +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), strip.text.y = element_text(angle = 0), legend.position = "None")

ngHelPoly1_nigon_location <- read_tsv("ngHelPoly1_nigon_location.tsv")
window_size = 3000000

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

p2 <- filter(ngHelPoly1_nigon_location) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(query_chr ~ ., scales = "free") +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), strip.text.y = element_text(angle = 0), legend.position = "None")

nxHelBake1_hcontortus_location <- read_tsv("nxHelBake1_hcontortus_location.tsv")
window_size = 3000000

p3 <- filter(nxHelBake1_hcontortus_location) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  #scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(query_chr ~ ., scales = "free") +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), strip.text.y = element_text(angle = 0), legend.position = "None")

ngHelPoly1_hcontortus_location <- read_tsv("ngHelPoly1_hcontortus_location.tsv")
window_size = 3000000

p4 <- filter(ngHelPoly1_hcontortus_location) %>%
  group_by(query_chr) %>%
  mutate(gene_count = n(), max_position = max(position)) %>%
  filter(gene_count > 30) %>%
  mutate(ints = as.numeric(as.character(cut(position,
                                            breaks = seq(0, max(position), window_size),
                                            labels = seq(window_size, max(position), window_size)))),
         ints = ifelse(is.na(ints), max(ints, na.rm = T) + window_size, ints)) %>%
  count(ints, assigned_chr, query_chr) %>%
  ungroup() %>%
  ggplot(aes(fill=assigned_chr, y=n, x=((ints-window_size)/1e6))) +
  #scale_fill_manual(values=cols, name = "Nigon") + 
  facet_grid(query_chr ~ ., scales = "free") +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  xlab("Position (Mb)") + ylab("BUSCO count (n)") +
  theme(panel.border = element_blank(), text = element_text(size=10), strip.text.y = element_text(angle = 0), legend.position = "None")

p <- (p1 + p2) / (p3 + p4)

ggsave("busco_painting.png", plot = p, width=10, height=10, units="in")
ggsave("busco_painting.pdf", plot = p, width=10, height=10, units="in")
