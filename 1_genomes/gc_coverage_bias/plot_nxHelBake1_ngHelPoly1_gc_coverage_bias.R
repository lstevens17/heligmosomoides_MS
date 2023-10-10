
# read in tables
HB_gc_cov <- read.table("nxHelBake1.1_gc_coverage.50N_filtered.tsv", col.names=c("chr", "start", "stop", "cov", "gc"))
HP_gc_cov <- read.table("ngHelPoly1.1_gc_coverage.50N_filtered.tsv", col.names=c("chr", "start", "stop", "cov", "gc"))

# get rid of NaN, caused by N containing-windows
HB_gc_cov <- HB_gc_cov %>% drop_na() %>% mutate(bin_dist = factor(gc%/%bin_size)) %>% mutate(window_size = stop - start) %>% filter(window_size == 100)
HP_gc_cov <- HP_gc_cov %>% drop_na() %>% mutate(bin_dist = factor(gc%/%bin_size)) %>% mutate(window_size = stop - start) %>% filter(window_size == 100)

# set bin size
bin_size <- 0.1

# plot box plot  of coverages in each bin
p1 <- HB_gc_cov %>% 
  ggplot(aes(x = bin_dist, y = cov)) +
  geom_boxplot(outlier.shape = NA, fill="#798cbd") + coord_cartesian(ylim=c(0,60)) +
  xlab("GC proportion") + ylab("Read coverage (x)") +  
  theme_bw()

# plot window count in each bin
p2 <- HB_gc_cov %>% 
  ggplot(aes(x = bin_dist)) + 
  geom_bar(fill="#798cbd") + 
  xlab("GC proportion") + ylab("Window count (n)") + 
  ggtitle(expression(paste(italic("H. bakeri"), " nxHelBake1"))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# plot box plot  of coverages in each bin
p3 <- HP_gc_cov %>% 
  ggplot(aes(x = bin_dist, y = cov)) +
  geom_boxplot(outlier.shape = NA, fill="#798cbd") + coord_cartesian(ylim=c(0,95)) +
  xlab("GC proportion") + ylab("Read coverage (x)") +  
  theme_bw()

# plot window count in each bin
p4 <- HP_gc_cov %>% 
  ggplot(aes(x = bin_dist)) + 
  geom_bar(fill="#798cbd") + 
  xlab("GC proportion") + ylab("Window count (n)") + 
  ggtitle(expression(paste(italic("H. polygyrus"), " ngHelPoly1"))) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

design <- "AAA
           BBB
           BBB
           BBB
           CCC
           DDD
           DDD
           DDD"

p <- p2 + p1 + p4 + p3 + plot_layout(design = design)

ggsave("nxHelBake1_ngHelPoly1_gc_coverage_bias.pdf", plot = p, width=8, height=10, units="in")
ggsave("nxHelBake1_ngHelPoly1_gc_coverage_bias.png", plot = p, width=8, height=10, units="in")