library(tidyverse)
library(cowplot)
load('data/microarray_NGS_objects.Rdata')
c <- qsmooth_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
  #filter(Sample %in% colData$Sample) %>%
  left_join(sample_meta) %>%
  mutate(Dataset = glue::glue("{Paper} ({Organism})")) %>%
  mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
  ggplot(aes(x=`log2(Counts)`, y = Dataset, fill = Technology)) +
  ggridges::geom_density_ridges() +
  #geom_density(alpha = 0.2, size = 1) +
  facet_wrap(~Fusion) +
  cowplot::theme_cowplot() +
  xlab('log2(quantile norm counts)') +
  geom_vline(xintercept = 3)


b <- same_norm %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
  #filter(Sample %in% colData$Sample) %>%
  left_join(sample_meta) %>%
  mutate(Dataset = glue::glue("{Paper} ({Organism})")) %>%
  mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
  ggplot(aes(x=`log2(Counts)`, y = Dataset, fill = Technology)) +
  ggridges::geom_density_ridges() +
  #geom_density(alpha = 0.2, size = 1) +
  facet_wrap(~Fusion) +
  cowplot::theme_cowplot() +
  xlab('log2(counts)')
  #geom_vline(xintercept = 3)


# a <- same_norm %>%
#   as_tibble(rownames = 'Gene') %>%
#   pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
#   #filter(Sample %in% colData$Sample) %>%
#   left_join(sample_meta) %>% ggplot(aes(x=`log2(Counts)`, group = Sample, color = Technology)) +
#   #ggridges::geom_density_ridges() +
#   geom_density(alpha = 0.2, size = 0.2) +
#   facet_wrap(~Fusion) +
#   cowplot::theme_cowplot() +
#   xlab('log2(quantile norm counts)') +
#   geom_vline(xintercept = 3)

svg(filename = 'figs/supFig01.svg', width = 8, height = 8)
plot_grid(b,c, align = 'hv', rel_heights = c(3,3), labels = 'auto', ncol = 1)
dev.off()
