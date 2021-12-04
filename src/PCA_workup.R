
run_PCA <- function(matrix, n_top_var = 2000){
  ntop = n_top_var
  same_vars <- rowVars(as.matrix(matrix))
  select <- order(same_vars, decreasing = TRUE)[seq_len(min(ntop,
                                                            length(same_vars)))]
  matrix_select <- matrix[select,]
  PCA <- prcomp(t(matrix_select), scale = F)
  #PCA$percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
  PCA
}

PCA_log <- run_PCA(same_norm[,sample_meta %>% filter(Section == 'OF') %>% pull(Sample) %>% unique()], n_top_var = 2000)
PCA_rank <- run_PCA(same_rank_norm[,sample_meta %>% filter(Section == 'OF') %>% pull(Sample) %>% unique()], n_top_var = 2000)
PCA_qsmooth <- run_PCA(ngs_counts[,sample_meta %>% filter(Section == 'OF') %>% pull(Sample) %>% unique()], n_top_var = 2000)
PCA_qnorm <- run_PCA(qnorm_counts[,sample_meta %>% filter(Section == 'OF') %>% pull(Sample) %>% unique()], n_top_var =2000)
#PCA_sva <- run_PCA(same_sva_norm, n_top_var = 2000)


plotter <- function(PCA, plot_title = NA){
  title <- ggdraw() +
    draw_label(
      plot_title,
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  # PC 1 2
  one <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta,
              by = 'Sample') %>%
    ggplot(aes(x=PC1,y=PC2, color = Organism, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot()

  two <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC1,y=PC2, color = Technology, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot()+
    ggsci::scale_color_aaas()

  twoB <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC1,y=PC2, color = Paper, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_jama()

  twoC <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC1,y=PC2, color = Fusion, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_d3()

  twoD <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Section = gsub('OFM','OF',Section)),
              by = 'Sample') %>%
    ggplot(aes(x=PC1,y=PC2, color = Section, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_futurama()

  # PC 3 4
  three <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC3,y=PC4, color = Organism, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot()

  four <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC3,y=PC4, color = Technology, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_aaas()

  fourB <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC3,y=PC4, color = Paper, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_jama()

  fourC <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC3,y=PC4, color = Fusion, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_d3()

  fourD <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Section = gsub('OFM','OF',Section)),
              by = 'Sample') %>%
    ggplot(aes(x=PC3,y=PC4, color = Section, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_futurama()

  # PC 5 6
  five <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC5,y=PC6, color = Organism, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot()

  six <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC5,y=PC6, color = Technology, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_aaas()

  sixB <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC5,y=PC6, color = Paper, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_jama()

  sixC <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC5,y=PC6, color = Fusion, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_d3()

  sixD <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Section = gsub('OFM','OF',Section)),
              by = 'Sample') %>%
    ggplot(aes(x=PC5,y=PC6, color = Section, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_futurama()

  # PC 7 8
  seven <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC7,y=PC8, color = Organism, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot()

  eight <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC7,y=PC8, color = Technology, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_aaas()

  eightB <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC7,y=PC8, color = Paper, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_jama()

  eightC <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq')),
              by = 'Sample') %>%
    ggplot(aes(x=PC7,y=PC8, color = Fusion, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_d3()

  eightD <- PCA$x %>% as_tibble(rownames = 'Sample') %>%
    mutate(Sample = case_when(grepl('CEL', Sample) ~ str_extract(Sample, 'GSM\\d+'),
                              TRUE ~ Sample)) %>%
    left_join(sample_meta %>% mutate(Section = gsub('OFM','OF',Section)),
              by = 'Sample') %>%
    ggplot(aes(x=PC7,y=PC8, color = Section, shape = Fusion)) +
    geom_point(size=4) +
    cowplot::theme_cowplot() +
    ggsci::scale_color_futurama()

  merge <- cowplot::plot_grid(one, two, twoB, twoC, twoD,
                              three, four, fourB, fourC, fourD,
                              five, six, sixB, sixC, sixD,
                              seven, eight, eightB, eightC, eightD,
                              ncol = 5)
  plot_grid(title, merge, ncol = 1, rel_heights = c(0.1,5))
}

system('mkdir -p analysis')
pdf('analysis/PCA_noNorm.pdf', width = 24, height = 12)
plotter(PCA_log, plot_title = 'No normalization')
dev.off()

pdf('analysis/PCA_rankNorm.pdf', width = 24, height = 12)
plotter(PCA_rank, plot_title = 'Rank normalization')
dev.off()

pdf('analysis/PCA_qsmoothNorm.pdf', width = 24, height = 12)
plotter(PCA_qsmooth, plot_title = 'qsmooth normalization')
dev.off()

pdf('analysis/PCA_qnorm.pdf', width = 24, height = 12)
plotter(PCA_qnorm, plot_title = 'Quantile Normalization')
dev.off()
