library(data.table)
library(tidyverse)
library(biomaRt)
library(tximport)
library(DESeq2)
library(ggrepel)
library(apeglm)
# Parallel
library(BiocParallel)
library(matrixStats)
library(cowplot)
library(qsmooth)

# OF Translational Time Points
# Info		Human	Mouse	Chick	Zebrafish
# Before Fusion	Before OFM fusion	CS15 (33)	E10.5	E5	36hpf
# During Fusion	Fusion @ Single Point	CS16 (37)	E11.5	E6	48hpf
# After Fusion	Completely Fused	CS18 (44)	E12.5	E7	56hpf

register(MulticoreParam(6))
#working_dir <- '/Volumes/ARC168/PROJECTS/hufnagel/macaque_fovea_RNA-seq/salmon_quant/'
working_dir <- '~/data/coloboma_meta_gene_expression/salmon_quant/'
files <- list.files(path=working_dir,recursive=TRUE,pattern='quant.sf', full.names = TRUE)

sample_meta <- read_tsv('data/sample_info2.tsv') %>%
  mutate(Paper = case_when(grepl('MTAB', Accession) ~ 'Patel - Sowden',
                           Accession == 'GSE109501' ~ 'Cao - Chen',
                           Accession == 'GSE13103' ~ 'Brown - Brooks',
                           Accession == 'GSE159822' ~ 'Richardson - Moosajee',
                           Accession == 'GSE84916' ~ 'Hardy - Brooks'),
         Technology = case_when(is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq'),
         Fusion = case_when(Stage == 'CS17' ~ 'During',
                            Stage == 'CS18' ~ 'After',
                            Stage == 'E10.5' ~ 'Before',
                            Stage == 'E11.5' ~ 'During',
                            Stage == 'E12.5' ~ 'After',
                            Stage == 'E5' ~ 'Before',
                            Stage == 'E6' ~ 'During',
                            Stage == 'E7' ~ 'After',
                            Stage == '32hpf' ~ 'Before',
                            Stage == '48hpf' ~ 'During',
                            Stage == '56hpf' ~ 'After'))
sample_meta_rnaseq <- read_tsv('data/sample_info2.tsv') %>% filter(!is.na(Run))%>%
  mutate(Paper = case_when(grepl('MTAB', Accession) ~ 'Patel - Sowden',
                           Accession == 'GSE109501' ~ 'Cao - Chen',
                           Accession == 'GSE13103' ~ 'Brown - Brooks',
                           Accession == 'GSE159822' ~ 'Richardson - Moosajee',
                           Accession == 'GSE84916' ~ 'Hardy - Brooks'),
         Technology = case_when(!is.na(Run) ~ 'Microarray', TRUE ~ 'RNA-Seq'),
         Fusion = case_when(Stage == 'CS17' ~ 'During',
                            Stage == 'CS18' ~ 'After',
                            Stage == 'E10.5' ~ 'Before',
                            Stage == 'E11.5' ~ 'During',
                            Stage == 'E12.5' ~ 'After',
                            Stage == 'E5' ~ 'Before',
                            Stage == 'E6' ~ 'During',
                            Stage == 'E7' ~ 'After',
                            Stage == '32hpf' ~ 'Before',
                            Stage == '48hpf' ~ 'During',
                            Stage == '56hpf' ~ 'After'))

rna_processor <- function(meta, organism, countsFromAbundance = 'no'){

  samples <- str_extract(files, 'SAMEA\\d+')
  sample_meta_rnaseq_org <- meta %>% filter(Sample %in% samples) %>% filter(Organism == organism)

  org_files <- files %>% enframe() %>%
    mutate(sample = str_extract(value, 'SAMEA\\d+')) %>%
    filter(sample %in% sample_meta_rnaseq_org$Sample) %>%
    pull(value)
  anno <- fread(org_files[1])
  anno$Gene <- sapply(anno$Name,function(x) strsplit(x,'\\|')[[1]][6])
  anno_tximport <- anno %>%
    dplyr::select(target_id = Name, Gene)


  txi <- tximport(org_files, type = "salmon", tx2gene = anno_tximport, countsFromAbundance = countsFromAbundance)
  txi.deseq2 <- data.frame(txi$counts)

  out <- list()
  out$txi <- txi
  out$counts <- txi.deseq2
  out$meta <- sample_meta_rnaseq_org
  out
}

human <- rna_processor(sample_meta_rnaseq, 'Human', 'lengthScaledTPM')
mouse <- rna_processor(sample_meta_rnaseq, 'Mouse', 'lengthScaledTPM')

colnames(human$counts) <- human$meta$Sample %>% unique()
colnames(mouse$counts) <- mouse$meta$Sample %>% unique()

same <- human$counts %>% as_tibble(rownames = 'Gene') %>%
  left_join(mouse$counts %>% as_tibble(rownames = 'Gene') %>%
              mutate(Gene = toupper(Gene))) %>%
  data.frame()
row.names(same) <- same$Gene
same <- same[,-1] %>% data.frame()
same <- same[complete.cases(same),]

same_log <- log2(same + 1)
# add microarray
load('data/microarray_table.Rdata')
same_log <- same_log %>% as_tibble(rownames = 'Gene') %>%
  left_join(microarray_table %>% as_tibble(rownames = 'Gene') %>%
              mutate(Gene = toupper(Gene)))
row.names(same_log) <- same_log$Gene
same_log <- same_log[,-1] %>% data.frame()
same_log <- same_log[complete.cases(same_log),]

# rank norm
same_rank_norm <- apply(same_log, 2, function(y) rank(y) / length(y))

# qsmooth norm
## smooth across species for now
qsmooth_factors <- same_log %>%
  colnames() %>%
  enframe() %>%
  mutate(Sample = case_when(grepl('CEL', value) ~ str_extract(value, 'GSM\\d+'),
                            TRUE ~ value)) %>%
  left_join(sample_meta %>% dplyr::select(Sample, Accession, Organism, Fusion) %>% unique(),
            by = 'Sample') %>%
  pull(Fusion)
qsmooth_batch <- same_log %>%
  colnames() %>%
  enframe() %>%
  mutate(Sample = case_when(grepl('CEL', value) ~ str_extract(value, 'GSM\\d+'),
                            TRUE ~ value)) %>%
  left_join(sample_meta %>% dplyr::select(Sample, Accession, Fusion) %>% unique(),
            by = 'Sample') %>%
  pull(Accession)

same_qsmooth <- qsmooth(same_log, group_factor = qsmooth_factors, batch = qsmooth_batch)
#####################################################################################


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

PCA_log <- run_PCA(same_log, n_top_var = 2000)
PCA_rank <- run_PCA(same_rank_norm, n_top_var = 2000)
PCA_qsmooth <- run_PCA(same_qsmooth@qsmoothData, n_top_var = 2000)

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
  merge <- cowplot::plot_grid(one, two, twoB,
                              three, four, fourB,
                              five, six, sixB,
                              seven, eight, eightB,
                              ncol = 3)
  plot_grid(title, merge, ncol = 1, rel_heights = c(0.1,5))
}

system('mkdir -p analysis')
pdf('analysis/PCA_noNorm.pdf', width = 12, height = 12)
plotter(PCA_log, plot_title = 'No normalization')
dev.off()

pdf('analysis/PCA_rankNorm.pdf', width = 12, height = 12)
plotter(PCA_rank, plot_title = 'Rank normalization')
dev.off()

pdf('analysis/PCA_qsmoothNorm.pdf', width = 12, height = 12)
plotter(PCA_qsmooth, plot_title = 'qsmooth normalization')
dev.off()

