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

register(MulticoreParam(2))
#working_dir <- '/Volumes/ARC168/PROJECTS/hufnagel/macaque_fovea_RNA-seq/salmon_quant/'
working_dir <- '~/data/coloboma_meta_gene_expression/salmon_quant/'
files <- list.files(path=working_dir,recursive=TRUE,pattern='quant.sf', full.names = TRUE)


# adjust metadata to classify stages in three categories:
## before fusion
## during fusion
## after fusion
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
                            Stage == '56hpf' ~ 'After'),
         Section = gsub('OFM','OF',Section))
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
                            Stage == '56hpf' ~ 'After'),
         Section = gsub('OFM','OF',Section))

rna_processor <- function(meta, organism, countsFromAbundance = 'no'){

  samples <- str_extract(files, 'SAMEA\\d+|SRS\\d+')
  sample_meta_rnaseq_org <- meta %>% filter(Sample %in% samples) %>% filter(Organism == organism)

  org_files <- files %>% enframe() %>%
    mutate(sample = str_extract(value,  'SAMEA\\d+|SRS\\d+')) %>%
    filter(sample %in% sample_meta_rnaseq_org$Sample) %>%
    pull(value)
  anno <- fread(org_files[1])
  # zebrafish fasta not from gencode (they don't make)
  ## so have to do a custom data ingress to get the tx <-> gene name table
  if (organism == 'Zebrafish'){
    zf_tx <- read_tsv('~/git/coloboma_meta_gene_expression_analysis/data/Danio_rerio.GRCz11.release104.cdna.info.txt.gz', col_names = FALSE)
    zf_tx_info <- zf_tx %>% separate(X1, into = c('tx','type','chromosome','gene','gene_biotype','transcript_biotype', 'gene_symbol', 'description'), sep =' ') %>%
      mutate(tx = str_extract(tx, 'ENSDART\\d+.\\d+'),
             gene = gsub('gene:','',gene),
             gene_symbol = gsub('gene_symbol:','',gene_symbol))
    anno <- anno %>% left_join(zf_tx_info %>% dplyr::select(tx, gene), by = c('Name' = 'tx')) %>%
      dplyr::rename(Gene = gene)
  } else {
    anno$Gene <- sapply(anno$Name,function(x) strsplit(x,'\\|')[[1]][6])
  }
  anno_tximport <- anno %>%
    dplyr::select(target_id = Name, Gene)


  txi <- tximport(org_files, type = "salmon", tx2gene = anno_tximport, countsFromAbundance = countsFromAbundance)
  txi.deseq2 <- data.frame(txi$counts)

  out <- list()
  out$txi <- txi
  out$counts <- txi.deseq2
  out$meta <- sample_meta_rnaseq_org
  colnames(out$counts) <- out$meta$Sample %>% unique()
  out
}

human <- rna_processor(sample_meta_rnaseq, 'Human', 'scaledTPM')
mouse <- rna_processor(sample_meta_rnaseq, 'Mouse', 'scaledTPM')
zebrafish <- rna_processor(sample_meta_rnaseq, 'Zebrafish', 'scaledTPM')

# merge mouse and human data
same <- human$counts %>% as_tibble(rownames = 'Gene') %>%
  left_join(mouse$counts %>% as_tibble(rownames = 'Gene') %>%
              mutate(Gene = toupper(Gene))) %>%
  data.frame()
row.names(same) <- same$Gene
same <- same[,-1] %>% data.frame()
same <- same[complete.cases(same),] # remove NA

########################################################
# now merge in zebrafish data
## grab human <-> zf gene table
hcop <- read_tsv('data/human_zebrafish_hcop_fifteen_column.txt.gz')
## collapse duplicate ZF gene names by summing
## e.g. abca4a and abca4b are all combined together into abca4
## that's what the group_by -> summarise steps are doing
zf_counts_human_gene_names <- zebrafish$counts %>%
  as_tibble(rownames = 'zebrafish_ensembl_gene') %>%
  mutate(zebrafish_ensembl_gene = gsub('\\.\\d+','', zebrafish_ensembl_gene)) %>%
  left_join(hcop %>% dplyr::select(human_symbol, zebrafish_ensembl_gene) %>%
              filter(human_symbol != '-')) %>%
  group_by(human_symbol) %>%
  summarise(across(where(is.numeric), sum))

## now join the human/mouse zebrafish data
same$Gene <- row.names(same)
same2 <- same %>% left_join(zf_counts_human_gene_names, by = c('Gene' = 'human_symbol'))
row.names(same2) <- same2$Gene
same2 <- same2[,!grepl('Gene',colnames(same2))] %>% data.frame()
same2 <- same2[complete.cases(same2),]

#####################################################################

same_log <- log2(same2 + 1)

## normalize with voom (replaces log norm)
# colData <- colnames(same2) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%
#   left_join(sample_meta_rnaseq %>% dplyr::select(Sample:Fusion, -Other) %>% unique())
# colData <- data.frame(colData)
# row.names(colData) <- colData$Sample
# colData$Section <- gsub('OFM','OF',colData$Section)
#
# mm <- model.matrix(~0  + colData$Fusion  + colData$Section + colData$Organism)
# colnames(mm) <- c('After','Before','During','DR','OF','Mouse','Zebrafish')
# y <- voom(same2, mm, plot = T)
#
# same_voom <- y$E

same_norm <- same_log

# add microarray
load('data/microarray_table.Rdata')
colnames(microarray_table) <- colnames(microarray_table) %>% gsub('_.*|.CEL.*','',.)
same_norm <- same_norm %>% as_tibble(rownames = 'Gene') %>%
  left_join(microarray_table %>% as_tibble(rownames = 'Gene') %>%
              mutate(Gene = toupper(Gene)))

same_norm <- same_norm %>% data.frame()
row.names(same_norm) <- same_norm$Gene
same_norm <- same_norm[,-1]
same_norm <- same_norm[complete.cases(same_norm),]
################################################################################



# rank norm
same_rank_norm <- apply(same_norm, 2, function(y) rank(y) / length(y))
########################################

# qsmooth norm
## smooth across each fusion / setion (OF vs not OF)
qsmooth_factors <- same_norm %>%
  colnames() %>%
  enframe() %>%
  mutate(Sample = case_when(grepl('CEL', value) ~ str_extract(value, 'GSM\\d+'),
                            TRUE ~ value)) %>%
  left_join(sample_meta %>% dplyr::select(Sample, Accession, Organism, Fusion, Section) %>% unique(),
            by = 'Sample') %>%
  mutate(S2 = case_when(Section == 'OF' ~ 'OF', TRUE ~ 'Retina'),
         FS = paste(Fusion, S2, sep = '_')) %>%
  pull(FS)
# qsmooth_batch <- same_norm %>%
#   colnames() %>%
#   enframe() %>%
#   mutate(Sample = case_when(grepl('CEL', value) ~ str_extract(value, 'GSM\\d+'),
#                             TRUE ~ value)) %>%
#   left_join(sample_meta %>% dplyr::select(Sample, Accession, Fusion, Technology) %>% unique(),
#             by = 'Sample') %>%
#   pull(Accession)
# same_qsmooth <- qsmooth(same_norm, group_factor = qsmooth_factors, batch = qsmooth_batch)
same_qsmooth <- qsmooth(same_norm, group_factor = qsmooth_factors)

qsmooth_counts <- same_qsmooth@qsmoothData
colnames(qsmooth_counts) <- colnames(qsmooth_counts) %>% gsub('_.*|.CEL.*','',.)
qsmooth_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
  #filter(Sample %in% colData$Sample) %>%
  left_join(sample_meta) %>% ggplot(aes(x=`log2(Counts)`, color = Technology, group = Sample)) +
  geom_density() +
  facet_wrap(~Fusion) +
  cowplot::theme_cowplot() +
  geom_vline(xintercept = 3)
#####################################################################################

# straight quantile norm across the entire dataset
qnorm <- preprocessCore::normalize.quantiles(as.matrix(same_norm)) %>% data.frame()
colnames(qnorm) <- colnames(same_norm)
row.names(qnorm) <- row.names(same_norm)

qnorm_counts <- qnorm
colnames(qnorm_counts) <- colnames(qnorm_counts) %>% gsub('_.*|.CEL.*','',.)
qnorm_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
  #filter(Sample %in% colData$Sample) %>%
  left_join(sample_meta) %>% ggplot(aes(x=`log2(Counts)`, color = Technology, group = Sample)) +
  geom_density() +
  facet_wrap(~Fusion)
#####################################################



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

save(same, same_norm, same_rank_norm, same_qsmooth, qnorm_counts, file = 'data/microarray_NGS_objects.Rdata')


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

