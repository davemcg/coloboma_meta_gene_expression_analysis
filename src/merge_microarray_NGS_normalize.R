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

# load NGS
load('data/NGS_processing.Rdata')

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
# qsmooth_counts %>%
#   as_tibble(rownames = 'Gene') %>%
#   pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
#   #filter(Sample %in% colData$Sample) %>%
#   left_join(sample_meta) %>% ggplot(aes(x=`log2(Counts)`, color = Technology, group = Sample)) +
#   geom_density() +
#   facet_wrap(~Fusion) +
#   cowplot::theme_cowplot() +
#   geom_vline(xintercept = 3)
#####################################################################################

# straight quantile norm across the entire dataset
qnorm <- preprocessCore::normalize.quantiles(as.matrix(same_norm)) %>% data.frame()
colnames(qnorm) <- colnames(same_norm)
row.names(qnorm) <- row.names(same_norm)

qnorm_counts <- qnorm
colnames(qnorm_counts) <- colnames(qnorm_counts) %>% gsub('_.*|.CEL.*','',.)
# qnorm_counts %>%
#   as_tibble(rownames = 'Gene') %>%
#   pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
#   #filter(Sample %in% colData$Sample) %>%
#   left_join(sample_meta) %>% ggplot(aes(x=`log2(Counts)`, color = Technology, group = Sample)) +
#   geom_density() +
#   facet_wrap(~Fusion)
#####################################################


save(same_norm, same_rank_norm, same_qsmooth, qsmooth_counts, qnorm_counts, sample_meta, file = 'data/microarray_NGS_objects.Rdata')


