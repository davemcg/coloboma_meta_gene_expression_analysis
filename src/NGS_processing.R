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

save(same, same_norm,sample_meta, file = 'data/NGS_processing.Rdata')
