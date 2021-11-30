#library(affy)
#library(mouse4302.db)
library(limma)
library(oligo)
library(biomaRt)
library(tidyverse)
## set up biomart
# ensembl_mart <- useMart("ensembl")
# ensembl_mart <- useDataset("mmusculus_gene_ensembl", mart = ensembl_mart)
#
# gene_probe_table <- getBM(mart = ensembl_mart,
#                           values = ID,
#                           attributes = c("ensembl_gene_id", "affy_mouse430_2", "mgi_symbol"))

gene_probe_table <- read_tsv('data/gene_mouse_probe.tsv.gz')
cel_files <- list.files('data/',pattern = 'CEL', full.names = TRUE)
data <- affy::ReadAffy(filenames = cel_files )

HNorm <- gcrma::gcrma(data)
ID <- featureNames(HNorm)




exprs(HNorm) %>% head()

sample_info <- read_tsv('data/sample_info2.tsv') %>% filter(is.na(Run)) %>% dplyr::select(CEL, Sample:Other)
sample_info %>% DT::datatable()

study_info <- read_tsv('data/study_info_2.tsv')
# PCA
library(matrixStats)
ntop = 1000
Pvars <- exprs(HNorm)
select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop,
                                                      length(Pvars)))]
PCA <- prcomp(t(Pvars), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

PCA$x %>% as_tibble(rownames = 'X1') %>%
  mutate() %>%
  left_join(sample_info, by = c('X1' = 'CEL')) %>%
  ggplot(aes(x=PC1,y=PC2, color = Stage)) +
  geom_point(size=4)

PCA$x %>% as_tibble(rownames = 'X1') %>%
  mutate() %>%
  left_join(sample_info, by = c('X1' = 'CEL')) %>%
  ggplot(aes(x=PC3,y=PC4, color = Stage)) +
  geom_point(size=4)

PCA$rotation %>%
  as_tibble(rownames = "affy_mouse430_2") %>%
  left_join(gene_probe_table) %>%
  dplyr::select(affy_mouse430_2, PC1, PC2, PC3, PC4, ensembl_gene_id, mgi_symbol) %>%
  arrange(-abs(PC2)) %>%
  head(20)


microarray_table <-
  exprs(HNorm) %>%
  as_tibble(rownames = 'affy_mouse430_2') %>%
  left_join(gene_probe_table) %>%
  filter(!is.na(mgi_symbol)) %>%
  group_by(mgi_symbol) %>%
  summarise(`GSM2944692_E11.5-OF-A.CEL.gz` = sum(`GSM2944692_E11.5-OF-A.CEL.gz`),
            `GSM2944693_E11.5-OF-B.CEL.gz` = sum(`GSM2944693_E11.5-OF-B.CEL.gz`),
            `GSM2944694_E11.5-OF-C.CEL.gz` = sum(`GSM2944694_E11.5-OF-C.CEL.gz`),
            `GSM2944695_E11.5-NR-A.CEL.gz` = sum(`GSM2944695_E11.5-NR-A.CEL.gz`),
            `GSM2944696_E11.5-NR-B.CEL.gz` = sum(`GSM2944696_E11.5-NR-B.CEL.gz`),
            `GSM2944697_E11.5-NR-C.CEL.gz` = sum(`GSM2944697_E11.5-NR-C.CEL.gz`),
            `GSM2944698_E11.5-TR-A.CEL.gz` = sum(`GSM2944698_E11.5-TR-A.CEL.gz`),
            `GSM2944699_E11.5-TR-B.CEL.gz` = sum(`GSM2944699_E11.5-TR-B.CEL.gz`),
            `GSM2944700_E11.5-TR-C.CEL.gz` = sum(`GSM2944700_E11.5-TR-C.CEL.gz`),
            `GSM327854.CEL.gz` = sum(`GSM327854.CEL.gz`),
            `GSM327855.CEL.gz` = sum(`GSM327855.CEL.gz`),
            `GSM327856.CEL.gz` = sum(`GSM327856.CEL.gz`),
            `GSM327857.CEL.gz` = sum(`GSM327857.CEL.gz`),
            `GSM327858.CEL.gz` = sum(`GSM327858.CEL.gz`),
            `GSM327859.CEL.gz` = sum(`GSM327859.CEL.gz`),
            `GSM327860.CEL.gz` = sum(`GSM327860.CEL.gz`),
            `GSM327861.CEL.gz` = sum(`GSM327861.CEL.gz`)) %>%
  data.frame()
row.names(microarray_table) <- microarray_table$mgi_symbol
microarray_table <- microarray_table[,-1]
save(microarray_table,  file = 'data/microarray_table.Rdata')
