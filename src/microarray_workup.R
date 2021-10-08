#library(affy)
#library(mouse4302.db)
library(limma)
library(oligo)
library(biomaRt)
library(tidyverse)
## set up biomart
ensembl_mart <- useMart("ensembl")
ensembl_mart <- useDataset("mmusculus_gene_ensembl", mart = ensembl_mart)

cel_files <- list.files('data/',pattern = 'CEL', full.names = TRUE)
data <- read.celfiles(cel_files )

HNorm <- rma(data)
ID <- featureNames(HNorm)



gene_probe_table <- getBM(mart = ensembl_mart,
      values = ID,
      attributes = c("ensembl_gene_id", "affy_mouse430_2", "mgi_symbol"))

exprs(HNorm) %>% head()

sample_info <- read_tsv('data/sample_info.txt', col_names = FALSE)
sample_info %>% DT::datatable()
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
  left_join(sample_info) %>%
  ggplot(aes(x=PC1,y=PC2, color = X5)) +
  geom_point(size=4)

PCA$rotation %>%
  as_tibble(rownames = "affy_mouse430_2") %>%
  left_join(gene_probe_table) %>%
  dplyr::select(affy_mouse430_2, PC1, PC2, PC3, PC4, ensembl_gene_id, mgi_symbol) %>%
  arrange(-abs(PC2)) %>%
  head(20)
e
