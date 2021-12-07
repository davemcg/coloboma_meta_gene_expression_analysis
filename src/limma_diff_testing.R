library(edgeR)
library(sva)
library(tidyverse)

load('data/microarray_NGS_objects.Rdata')

norm_counts <- qsmooth_counts # works best at producing most padj in OF temporal testing (After - During)

sample_meta_D <- sample_meta %>% filter(Sample %in% colnames(norm_counts)) %>%
  dplyr::select(Sample:Section, Layout:Fusion) %>%
  mutate(S2 = case_when(Section == 'OF' ~ 'OF', TRUE ~ 'OC')) %>%
  unique()

# remove low expression genes from consideration
cutoff <- 3
drop <- which(apply((norm_counts), 1, max) < cutoff) #which(apply(cpm(d_zed), 1, max) < cutoff)
d <- norm_counts[-drop,]
dim(d) # number of genes left


sample_meta_filter <- sample_meta_D #%>% filter(S2 %in% 'OF')
d_filtered <- d[,sample_meta_filter$Sample]
colData <- colnames(d_filtered) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%  left_join(sample_meta_filter)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample

# build model, with organism and tech as covariates to regress out
# make fusion_Section (OF, OC) types for contrast in limma
colData$Fusion <- paste(colData$Fusion, colData$S2, sep = '_')
mm <- model.matrix(~0 + colData$Fusion + colData$Organism + colData$Technology)
colnames(mm) <- c('AfterOC', 'AfterOF','BeforeOC', 'BeforeOF','DuringOC','DuringOF', 'Mouse','Zebrafish','RNAseq')

# sva batch correction covariates (num_sv finds two for the OF datasets)
num_sv <- num.sv(d_filtered,mm,method="leek")
print(num_sv)
mod0 = model.matrix(~1, data=colData)
sv_obj <- sva(as.matrix(d_filtered), mm, mod0, n.sv=num_sv)

## (rough) batch corrected expression
### do NOT use for diff testing
### only (maybe) for vis?
sva_counts <- removeBatchEffect(d_filtered, covariates = sv_obj$sv)
write_tsv(sva_counts %>% as_tibble(rownames = 'Gene'), file = 'data/sva_counts.tsv.gz')

colnames(sv_obj$sv) <- c('SVA1','SVA2')
modSv = cbind(mm,sv_obj$sv)

fit <- lmFit(d_filtered, modSv)
contrast.matrix = makeContrasts(DuringOF-BeforeOF, AfterOF-DuringOF, DuringOC-DuringOF, levels=modSv)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)

###################################################################

efit <- eBayes(fit_contrasts)
top.table_OF_AD <- topTable(efit, sort.by = "p", n = Inf, coef="AfterOF - DuringOF", adjust.method="BH")
head(top.table_OF_AD, 20)
write_tsv(top.table_OF_AD %>% as_tibble(rownames = 'Gene'), file = 'data/top.table_OF_AD.tsv.gz')
top.table_OF_AD %>% filter(adj.P.Val<0.05) %>% dim()

top.table_OF_DB <- topTable(efit, sort.by = "p", n = Inf, coef="DuringOF - BeforeOF")
head(top.table_OF_DB, 20)
write_tsv(top.table_OF_DB %>% as_tibble(rownames = 'Gene'), file = 'data/top.table_OF_DB.tsv.gz')
top.table_OF_DB %>% filter(adj.P.Val<0.05) %>% dim()

top.table_During <- topTable(efit, sort.by = "p", n = Inf, coef="DuringOC - DuringOF", adjust.method="BH")
head(top.table_During, 20)
write_tsv(top.table_During %>% as_tibble(rownames = 'Gene'), file = 'data/top.table_OF_OC_During.tsv.gz')
top.table_During %>% filter(adj.P.Val<0.05) %>% dim()
####################################################################

save(top.table_OF_AD, top.table_OF_DB, top.table_During, file = 'data/top_tables.Rdata')

##########
# create counts by section | stage | paper | organism | tech
#########
sva_counts_merge_tidy <- sva_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(norm counts)') %>%
  left_join(sample_meta_D) %>%
  mutate(Group = glue::glue({'{S2}_{Fusion}_{Paper}_{Organism}_{Technology}'})) %>%
  group_by(Gene, Group) %>%
  summarise(`log2(norm counts)` = mean(`log2(norm counts)`))

sva_counts_merge <- sva_counts_merge_tidy %>%
  pivot_wider(names_from = Group, values_from = `log2(norm counts)`)
write_tsv(sva_counts_merge %>% as_tibble(rownames = 'Gene'),'data/sva_counts_merge.tsv.gz')
