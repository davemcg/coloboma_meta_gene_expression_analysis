library(edgeR)
library(sva)
ngs_counts <- qsmooth_counts # works best at producing most padj in OF temporal testing (After - During)

sample_meta_D <- sample_meta %>% filter(Sample %in% colnames(ngs_counts)) %>%
  dplyr::select(Sample:Section, Layout:Fusion) %>%
  mutate(S2 = case_when(Section == 'OF' ~ 'OF', TRUE ~ 'Retina')) %>%
  unique()

# remove low expression genes from consideration
cutoff <- 3
drop <- which(apply((ngs_counts), 1, max) < cutoff) #which(apply(cpm(d_zed), 1, max) < cutoff)
d <- ngs_counts[-drop,]
dim(d) # number of genes left

# filter down to OF (optic fissure) samples
sample_meta_filter <- sample_meta_D %>% filter(S2 %in% 'OF')
d_filtered <- d[,sample_meta_filter$Sample]
colData <- colnames(d_filtered) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%  left_join(sample_meta_filter)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample

# build model, with organism and tech as covariates to regress out
mm <- model.matrix(~0 + colData$Fusion + colData$Organism + colData$Technology)
colnames(mm) <- c('After', 'Before','During', 'Mouse','Zebrafish','RNAseq')

# sva batch correction covariates (num_sv finds two for the OF datasets)
num_sv <- num.sv(d_filtered,mm,method="leek")
print(num_sv)
mod0 = model.matrix(~1, data=colData)
sv_obj <- svaseq(as.matrix(d_filtered), mm, mod0, n.sv=num_sv)

## rough batch corrected expression
cor_vals_OF <- removeBatchEffect(d_filtered, covariates = sv_obj$sv)

modSv = cbind(mm,sv_obj$sv)
colnames(modSv)[7:8] <- c('SVA1','SVA2')

fit <- lmFit(d_filtered, modSv)
contrast.matrix = makeContrasts(During-Before, After-During,levels=modSv)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)

###################################################################

efit <- eBayes(fit_contrasts)
top.table_OF_AD <- topTable(efit, sort.by = "p", n = Inf, coef="After - During", adjust.method="BH")
head(top.table_OF_AD, 20)
top.table_OF_AD %>% filter(adj.P.Val<0.05) %>% dim()

top.table_OF_DB <- topTable(efit, sort.by = "p", n = Inf, coef="During - Before")
head(top.table_OF_DB, 20)
####################################################################



# AGAIN!
# filter down to retina samples (not the OF)
sample_meta_filter <- sample_meta_D %>% filter(!S2 %in% 'OF')
d_filtered <- d[,sample_meta_filter$Sample]
colData <- colnames(d_filtered) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%  left_join(sample_meta_filter)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample

# build model, with organism and tech as covariates to regress out
mm <- model.matrix(~0 + colData$Fusion + colData$Organism + colData$Technology)
colnames(mm) <- c('After', 'Before','During', 'Mouse','Zebrafish','RNAseq')

# sva batch correction covariates (num_sv finds three for the retina datasets)
num_sv <- num.sv(d_filtered,mm,method="leek")
print(num_sv)
mod0 = model.matrix(~1, data=colData)
sv_obj <- svaseq(as.matrix(d_filtered), mm, mod0, n.sv=num_sv)

## rough batch corrected expression
cor_vals_retina <- removeBatchEffect(d_filtered, covariates = sv_obj$sv)

modSv = cbind(mm,sv_obj$sv)
colnames(modSv)[7:9] <- c('SVA1','SVA2', 'SVA3')

fit <- lmFit(d_filtered, modSv)
contrast.matrix = makeContrasts(During-Before, After-During,levels=modSv)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)

###################################################################

efit <- eBayes(fit_contrasts)
top.table_RETINA_AD <- topTable(efit, sort.by = "p", n = Inf, coef="After - During", adjust.method="BH")
head(top.table_RETINA_AD, 20)
top.table_RETINA_AD %>% filter(adj.P.Val<0.05) %>% dim()

top.table_RETINA_DB <- topTable(efit, sort.by = "p", n = Inf, coef="During - Before")
head(top.table_RETINA_DB, 20)
####################################################################








qsmooth_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'Expression') %>%
  mutate(Sample = gsub('_.*|.CEL.*','',Sample)) %>%
  left_join(sample_meta_D) %>%
  mutate(S2 = case_when(Section == 'OF' ~ 'OF',TRUE ~ 'Retina')) %>%
  #filter(Gene %in% c('CRYM')) %>%
  filter(Gene %in% row.names(top.table_OF_DB %>% head(4))) %>%
  mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
  ggplot(aes(x=Fusion, y=Expression, color = Organism, shape = Technology)) +
  # geom_boxplot(aes(group = Fusion), color = 'Black', outlier.colour = NA) +
  # geom_boxplot(aes(group = interaction(Organism,Fusion)), outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom(size = 3, alpha = 0.7) +
  cowplot::theme_cowplot() +
  facet_grid(~Gene + S2, scales = 'free_y') +
  ggsci::scale_color_aaas() +
  ylab('log2 (qn counts)') +
  stat_summary(fun=mean, geom = 'line', aes(group = Organism))



# volcano
volcano_maker <- function(df, title="Volcano Plot", pvalue='P.Value', padj='adj.P.Val', logFC='logFC', gene_list = ''){
  df$pvalue <- df[,pvalue]
  df$log2FoldChange <- df[,logFC]
  df$padj <- df[,padj]
  df$Gene <- row.names(df)
  df <- df[!is.na(df$pvalue),]
  print(dim(df))

  df <- df %>% mutate(Class = case_when(padj < 0.05 & abs(logFC) > 1~ "FDR < 0.05 & logFC > 1",
                                        padj < 0.1 & abs(logFC) > 1 ~ 'FDR < 0.1 & logFC > 1',
                                        TRUE ~ 'Not significant'))
  df$GeneT <- df$Gene

  df$Gene[!df$Gene %in% gene_list] <- ''

  plot <- ggplot(data=df,aes(label=Gene, x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(colour=Class)) +
    scale_colour_manual(values=c("darkred", "red", "grey")) +
    cowplot::theme_cowplot() +
    geom_vline(aes(xintercept=-1),linetype="dotted") +
    geom_vline(aes(xintercept=1),linetype="dotted") +
    geom_vline(aes(xintercept=-2),linetype="dotted") +
    geom_vline(aes(xintercept=2),linetype="dotted") +
    geom_label_repel(data=bind_rows(subset(df, padj<0.01)),
                     aes(label=GeneT)) +
    xlab('logFC') + ylab('-log10(p value)') +
    ggtitle(title) + cowplot::theme_cowplot()

    plot
}

volcano_maker(top.table_OF_AD)



# genes that are diff in After - During exclusive (by padj) in OF
retina_AD_diff_genes <- top.table_RETINA_AD %>% as_tibble(rownames = 'Gene') %>% filter(adj.P.Val < 0.05) %>% pull(Gene)
## OF specific
top.table_OF_AD %>% as_tibble(rownames = 'Gene') %>% filter(adj.P.Val < 0.05, !Gene %in% retina_AD_diff_genes)
## shared
top.table_OF_AD %>% as_tibble(rownames = 'Gene') %>% filter(adj.P.Val < 0.05, Gene %in% retina_AD_diff_genes)


# genes that are diff in During - Before exclusive (by padj) in OF
retina_DB_diff_genes <- top.table_RETINA_DB %>% as_tibble(rownames = 'Gene') %>% filter(adj.P.Val < 0.05) %>% pull(Gene)
top.table_OF_DB %>% as_tibble(rownames = 'Gene') %>% filter(adj.P.Val < 0.05, !Gene %in% retina_DB_diff_genes)
top.table_OF_DB %>% as_tibble(rownames = 'Gene') %>% filter(adj.P.Val < 0.05, Gene %in% retina_DB_diff_genes)
