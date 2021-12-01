library(edgeR)
ngs_counts <- qsmooth_counts
colnames(ngs_counts) <- colnames(ngs_counts) %>% gsub('_.*|.CEL.*','',.)
#d_zed <- DGEList(ngs_counts)
#d_zed <- calcNormFactors(d_zed)

sample_meta_D <- sample_meta %>% filter(Sample %in% colnames(ngs_counts)) %>%
  dplyr::select(Sample:Section, Layout:Fusion) %>%
  unique()

colData <- colnames(ngs_counts) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%  left_join(sample_meta_D)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample

# d_zed$samples <- colData
#
cutoff <- 1
drop <- which(apply((ngs_counts), 1, max) < cutoff) #which(apply(cpm(d_zed), 1, max) < cutoff)
d <- ngs_counts[-drop,]
dim(d) # number of genes left
d <- ngs_counts


run_limma <- function(counts_matrix,
                      sample_meta_filtered,
                      model = )

# reduce down to OF and OFM (optic fissure (margin))
# ofm_meta <- sample_meta_D %>% filter(Section %in% c('OF','OFM'))
counts_filtered <- counts[,sample_meta_filter$Sample]
colData <- colnames(counts_filtered) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%  left_join(sample_meta_filter)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample


mm <- model.matrix(~0  + colData$Fusion + colData$Organism + colData$Technology)
colnames(mm) <- c('After', 'Before','During', 'Mouse','Zebrafish','RNAseq')
#y <- voom(d, mm, plot = T)
fit <- lmFit(d_ofm, mm)
contrast.matrix = makeContrasts(During-Before, After-During,levels=mm)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)



efit <- eBayes(fit_contrasts, proportion = 0.1)
top.table <- topTable(efit, sort.by = "p", n = Inf, coef="After - During", adjust.method="BH")
head(top.table, 20)

top.table <- topTable(efit, sort.by = "p", n = Inf, coef="During - Before")
head(top.table, 20)




qsmooth_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
  filter(Sample %in% colData$Sample) %>%
  left_join(sample_meta_D) %>%
  #filter(Gene %in% row.names(top.table %>% head(5))) %>%
  filter(Gene %in% c('ATOH7','NEUROG2','ABCA4')) %>%
  mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
  ggplot(aes(x=Fusion, y=`log2(Counts)`, color = Organism, shape = Technology)) +
  # geom_boxplot(aes(group = Fusion), color = 'Black', outlier.colour = NA) +
  # geom_boxplot(aes(group = interaction(Organism,Fusion)), outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom(size = 3, alpha = 0.7) +
  cowplot::theme_cowplot() +
  facet_wrap(~Gene, scales = 'free_y') +
  ggsci::scale_color_aaas() +
  stat_summary(fun=median, geom = 'line', aes(group = Organism))


# volcano
volcano_maker <- function(df, title="Volcano Plot", pvalue='P.Value', padj='adj.P.Val', logFC='logFC', gene_list = ''){
  df$pvalue <- df[,pvalue]
  df$log2FoldChange <- df[,logFC]
  df$padj <- df[,padj]
  df$Gene <- row.names(df)
  df <- df[!is.na(df$pvalue),]
  print(dim(df))
  df$Class <- 'Not Significant'
  df$Class[df[,'padj'] < 0.05] <- "FDR < 0.05"
  df$GeneT <- df$Gene

  df$Gene[!df$Gene %in% gene_list] <- ''

  plot <- ggplot(data=df,aes(label=Gene, x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(colour=Class)) +
    scale_colour_manual(values=c("darkred","grey")) +
    cowplot::theme_cowplot() +
    geom_vline(aes(xintercept=-1),linetype="dotted") +
    geom_vline(aes(xintercept=1),linetype="dotted") +
    geom_vline(aes(xintercept=-2),linetype="dotted") +
    geom_vline(aes(xintercept=2),linetype="dotted") +
    geom_label_repel(data=bind_rows(subset(df, padj<0.05)),
                     aes(label=GeneT)) +
    xlab('logFC') + ylab('-log10(p value)') +
    ggtitle(title) + cowplot::theme_cowplot()

    plot
}


