library(edgeR)
ngs_counts <- qnorm_counts
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

cutoff <- 5
drop <- which(apply((ngs_counts), 1, max) < cutoff) #which(apply(cpm(d_zed), 1, max) < cutoff)
d <- ngs_counts[-drop,]
dim(d) # number of genes left


# reduce down to OF and OFM (optic fissure (margin))
ofm_meta <- sample_meta_D %>% filter(Section %in% c('OF','OFM'))
d_ofm <- d[,ofm_meta$Sample]
colData <- colnames(d_ofm) %>% as_tibble() %>% dplyr::rename(Sample = value) %>%  left_join(ofm_meta)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample


mm <- model.matrix(~0  + colData$Fusion + colData$Organism + colData$Technology)
colnames(mm) <- c('After','Before','During', 'Mouse','Zebrafish', 'RNAseq')
#y <- voom(d, mm, plot = T)
fit <- lmFit(d_ofm, mm)
contrast.matrix = makeContrasts(After-During, After-Before, During-Before, levels=mm)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)



efit <- eBayes(fit_contrasts)
top.table <- topTable(efit, sort.by = "p", n = Inf, coef="After - During")
head(top.table, 20)

top.table <- topTable(efit, sort.by = "p", n = Inf, coef="During - Before")
head(top.table, 20)

top.table <- topTable(efit, sort.by = "p", n = Inf, coef="After - Before")
head(top.table, 20)




qnorm_counts %>%
  as_tibble(rownames = 'Gene') %>%
  pivot_longer(-Gene, names_to = 'Sample', values_to = 'log2(Counts)') %>%
  filter(Sample %in% colData$Sample) %>%
  left_join(sample_meta_D) %>%
  filter(Gene %in% row.names(top.table %>% head(5))) %>%
  #filter(Gene %in% c('MITF','NEUROG2')) %>%
  mutate(Fusion = factor(Fusion, levels = c('Before','During','After'))) %>%
  ggplot(aes(x=Fusion, y=`log2(Counts)`, color = Organism, shape = Technology)) +
  # geom_boxplot(aes(group = Fusion), color = 'Black', outlier.colour = NA) +
  geom_boxplot(aes(group = interaction(Organism,Fusion)), outlier.colour = NA) +
  ggbeeswarm::geom_quasirandom(size = 3, alpha = 0.7) +
  cowplot::theme_cowplot() +
  facet_wrap(~Gene, scales = 'free_y') +
  ggsci::scale_color_aaas()
