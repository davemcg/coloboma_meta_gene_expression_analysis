library(edgeR)
ngs_counts <- same

d_zed <- DGEList(ngs_counts)
d_zed <- calcNormFactors(d_zed)

sample_meta_HomoMus <- sample_meta %>% filter(Sample %in% colnames(ngs_counts)) %>%
  dplyr::select(Sample:Section, Layout:Fusion) %>%
  unique()

colData <- d_zed$samples %>% as_tibble(rownames = 'Sample') %>% left_join(sample_meta_HomoMus)
colData <- data.frame(colData)
row.names(colData) <- colData$Sample

d_zed$samples <- colData

cutoff <- 2
drop <- which(apply(cpm(d_zed), 1, max) < cutoff)
d <- d_zed[-drop,]
dim(d) # number of genes left


mm <- model.matrix(~0 + d$samples$Fusion)
y <- voom(d, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
