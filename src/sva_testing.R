

num_sv <- num.sv(d_filtered,mm,method="leek")
mod0 = model.matrix(~1, data=colData)
sv_obj <- sva(as.matrix(d_filtered), mm, mod0,n.sv=num_sv)

cor_vals <- removeBatchEffect(d_filtered, covariates = sv_obj$sv)


modSv = cbind(mm,sv_obj$sv)
colnames(modSv)[7:8] <- c('SVA1','SVA2')


#y <- voom(d, mm, plot = T)
fit <- lmFit(d_filtered, modSv)
contrast.matrix = makeContrasts(During-Before, After-During,levels=modSv)
fit_contrasts <- contrasts.fit(fit, contrast.matrix)



efit <- eBayes(fit_contrasts)
top.table <- topTable(efit, sort.by = "p", n = Inf, coef="After - During", adjust.method="BH")
head(top.table, 20)
top.table %>% filter(adj.P.Val<0.05) %>% dim()


top.table <- topTable(efit, sort.by = "p", n = Inf, coef="During - Before")
head(top.table, 20)
top.table %>% filter(adj.P.Val<0.05) %>% dim()
