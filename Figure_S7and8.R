```{r init}
library(Matrix)
library(tidyverse)
library(biomaRt)
library(monocle3)
library(rgl)
library(tricycle)
library(SignacX)
library(princurve)
library(ggplot2)
library(Hmisc)
library(pheatmap)
library("ComplexHeatmap")
library(circlize)
library(dplyr)
library(ggrepel)
library(latex2exp)
library(scales)
```

```{r loadData}
cds<-readRDS("BRN1-2.rds")
fname <- 'Figure_S7'
```

```{r Pannel E}
genes <- c('Ezh2','Ngn2')
indx<-which(fData(cds)$gene_short_name %in% genes)
cdsP <- cds[indx,pData(cds)$celltype %in% c("AP")]

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 12.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 12.5]))
h12 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")])))/total 
rownames(h12) <- fData(cdsP)$gene_short_name
colnames(h12) <- c('Control')
h12$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/total
h12[is.na(h12)] <- 0
h12 <- (h12-min(h12))/(max(h12)-min(h12))

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 14.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 14.5]))
h14 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")])))/total
rownames(h14) <- fData(cdsP)$gene_short_name
colnames(h14) <- c('Control')
h14$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/total
h14[is.na(h14)] <- 0
h14 <- (h14-min(h14))/(max(h14)-min(h14))

expression12 <- reshape2::melt(as.matrix(h12),varnames = c('genes','condition'),value.name = 'mean')
percent12 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))[2])
rownames(percent12) <- fData(cdsP)$gene_short_name
colnames(percent12) <- c('Control')
percent12$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))[2]
percent12 <- reshape2::melt(as.matrix(percent12),varnames = c('genes','condition'),value.name = 'percent')

expression14 <- reshape2::melt(as.matrix(h14),varnames = c('genes','condition'),value.name = 'mean')
percent14 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))[2])
rownames(percent14) <- fData(cdsP)$gene_short_name
colnames(percent14) <- c('Control')
percent14$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))[2]
percent14 <- reshape2::melt(as.matrix(percent14),varnames = c('genes','condition'),value.name = 'percent')

pdf(paste0("plots/",fname,"_PannelE_AP.pdf"))

  ggplot(reshape::melt(as.matrix(h12)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+     
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E12.5 AP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  
  ggplot(reshape::melt(as.matrix(h14)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+ 
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E14.5 AP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  dot12 <- merge(expression12,percent12)
  ggplot(dot12,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E12.5 AP')
  
  dot14 <- merge(expression14,percent14)
  ggplot(dot14,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E14.5 AP')
  
dev.off()

cdsP <- cds[indx,pData(cds)$celltype %in% c("BP")]

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 12.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 12.5]))
h12 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")])))/total 
rownames(h12) <- fData(cdsP)$gene_short_name
colnames(h12) <- c('Control')
h12$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/total
h12[is.na(h12)] <- 0
h12 <- (h12-min(h12))/(max(h12)-min(h12))

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 14.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 14.5]))
h14 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")])))/total
rownames(h14) <- fData(cdsP)$gene_short_name
colnames(h14) <- c('Control')
h14$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/total
h14[is.na(h14)] <- 0
h14 <- (h14-min(h14))/(max(h14)-min(h14))

expression12 <- reshape2::melt(as.matrix(h12),varnames = c('genes','condition'),value.name = 'mean')
percent12 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))[2])
rownames(percent12) <- fData(cdsP)$gene_short_name
colnames(percent12) <- c('Control')
percent12$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))[2]
percent12 <- reshape2::melt(as.matrix(percent12),varnames = c('genes','condition'),value.name = 'percent')

expression14 <- reshape2::melt(as.matrix(h14),varnames = c('genes','condition'),value.name = 'mean')
percent14 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))[2])
rownames(percent14) <- fData(cdsP)$gene_short_name
colnames(percent14) <- c('Control')
percent14$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))[2]
percent14 <- reshape2::melt(as.matrix(percent14),varnames = c('genes','condition'),value.name = 'percent')

pdf(paste0("plots/",fname,"_PannelE_BP.pdf"))

  ggplot(reshape::melt(as.matrix(h12)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+     
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E12.5 BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  
  ggplot(reshape::melt(as.matrix(h14)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+ 
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E14.5 BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  dot12 <- merge(expression12,percent12)
  ggplot(dot12,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E12.5 BP')
  
  dot14 <- merge(expression14,percent14)
  ggplot(dot14,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E14.5 BP')
  
dev.off()
```

```{r Pannel C}
genes<-c('Ccnd2', 'Ccna2', 'Cdk1', 'Cdk2ap1', 'Cdkn2d', 'Cdkn2c', 'Cdk4', 'Cdk6', 'Cdc45', 'Kif11', 'Gmnn', 
           'Clspn', 'Cdc25a', 'Cdc25b', 'Top2b', 'Rrm1', 'Rrm2', 'Fau', 'Btg2', 'Ppp2r2b')
indx<-which(fData(cds)$gene_short_name %in% genes)
cdsP <- cds[indx,pData(cds)$celltype %in% c("BP")]

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 12.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 12.5]))
h12 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")])))/total 
rownames(h12) <- fData(cdsP)$gene_short_name
colnames(h12) <- c('Control')
h12$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/total
h12[is.na(h12)] <- 0
h12 <- (h12-min(h12))/(max(h12)-min(h12))

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 14.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 14.5]))
h14 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")])))/total
rownames(h14) <- fData(cdsP)$gene_short_name
colnames(h14) <- c('Control')
h14$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/total
h14[is.na(h14)] <- 0
h14 <- (h14-min(h14))/(max(h14)-min(h14))

expression12 <- reshape2::melt(as.matrix(h12),varnames = c('genes','condition'),value.name = 'mean')
percent12 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))[2])
rownames(percent12) <- fData(cdsP)$gene_short_name
colnames(percent12) <- c('Control')
percent12$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))[2]
percent12 <- reshape2::melt(as.matrix(percent12),varnames = c('genes','condition'),value.name = 'percent')

expression14 <- reshape2::melt(as.matrix(h14),varnames = c('genes','condition'),value.name = 'mean')
percent14 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))[2])
rownames(percent14) <- fData(cdsP)$gene_short_name
colnames(percent14) <- c('Control')
percent14$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))[2]
percent14 <- reshape2::melt(as.matrix(percent14),varnames = c('genes','condition'),value.name = 'percent')

pdf(paste0("plots/",fname,"_PannelC_BP.pdf"))

  ggplot(reshape::melt(as.matrix(h12)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+     
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E12.5 BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  
  ggplot(reshape::melt(as.matrix(h14)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+ 
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E14.5 BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  dot12 <- merge(expression12,percent12)
  ggplot(dot12,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E12.5 BP')
  
  dot14 <- merge(expression14,percent14)
  ggplot(dot14,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E14.5 BP')
  
dev.off()

cdsP <- cds[indx,pData(cds)$celltype %in% c("AP")]

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 12.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 12.5]))
h12 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.control")])))/total 
rownames(h12) <- fData(cdsP)$gene_short_name
colnames(h12) <- c('Control')
h12$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("12.5.ko")]))/total
h12[is.na(h12)] <- 0
h12 <- (h12-min(h12))/(max(h12)-min(h12))

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 14.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 14.5]))
h14 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")])))/total
rownames(h14) <- fData(cdsP)$gene_short_name
colnames(h14) <- c('Control')
h14$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/total
h14[is.na(h14)] <- 0
h14 <- (h14-min(h14))/(max(h14)-min(h14))

expression12 <- reshape2::melt(as.matrix(h12),varnames = c('genes','condition'),value.name = 'mean')
percent12 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.control')]))[2])
rownames(percent12) <- fData(cdsP)$gene_short_name
colnames(percent12) <- c('Control')
percent12$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('12.5.ko')]))[2]
percent12 <- reshape2::melt(as.matrix(percent12),varnames = c('genes','condition'),value.name = 'percent')

expression14 <- reshape2::melt(as.matrix(h14),varnames = c('genes','condition'),value.name = 'mean')
percent14 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))[2])
rownames(percent14) <- fData(cdsP)$gene_short_name
colnames(percent14) <- c('Control')
percent14$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))[2]
percent14 <- reshape2::melt(as.matrix(percent14),varnames = c('genes','condition'),value.name = 'percent')

pdf(paste0("plots/",fname,"_PannelC_AP.pdf"))

  ggplot(reshape::melt(as.matrix(h12)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+     
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E12.5 AP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  
  ggplot(reshape::melt(as.matrix(h14)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+ 
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E14.5 AP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  dot12 <- merge(expression12,percent12)
  ggplot(dot12,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E12.5 AP')
  
  dot14 <- merge(expression14,percent14)
  ggplot(dot14,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E14.5 AP')
  
dev.off()
```

```{r Pannel 7D and 8F}
cdsP <- cds[,pData(cds)$celltype %in% c('AP')]

gene_fits <- fit_models(cdsP[,pData(cdsP)$age==14.5], model_formula_str = "~Batch+condition",cores = 10)
fit_coefs <- coefficient_table(gene_fits)
unique(coefficient_table(gene_fits)$term)

sig_terms <- DEcondition %>% filter (q_value < 0.001 & std_err > 0)
sig_genes <- sig_terms %>% pull(gene_short_name)

s2 <- sig_terms
s2 <- s2[s2$std_err < 0.15,]
s2$label <- s2$gene_short_name
s2$q_value[s2$q_value == 0] <- 1E-310

s2 <- s2[!grepl('Rpl',s2$gene_short_name),]
s2 <- s2[!grepl('Rps',s2$gene_short_name),]
s2 <- s2[!grepl('mt-',s2$gene_short_name),]

highlight <- c('') #fill in the gene names
pdf(paste0("plots/",fname,"_PannelD_AP.pdf"))
  ggplot(s2,aes(x=normalized_effect,y=-log10(q_value),label=label)) +
    geom_point(data = s2[!s2$gene_short_name %in%(highlight),],aes(x=normalized_effect,y=-log10(q_value),label=label), color = "black")+
    geom_point(data = s2[s2$gene_short_name %in%c(highlight),],aes(x=normalized_effect,y=-log10(q_value),label=label), color = "red")+
    geom_text_repel(data =s2[s2$gene_short_name %in%c(highlight),], max.overlaps = Inf,direction = 'y') +
    xlab(paste("Normalized effect size conditionko")) +
    geom_vline(xintercept = 0,linetype="dashed") +
    ylab(TeX(r"($-log_{10}$ q-value)")) +
    theme(legend.position = "none") + 
    monocle3:::monocle_theme_opts()
dev.off()

cdsP <- cds[,pData(cds)$celltype %in% c('BP')]

gene_fits <- fit_models(cdsP[,pData(cdsP)$age==14.5], model_formula_str = "~Batch+condition",cores = 10)
fit_coefs <- coefficient_table(gene_fits)
unique(coefficient_table(gene_fits)$term)

sig_terms <- DEcondition %>% filter (q_value < 0.001 & std_err > 0)
sig_genes <- sig_terms %>% pull(gene_short_name)

s2 <- sig_terms
s2 <- s2[s2$std_err < 0.15,]
s2$label <- s2$gene_short_name
s2$q_value[s2$q_value == 0] <- 1E-310

s2 <- s2[!grepl('Rpl',s2$gene_short_name),]
s2 <- s2[!grepl('Rps',s2$gene_short_name),]
s2 <- s2[!grepl('mt-',s2$gene_short_name),]

highlight <- c(' ') #fill in the gene names
pdf(paste0("plots/",fname,"_PannelD_BP.pdf"))
  ggplot(s2,aes(x=normalized_effect,y=-log10(q_value),label=label)) +
    geom_point(data = s2[!s2$gene_short_name %in%(highlight),],aes(x=normalized_effect,y=-log10(q_value),label=label), color = "black")+
    geom_point(data = s2[s2$gene_short_name %in%c(highlight),],aes(x=normalized_effect,y=-log10(q_value),label=label), color = "red")+
    geom_text_repel(data =s2[s2$gene_short_name %in%c(highlight),], max.overlaps = Inf,direction = 'y') +
    xlab(paste("Normalized effect size conditionko")) +
    geom_vline(xintercept = 0,linetype="dashed") +
    ylab(TeX(r"($-log_{10}$ q-value)")) +
    theme(legend.position = "none") + 
    monocle3:::monocle_theme_opts()
dev.off()
```

```{r Pannel F2}
genes <- c('Robo2', 'Rorb', 'Cux2', 'Ptn', 'Sema6a', 'Rprm', 'Fezf2', 'Foxp2')
indx<-which(fData(cds)$gene_short_name %in% genes)
cdsP <- cds[indx,pData(cds)$celltype %in% c("BP")]

total <- rowSums(exprs(cdsP[,pData(cdsP)$age == 14.5]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age == 14.5]))
h14 <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.control")])))/total
rownames(h14) <- fData(cdsP)$gene_short_name
colnames(h14) <- c('Control')
h14$dKO <-rowSums(exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c("14.5.ko")]))/total
h14[is.na(h14)] <- 0
h14 <- (h14-min(h14))/(max(h14)-min(h14))

expression14 <- reshape2::melt(as.matrix(h14),varnames = c('genes','condition'),value.name = 'mean')
percent14 <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.control')]))[2])
rownames(percent14) <- fData(cdsP)$gene_short_name
colnames(percent14) <- c('Control')
percent14$dKO <- rowSums(!!exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))/dim(exprs(cdsP[,pData(cdsP)$age.cond %in% c('14.5.ko')]))[2]
percent14 <- reshape2::melt(as.matrix(percent14),varnames = c('genes','condition'),value.name = 'percent')

pdf(paste0("plots/FigureS8_PannelF2.pdf"))

  ggplot(reshape::melt(as.matrix(h14)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+ 
    scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
    scale_y_discrete(limits = rev(genes)) + 
    labs(title = 'E14.5 BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.background = element_blank())
  dot14 <- merge(expression14,percent14)
  ggplot(dot14,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E14.5 BP')
  
dev.off()
```
