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
fname<-"Figure2"
```

```{r Pannel A}
hue.colors = c("#2E22EA","#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
hue.n = 50

pdf(paste0("plots/",fname,"PannelA_tricycle.pdf"),width=12,height=10)
  plot_cells(cds,color_cells_by = 'tricyclePosition',label_cell_groups = FALSE)+
    scale_color_gradientn(limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = "none")
dev.off()
```

```{r Pannel B}
cdsP <- cds[,pData(cds)$celltype %in% c('AP')]
pdf(paste0("plots/",fname,"PannelB_AP.pdf"),width=12,height=10)
 plot_ccposition_den(pData(cdsP)$tricyclePosition,pData(cdsP)$age.cond,'age.cond',bw = 10) + ggpubr::theme_pubr()
dev.off()

cdsP <- cds[,pData(cds)$celltype %in% c('BP')]
pdf(paste0("plots/",fname,"PannelB_BP.pdf"),width=12,height=10)
 plot_ccposition_den(pData(cdsP)$tricyclePosition,pData(cdsP)$age.cond,'age.cond',bw = 10) + ggpubr::theme_pubr()
dev.off()
```

```{r Pannel I}
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

pdf(paste0("plots/",fname,"_PannelI_AP.pdf"))

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

pdf(paste0("plots/",fname,"_PannelI_BP.pdf"))

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
