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
```

```{r loadData}
cds<-readRDS("BRN1-2.rds")
fname <- 'Figure2'
```

```{r Pannel B}
#Eomes is the alias of Tbr2 in the genome assembly used for alignment
markers <- c('Eomes','Neurod2','Mki67')
pdf(paste0("plots/",fname,"PannelD_marker.pdf"),width=12,height=10)
 for (i in markers){
   print(plot_cells(cds,genes = i,min_expr = 0)+ 
           scale_color_gradientn(colors = brewer.pal(9,"YlGnBu"),na.value = "grey90"))
 }
dev.off()
```

```{r Pannel C}
indirBPs<-tibble::as_tibble(colData(cds)) %>%
  filter(pData(cds)$tricyclePosition < 1.8*pi & pData(cds)$tricyclePosition > .4*pi & pData(cds)$celltype %in% c('BP')|clusters(cds) == 7)

pData(cds)$BP_type <- "APorNeuron"
pData(cds)[pData(cds)$celltype %in% c('BP'),]$BP_type <- 'direct'
pData(cds)[indirBPs$cell_barcode,]$BP_type <- 'indirect'

pdf(paste0("plots/",fname,"PannelC_BPtype.pdf"),width=12,height=10)
   plot_cells(cds,color_cells_by="BP_type",label_cell_groups=FALSE) + ggtitle("BP_type")+
   	scale_colour_manual(values=c('indirect'='#F8766D','direct'='#00BFC4','APorNeuron'='#E5E5E5'))
dev.off()

CTfeq<-tibble::as_tibble(pData(cds)[pData(cds)$BP_type!='APorNeuron',]) %>%
  group_by(condition,age,BP_type) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

pdf(paste0("plots/",fname,"PannelC_BPtype_proportion.pdf"),width = 10,height = 10)
  ggplot(data = CTfeq,aes(x = freq, y = condition, fill = BP_type))+
         geom_bar(position="fill",stat="identity") + 
         scale_y_discrete(limits=rev,labels = c('dKO','Control'))+
         scale_fill_manual(values=c('indirect'='#F8766D','direct'='#00BFC4'))+
         facet_wrap(~age,labeller = as_labeller(c(`12.5` = "E12.5",`14.5` = "E14.5")))+
         labs(x = 'cell proportion',y = '')+ 
         theme(panel.grid.major = element_blank(),strip.background = element_blank(), 
               strip.text.x = element_text(size = 18),
               panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()
```

```{r Pannel F}
genes <- c('Hes1','Notch1','Dll1')
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

pdf(paste0("plots/",fname,"_PannelF.pdf"))
  dot14 <- merge(expression14,percent14)

  ggplot(reshape::melt(as.matrix(h14)), aes(x = X2, y =X1 , fill = value))+ geom_tile(color = "white")+     
      scale_fill_gradientn(colors=viridis(16),guide="colorbar",limits = c(0,1))+ coord_fixed(ratio = 0.8) + 
      scale_y_discrete(limits = rev(genes)) + 
      labs(title = 'E14.5 BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                   panel.grid.minor = element_blank(),
                                                                   panel.background = element_blank())

  ggplot(dot14,aes(x=condition, y = genes, color = mean, size = percent)) + 
    geom_point()  + scale_y_discrete(limits = rev(genes)) +
    scale_color_gradientn(colors = viridis(16),name = 'mean expression',limits = c(0,1))+
    scale_size(name = 'percentage per condition', range = c(0, 10))+ monocle3:::monocle_theme_opts() +
    labs(title = 'E14.5 BP')
dev.off()
```
