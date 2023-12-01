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

```{r pannel 2F and 3A}
cdsP <- cds[,pData(cds)$celltype %in% c('AP','BP')]

gene_fits <- fit_models(cdsN[,pData(cdsN)$age==14.5], model_formula_str = "~Batch+condition",cores = 10)
fit_coefs <- coefficient_table(gene_fits)
unique(coefficient_table(gene_fits)$term)

sig_terms <- DEcondition %>% filter (q_value < 0.05 & std_err > 0)
sig_genes <- sig_terms %>% pull(gene_short_name)

s2 <- sig_terms
s2 <- s2[s2$std_err < 0.15,]
s2$label <- s2$gene_short_name
s2$q_value[s2$q_value == 0] <- 1E-310

s2 <- s2[!grepl('Rpl',s2$gene_short_name),]
s2 <- s2[!grepl('Rps',s2$gene_short_name),]
s2 <- s2[!grepl('mt-',s2$gene_short_name),]

fname <- 'Figure2'
highlight <- c('Notch1','Dll1','Hes1')
pdf(paste0("plots/",fname,"_PannelF.pdf"))
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

fname <- 'Figure3'
highlight <- c('Tacc1','Afdn','Pak3','Smarca5','Zbtb18','Kif11','Msmo1','Tpx2','Haus1','Bub1b','Ckap2l','Plk4','Ccnd2','Auts2','Wls','Nav2','Ier3ip1',
               'Col4a1','Sgce','Mid1','Phgdh','Dab1','Col4a2','Zeb2','Celsr1','Cdk5rap2','Knl1','Aspm','Cenpj','Stil','Cep135','Cdk6','Sass6',
               'Mfsd2a','Cit','Copb2')
pdf(paste0("plots/",fname,"_PannelA.pdf"))
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

