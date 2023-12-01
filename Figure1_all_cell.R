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
fname<-"Figure1"
```

```{r Pannel E}
pdf(paste0("plots/",fname,"PannelE_signature.pdf"),width=12,height=10)
 plot_cells(cds,color_cells_by="AP_signature",label_cell_groups=FALSE) + 
  scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits = c(NA,20)) + ggtitle("AP signature")+
  guides(color = guide_colorbar(title = "AP signature"))

 plot_cells(cds,color_cells_by="BP_signature",label_cell_groups=FALSE) + 
  scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits = c(NA,30)) + ggtitle("BP signature")+
  guides(color = guide_colorbar(title = "BP signature"))

 plot_cells(cds,color_cells_by="Neuron_signature",label_cell_groups=FALSE) + 
  scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits = c(-5,15),oob = squish) + ggtitle("Neuron signature")+
  guides(color = guide_colorbar(title = "Neuron signature"))
dev.off()
```

```{r Pannel F}
pdf(paste0("plots/",fname,"PannelF_celltype.pdf"),width=12,height=10)
    plot_cells(cds,color_cells_by="celltype",label_cell_groups=FALSE) + 
       ggtitle("Cell Type") + guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))
dev.off()
```

```{r Pannel J}
hue.colors = c("#2E22EA","#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
hue.n = 50

pdf(paste0("plots/",fname,"PannelJ_tricycle.pdf"),width=12,height=10)
  plot_cells(cds,color_cells_by = 'tricyclePosition',label_cell_groups = FALSE)+
    scale_color_gradientn(limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = "none")
dev.off()
```

```{r Pannel K}
cdsP <- cds[,pData(cds)$celltype %in% c('AP')]
pdf(paste0("plots/",fname,"PannelK_AP.pdf"),width=12,height=10)
 plot_ccposition_den(pData(cdsP)$tricyclePosition,pData(cdsP)$age.cond,'age.cond',bw = 10) + ggpubr::theme_pubr()
dev.off()

cdsP <- cds[,pData(cds)$celltype %in% c('BP')]
pdf(paste0("plots/",fname,"PannelK_BP.pdf"),width=12,height=10)
 plot_ccposition_den(pData(cdsP)$tricyclePosition,pData(cdsP)$age.cond,'age.cond',bw = 10) + ggpubr::theme_pubr()
dev.off()
```

