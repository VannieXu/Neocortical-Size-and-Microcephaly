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
fname <- 'FigureS5'
```

```{r Pannel B}
#Eomes is the alias of Tbr2 in the genome assembly used for alignment
markers <- c('Nes','Tbr2','Neurod2','Sox2','Btg2','Tbr1')
pdf(paste0("plots/",fname,"PannelB_marker.pdf"),width=12,height=10)
 for (i in markers){
   print(plot_cells(cds,genes = i,min_expr = 0)+ 
           scale_color_gradientn(colors = brewer.pal(9,"YlGnBu"),na.value = "grey90"))
 }
dev.off()
```

```{r Pannel C}
pdf(paste0("plots/",fname,"_PanelC_age_condition.pdf"),width=12,height=10)
  plot_cells(cds,color_cells_by = 'age',label_cell_groups = FALSE)+
    scale_colour_manual(values = c("#66C2A5","#8DA0CB")) + 
    ggtitle("Age E12.5 vs E14.5") 
  
  plot_cells(cds,color_cells_by = 'condition',label_cell_groups = FALSE)+
    scale_colour_manual(values = c("#FDC086", "#386CB0")) + 
    ggtitle("Condition CT vs dKO") 
dev.off()
```

```{r Pannel D}
