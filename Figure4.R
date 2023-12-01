```{r init}
library(Matrix)
library(monocle3)
library(tricycle)
library(dplyr)
library(RColorBrewer)
library(scales)
library(ggplot2)
```

```{r load data}
cds <- readRDS('E36_v2.rds')
fname <- 'Figure4'

S_matrix <- SingleCellExperiment::reducedDims(cds)[['UMAP']]
data_df <- data.frame(S_matrix)

colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- row.names(data_df)

data_df <- as.data.frame(cbind(data_df, colData(cds)))
```{r}

```{r Pannel J condition}
pdf(paste0("plots/",fname,"_PannelJ.pdf"),width = 12,height = 10)
ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =condition), size=1)+
  guides(color = guide_legend(title = "Condition",override.aes = list(size = 4)))+
  scale_colour_manual(values=c('WT'='#FDC086','KO'='#386CB0'))+
  ggtitle('Condition') + 
  xlab('UMAP1') +
  ylab('UMAP2') +
  monocle3:::monocle_theme_opts()
dev.off()
```

```{r Pannel K celltype}
pdf(paste0("plots/",fname,"_PannelK.pdf"),width=12,height=10) 
ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =celltype), size=0.5)+
  guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))+
  ggtitle('Cell Type') + 
  xlab('UMAP1') +
  ylab('UMAP2') +
  monocle3:::monocle_theme_opts()
dev.off()
```

```{r Pannel L}
pdf(paste0("plots/",fname,"_PannelL_signature.pdf"),width = 12,height = 10)
  ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =BP_signature), size=1)+
    #guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))+
    scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits =c(-2,10),oob=squish)+
    ggtitle('BP signature') + 
    xlab('UMAP1') +
    ylab('UMAP2') +
    monocle3:::monocle_theme_opts()

  ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =AP_signature), size=1)+
    #guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))+
    scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits =c(-2,10),oob=squish)+
    ggtitle('AP signature') + 
    xlab('UMAP1') +
    ylab('UMAP2') +
    monocle3:::monocle_theme_opts()

  ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =Neuron_signature), size=1)+
    #guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))+
    scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits =c(NA,12),oob=squish)+
    ggtitle('Neuron signature') + 
    xlab('UMAP1') +
    ylab('UMAP2') +
    monocle3:::monocle_theme_opts()
dev.off()

genes <- c('NES','BTG2','DCX')
pdf(paste0("plots/",fname,"_PannelL_marker.pdf"),width = 12,height = 10)
  for (i in genes){
  print(monocle3::plot_cells(cds,show_trajectory_graph = FALSE,label_cell_groups = FALSE,
                       genes = i,min_expr = 0,scale_to_range = FALSE,cell_size = 1) + 
        scale_color_gradientn(colors = brewer.pal(9,"YlGnBu"),na.value = "grey90") + 
        ggtitle(paste('E36',i)))
    }
dev.off()
```

```{r Pannel M}
pdf(paste0("plots/",fname,"_PannelM_BPsignature.pdf"),width = 12,height = 10)
  ggplot(data=data_df, aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =BP_signature), size=1)+
    #guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))+
    scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits =c(-2,10),oob=squish)+
    ggtitle('BP signature') + 
    xlab('UMAP1') +
    ylab('UMAP2') +
    monocle3:::monocle_theme_opts()

  for (i in unique(pData(cds)$condition))
  {
    g <- ggplot(data=data_df[data_df$condition == i,], aes(x=data_dim_1, y=data_dim_2)) + geom_point(aes(color =BP_signature), size=1)+
      #guides(color = guide_legend(title = "Cell Type",override.aes = list(size = 4)))+
      scale_color_gradientn(colors = rev(brewer.pal(10,"Spectral")),limits =c(-2,10),oob=squish)+
      ggtitle(paste(i,'BP signature')) + 
      xlab('UMAP1') +
      ylab('UMAP2') +
      monocle3:::monocle_theme_opts()
    print(g)
  }
dev.off()

pdf(paste0("plots/",fname,"_PannelM_MKI67.pdf"),width = 12,height = 10)
  monocle3::plot_cells(cds,show_trajectory_graph = FALSE,label_cell_groups = FALSE,
                       genes = 'EOMES',min_expr = 0,scale_to_range = FALSE,cell_size = 1) + 
    scale_color_gradientn(colors = brewer.pal(9,"YlGnBu"),na.value = "grey90") + 
    ggtitle('E36')

  for (i in unique(pData(cds)$condition))
  {
    g <- monocle3::plot_cells(cds[,pData(cds)$condition ==i],show_trajectory_graph = FALSE,label_cell_groups = FALSE,
                              genes = 'EOMES',min_expr = 0,scale_to_range = FALSE,cell_size = 1) + 
      scale_color_gradientn(colors = brewer.pal(9,"YlGnBu"),na.value = "grey90")+
      ggtitle(i)
    print(g)
  }
dev.off()

hue.colors = c("#2E22EA","#9E3DFB", "#F86BE2", "#FCCE7B", "#C4E416", "#4BBA0F", "#447D87", "#2C24E9")
hue.n = 50

pdf(paste0("plots/",fname,"PannelM_tricycle.pdf"),width=12,height=10)
  plot_cells(cds,color_cells_by = 'tricyclePosition',label_cell_groups = FALSE)+
    scale_color_gradientn(limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = "none")
  for (i in unique(pData(cds)$condition))
  {
    g <- plot_cells(cds,color_cells_by = 'tricyclePosition',label_cell_groups = FALSE)+
    scale_color_gradientn(limits = range(0, 2 * pi), breaks = seq(from = 0, to = 2 * pi, length.out = hue.n), colors = hue.colors, guide = "none")+
      ggtitle(i)
    print(g)
  }
dev.off()
```

```{r Pannel N and O}
pdf(paste0("plots/",fname,"PannelN.pdf"),width=12,height=10)
  plot_cells(cds,color_cells_by="in_dir",label_cell_groups=FALSE,cell_size = 0.5) + 
    ggtitle("BP type") + guides(color = guide_legend(title = "BP type",override.aes = list(size = 4)))
dev.off()

CTfeq<-tibble::as_tibble(pData(cds[,pData(cds)$in_dir%in%c('direct','indirect')])) %>%
  group_by(condition,in_dir) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))


pdf(paste0("plots/",fname,"PannelO.pdf"),width = 12,height = 10)
  ggplot(data = CTfeq,aes(x = freq, y = condition, fill = in_dir))+
    geom_bar(position="fill",stat="identity") + 
    scale_y_discrete(limits=rev)+
    labs(x = 'cell proportion',y = '')+ 
    ggtitle('E36') + 
    theme(panel.grid.major = element_blank(),strip.background = element_blank(), 
          strip.text.x = element_text(size = 18),
          panel.grid.minor = element_blank(),panel.background = element_blank())
dev.off()
```

```{r Pannel P}
genes <- c('HES1','NOTCH1','DLL1')
indx<-which(fData(cds)$gene_short_name %in% genes)
cdsP <- cds[indx,pData(cds)$celltype %in% c("AP","BP") ]

total <- rowSums(exprs(cdsP))/rowSums(!!exprs(cdsP))
h <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')]))/(rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')])))/total) 
rownames(h) <- fData(cdsP)$gene_short_name
colnames(h) <- c('Control')
h$KO <-rowSums(exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/total
h[is.na(h)] <- 0
h <- (h-min(h))/(max(h)-min(h))

pdf(paste0("plots/",fname,"PannelP.pdf"),width=10,height=10)
  ggplot(reshape2::melt(as.matrix(h)), aes(x = Var2, y =Var1 , fill = value))+ geom_tile(color = "white")+     
    #scale_fill_gradientn(colors=hcl.colors(16, palette = "GnBu", alpha = NULL, rev = FALSE, fixup = TRUE),guide="colorbar")+ coord_fixed(ratio = 0.8) +
    scale_fill_gradientn(colors=viridis(16),guide="colorbar")+ coord_fixed(ratio = 0.8) +
    scale_y_discrete(limits = genes) + 
    labs(title = 'E36 AP BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                                 panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank())

  expression <- reshape2::melt(as.matrix(h),varnames = c('genes','condition'),value.name = 'mean')
  percent <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')]))/dim(cdsP[,pData(cdsP)$sample %in% c('E36_WT')])[2])
  rownames(percent) <- fData(cdsP)$gene_short_name
  colnames(percent) <- c('Control')
  percent$KO <- rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/dim(cdsP[,pData(cdsP)$sample %in% c('E36_KO')])[2]
  percent <- reshape2::melt(as.matrix(percent),varnames = c('genes','condition'),value.name = 'percent')

  dot <- merge(expression,percent)
  ggplot(dot,aes(x=condition, y = genes, color = mean, size = percent)) + 
          geom_point()  + scale_y_discrete(limits = genes) +
          scale_color_gradientn(colors = viridis(16),name = 'mean expression')+
          scale_size(name = 'proportion per condition', range = c(0, 10),limits = c(0,NA))+ monocle3:::monocle_theme_opts() +
          labs(title ='E36 AP BP')
dev.off()
```

```{r Pannel Q}
genes <- c('RPRM', 'LHX2', 'ISLR2', 'CUX2', 'RORB', 'ROBO2')
indx<-which(fData(cds)$gene_short_name %in% genes)
cdsP <- cds[indx,pData(cds)$celltype %in% c("Neuron") ]

total <- rowSums(exprs(cdsP))/rowSums(!!exprs(cdsP))
h <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')]))/(rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')])))/total) 
rownames(h) <- fData(cdsP)$gene_short_name
colnames(h) <- c('Control')
h$KO <-rowSums(exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/total
h[is.na(h)] <- 0
h <- (h-min(h))/(max(h)-min(h))

pdf(paste0("plots/",fname,"PannelM.pdf"),width=10,height=10)
  ggplot(reshape2::melt(as.matrix(h)), aes(x = Var2, y =Var1 , fill = value))+ geom_tile(color = "white")+     
    #scale_fill_gradientn(colors=hcl.colors(16, palette = "GnBu", alpha = NULL, rev = FALSE, fixup = TRUE),guide="colorbar")+ coord_fixed(ratio = 0.8) +
    scale_fill_gradientn(colors=viridis(16),guide="colorbar")+ coord_fixed(ratio = 0.8) +
    scale_y_discrete(limits = genes) + 
    labs(title = 'E36 Neuron', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                                 panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank())

  expression <- reshape2::melt(as.matrix(h),varnames = c('genes','condition'),value.name = 'mean')
  percent <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')]))/dim(cdsP[,pData(cdsP)$sample %in% c('E36_WT')])[2])
  rownames(percent) <- fData(cdsP)$gene_short_name
  colnames(percent) <- c('Control')
  percent$KO <- rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/dim(cdsP[,pData(cdsP)$sample %in% c('E36_KO')])[2]
  percent <- reshape2::melt(as.matrix(percent),varnames = c('genes','condition'),value.name = 'percent')

  dot <- merge(expression,percent)
  ggplot(dot,aes(x=condition, y = genes, color = mean, size = percent)) + 
          geom_point()  + scale_y_discrete(limits = genes) +
          scale_color_gradientn(colors = viridis(16),name = 'mean expression')+
          scale_size(name = 'proportion per condition', range = c(0, 10),limits = c(0,NA))+ monocle3:::monocle_theme_opts() +
          labs(title ='E36 Neuron')
dev.off()
```

```{r Pannel L}
genes <- c('MCPH1','CDK5RAP2','ASPM','CENPJ','STIL','CEP135','SASS6','CIT','COPB2')
indx<-which(fData(cds)$gene_short_name %in% genes)
cdsP <- cds[indx,pData(cds)$celltype %in% c("AP","BP") ]

total <- rowSums(exprs(cdsP))/rowSums(!!exprs(cdsP))
h <- data.frame(rowSums(exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')]))/(rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')])))/total) 
rownames(h) <- fData(cdsP)$gene_short_name
colnames(h) <- c('Control')
h$KO <-rowSums(exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/total
h[is.na(h)] <- 0
h <- (h-min(h))/(max(h)-min(h))

pdf(paste0("plots/",fname,"PannelL.pdf"),width=10,height=10)
  ggplot(reshape2::melt(as.matrix(h)), aes(x = Var2, y =Var1 , fill = value))+ geom_tile(color = "white")+     
    #scale_fill_gradientn(colors=hcl.colors(16, palette = "GnBu", alpha = NULL, rev = FALSE, fixup = TRUE),guide="colorbar")+ coord_fixed(ratio = 0.8) +
    scale_fill_gradientn(colors=viridis(16),guide="colorbar")+ coord_fixed(ratio = 0.8) +
    scale_y_discrete(limits = genes) + 
    labs(title = 'E36 AP BP', x = 'Condition',y = 'Genes')+ theme(panel.grid.major = element_blank(),   
                                                                                 panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank())

  expression <- reshape2::melt(as.matrix(h),varnames = c('genes','condition'),value.name = 'mean')
  percent <- data.frame(rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_WT')]))/dim(cdsP[,pData(cdsP)$sample %in% c('E36_WT')])[2])
  rownames(percent) <- fData(cdsP)$gene_short_name
  colnames(percent) <- c('Control')
  percent$KO <- rowSums(!!exprs(cdsP[,pData(cdsP)$sample %in% c('E36_KO')]))/dim(cdsP[,pData(cdsP)$sample %in% c('E36_KO')])[2]
  percent <- reshape2::melt(as.matrix(percent),varnames = c('genes','condition'),value.name = 'percent')

  dot <- merge(expression,percent)
  ggplot(dot,aes(x=condition, y = genes, color = mean, size = percent)) + 
          geom_point()  + scale_y_discrete(limits = genes) +
          scale_color_gradientn(colors = viridis(16),name = 'mean expression')+
          scale_size(name = 'proportion per condition', range = c(0, 10),limits = c(0,NA))+ monocle3:::monocle_theme_opts() +
          labs(title ='E36 AP BP')
dev.off()
```
