---
title: "Differential Expression Analysis of scRNAseq Data"
author: "Christopher O'Connor"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
# Chunk Timer capture
library(tictoc)
source("E:/timer_hook.R")
```

## Library

```{r message=FALSE, warning=FALSE}
# Load packages
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(scater)
library(ggthemes)
library(corrplot)
library(Polychrome)
library(ggbeeswarm)
library(ggridges)
library(igraph)

library(Seurat) # scRNAseq analysis
library(SeuratWrappers)
library(SeuratDisk)
library(SingleR) # Used for annotating cell type clusters
library(monocle3) # Trajectory analysis
library(SingleCellExperiment)
library(celldex) # Supports cell annotation moving forward as SingleR's database becomes depreciated 
```

## Read in Data

```{r}
# Read in data object from previous script (02_Doublet.R)
setwd("E:")
sobj <- LoadH5Seurat("../new_02_k23_object.h5Seurat")
```

```{r}
imm <- sobj
head(imm@meta.data)

# Remove excess object
rm(sobj)
```

## Summary Plots

```{r}
# Visualise cluster markings
DimPlot(imm, reduction = "umap", group.by = "sample") # grouped by samples
DimPlot(imm, reduction = "umap", split.by = "sample", label = TRUE)
DimPlot(imm, reduction = "umap", label = TRUE)

# Target Genes
# FeaturePlot(imm, features = c("Ncr1", "Cd3d", "Cd8a", "Cd4", "Tcf7", "Ctla4", "Lag3", "Pdcd1", "Foxp3", "Icos", "Lef1", "Gzmb"))

```

# Differential Expression Analysis

# 1. Marker Identification

## 1.2 Find All Markers

```{r FindAllMarkers, start_timer=TRUE, name = "FindAllMarkers_chunk"}
# Identify marker genes between all clusters
imm.markers <- FindAllMarkers(object = imm, 
                               only.pos = TRUE, # sig. upregulated 
                               min.pct = 0.25,  # min. expression (25%) 
                               thresh.use = 0.25,
                               assay = "RNA")
# marker must be expressed atleast 25% in a cluster compared to others
imm.markers
```

```{r}
# Top 'X' markers where 'n' is x
imm.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> top3

# imm <- ScaleData(imm, 
#                     features = as.character(unique(top3$gene)), 
#                     assay = "RNA")
```

```{r}
# Visualise top three markers for each clusters against each other
DoHeatmap(imm, 
          features = as.character(unique(top3$gene)),
          group.by = "RNA_snn_res.0.5",
          assay = "RNA")

DotPlot(imm,
        features = rev(as.character(unique(top3$gene))), 
        group.by = "RNA_snn_res.0.5",   
        assay = "RNA") + 
        coord_flip()

DotPlot(imm,
        features = rev(as.character(unique(top3$gene))), 
        group.by = "sample",   
        assay = "RNA") + 
        coord_flip()
```

## 1.3 Preliminary Cell Type Investigation

```{r message=FALSE, warning=FALSE}
feFun <- function(imm, feature, title) {
      plot <- FeaturePlot(imm, features = feature) +
              scale_colour_gradientn(
              colours = rev(brewer.pal(n = 11, name = "Spectral"))) +
              ggtitle(title)
  return(plot)
}

features <- c(
            "Cd3d", "Cd3e", "Cd4", # T cell 
            "Lef1",                # Memory T cell
            "Gzmb", "Cd8a",        # Cytotoxic / CD8 T cells 
            "Lag3",                # Exhausted CD4+ and CD8+ T cells
            "Cd79a", "Cd79b",      # B cells
            "Cd14",                # Monocytes
            "Ncr1")                # Natural Killer (NK cells)
        
titles <- c(
            "Cd3d: T cell", 
            "Cd3e: T cell", 
            "Cd4: T cell",
            "Lef1: T cell (Memory T cell)",
            "Gzmb: CD8+ T cells (Cytotoxic T cells)",
            "Cd8a: CD8+ T cells (Cytotoxic T cells)", 
            "Lag3: Exhausted CD4+ and CD8+ T cells",
            "Cd79a: B cell", 
            "Cd79b: B cell",
            "Cd14: Monocytes",
            "Ncr1: Natural Killer (NK) cells")
              

for (i in seq_along(features)) {
  plot_obj <- feFun(imm, feature = features[i], title = titles[i])
  print(plot_obj)
}
```

## 1.4 Preliminary NK and Myeloid Markers

```{r message=FALSE, warning=FALSE}
# Study Marker Genes
features <- c(
            "Gzmk", "Klrb1c", #"Nkg2d", # NK cell 
            "Csf1r", "Adgre1") #"Il6r") # Myeloid

        
titles <- c(
            "NK ", 
            "NK", 
            "Myeloid",
            "Myeloid")
              

for (i in seq_along(features)) {
  plot_obj <- feFun(imm, feature = features[i], title = titles[i])
  print(plot_obj)
}
```

```{r}
# No condition
FeaturePlot(imm, features = c(
            "Gzmk", "Klrb1c",  # NK cell 
            "Csf1r", "Adgre1"), # Myeloid
            max.cutoff = 3, 
            cols = c("grey", "red"))

# Split by sample condition
FeaturePlot(imm, features = c(
            "Gzmk", "Klrb1c",  # NK cell 
            "Csf1r", "Adgre1"), # Myeloid
            split.by = "sample", 
            max.cutoff = 3, 
            cols = c("grey", "red"))
```

```{r}
features <- c("Gzmk", "Klrb1c",  # NK cell 
            "Csf1r", "Adgre1")

RidgePlot(imm, features = features, ncol = 2)
```

# 2. Annotate Clusters

```{r Annotation, start_timer=TRUE, name = "Annotate_Cell_Type_chunk", message=FALSE, warning=FALSE}
# Define reference dataset
imm_ref <- SingleR::ImmGenData()

# Create object
sce <- as.SingleCellExperiment(DietSeurat(imm)) 

# Perform reference and extract label 
# Performs a fine and main annotation call. Main offers a strong annotation of cell types in clusters without over-annotating and breaking up clusters
label.main <- SingleR(test = sce, ref = imm_ref, labels = imm_ref$label.main)
label.fine <- SingleR(test = sce, ref = imm_ref, labels = imm_ref$label.fine)
```

```{r}
table(Label=label.main$labels,
      Lost=is.na(label.main$pruned.labels))
```

## 2.1 Visuals

### 2.1a Score Heatmap

```{r message=FALSE, warning=FALSE}
# Generate heatmap of cluster cell types
plotScoreHeatmap(label.main)
```

### 2.1b DeltaDistribution

```{r message=FALSE, warning=FALSE}
plotDeltaDistribution(label.main)
```

### 2.1c Heatmap

```{r}
tab <- table(cluster=imm$RNA_snn_res.0.3, label=label.main$labels) 
pheatmap::pheatmap(log10(tab+10))
```

## 2.2 Extract Meta Labels

```{r}
imm$cell.type.main <- label.main$pruned.labels
imm$cell.type.fine <- label.fine$pruned.labels
```

```{r}
# Set idents of original seurat object with annotated labels
imm <- SetIdent(imm, value = label.main$pruned.labels)
levels(imm)
```

```{r}
# Visualise
options(ggrepel.max.overlaps = Inf)
DimPlot(imm, label = T , repel = T, label.size = 3)
UMAPPlot(imm, split.by = 'sample', label = T, repel = T, label.size = 3) & NoAxes() & NoLegend()
```

## 2.3 Subset Cell Groups

```{r}
# Define lymphoid cells
Lymphoid_list <- c("T cells", "B cells", "NK cells", "NKT", "Tgd", "ILC")

Lymphoid <- subset(x = imm, 
                       idents = c(Lymphoid_list))

Lymphoid <- RunUMAP(Lymphoid, 
                       dims = 1:10, 
                       n.neighbors=5)

DimPlot(Lymphoid, label = T, repel = T)
## Subset Lymphocytes
# T cells ( T lymphocytes)
# NK cells (Natural Killer)
# ILC (Innate Lymphoid Cells)
# NKT (Natural Killer T cells)
# B cells (B lymphocytes)
# Tgd (Gamma-Delta T cells)

# Define myeloid cells
Myeloid_list <- c("Monocytes", "Neutrophils", "Eosinophils", "DC", "Mast cells", 
                  "Macrophages") 

Myeloid <- subset(x = imm, 
                       idents = c(Myeloid_list))

Myeloid <- RunUMAP(Myeloid, 
                       dims = 1:10, 
                       n.neighbors=5)

DimPlot(Myeloid, label = T, repel = T)
## Subset Myeloid
# Macrophages
# Neutrophils
# Monocytes
# Mast cells
# Dendritic cells (DC)
# Eosinophils

# Define non_immune cells
non_immune_list <- c("Fibroblasts", "Stem cells", "Endothelial cells", 
                     "Stromal cells", "Microglia")

non_immune <- subset(x = imm, 
                       idents = c(non_immune_list))

DimPlot(non_immune, label = T, repel = T)
## Non-Immune Cells
# Fibroblasts
# Endothelial cells
# Stromal cells
# Microglia (a type of macrophage found in the brain)
# Stem cells
```

# DEG Analysis of NK and Myeloid cells

## NK and Myeloid Markers

This section contains the workflow process for finding markers and then applying them to a subset to visualise co-expression.

```{r}
# Find co-expressed cells for NK and Myeloid markers
genes <- c("Gzmk", "Klrb1c", "Nkg2d", # NK cell
          "Csf1r", "Adgre1")          # Myeloid

Gene_Markers <- FindAllMarkers(object = imm, assay = "RNA") %>%
                filter(p_val_adj < 0.05) %>%
                filter(gene %in% genes)

head(Gene_Markers)
```

```{r}
# Subset data to specified columns
subset_results <- Gene_Markers[, c("avg_log2FC", "p_val_adj", "cluster")]
head(subset_results)
```

## NK and Myeloid Subset

```{r}
NK_Myeloid <- subset(x = imm, 
                       idents = c(Myeloid_list, "NK cells"))

NK_Myeloid <- RunPCA(NK_Myeloid,
                     dims = 1:10, 
                     n.neighbors=5)

DimPlot(NK_Myeloid, label = T, repel = T)

nk_markers <- c("Gzmk", "Klrb1c", "Nkg2d")
myeloid_markers <- c("Csf1r", "Adgre1", "Il6r")

VlnPlot(NK_Myeloid, features = nk_markers, slot = "counts", log = TRUE)
VlnPlot(NK_Myeloid, features = myeloid_markers, slot = "counts", log = TRUE)
```

## Visuals

```{r}
# Define gene markers 
features <- c("Gzmk", "Klrb1c", "Nkg2d",  # NK cell 
              "Csf1r", "Adgre1", "Il6r")  # Myeloid

# Visualise
RidgePlot(NK_Myeloid, features = features, ncol = 2)
VlnPlot(NK_Myeloid, features = features)
DotPlot(NK_Myeloid, features = features)
DotPlot(NK_Myeloid, features = features, split.by = "sample")
FeaturePlot(imm, features = c("Gzmk", "Adgre1"), blend = TRUE)
```

## Conditional Based on Sample

```{r}
# Markers between conditions (sample)
imm$celltype.cnd <- paste0(imm$cell.type.main,'_', imm$sample)
Idents(imm) <- imm$celltype.cnd

# Find markers betwen samples
response <- FindMarkers(imm, 
                        ident.1 = 'NK cells_PI1', 
                        ident.2 = 'NK cells_AP1', 
                        assay = "RNA")
head(response)

FeaturePlot(imm, features = nk_markers, 
            split.by = "sample", 
            max.cutoff = 3, 
            cols = c("grey", "red"))

FeaturePlot(imm, features = myeloid_markers, 
            split.by = "sample", 
            max.cutoff = 3, 
            cols = c("grey", "red"))
```

# 3. Save Object

```{r}
# Define directory
# Rename output 
setwd("E:")
SaveH5Seurat(imm, filename = "03_k23_object.h5Seurat", overwrite = TRUE)
```

```{r}
# Nk and Myeloid subset
SaveH5Seurat(NK_Myeloid, 
             filename = "NK_Myeloid_sub_object.h5Seurat", 
             overwrite = TRUE)

# Myeloid cells subset
SaveH5Seurat(Myeloid, 
             filename = "Myeloid_sub_object.h5Seurat", 
             overwrite = TRUE)

# Lymphoid cells subset
SaveH5Seurat(Lymphoid, 
             filename = "Lymphoid_sub_object.h5Seurat",
             overwrite = TRUE)
```

# Additional Content (Optional)

## Example: Pseudo / Trajectory Analysis

```{r}
cds <- as.cell_data_set(NK_Myeloid)
fData(cds)$gene_short_name <- rownames(fData(cds))
```

```{r}
# Cluster cells
cds <- cluster_cells(cds)

# Learn graph
cds <- learn_graph(cds)
```

```{r}
# cds@colData
# Plot trajectories
p1 <- plot_cells(cds, 
                 color_cells_by = "cell.type.main", 
                 show_trajectory_graph = FALSE) +
                 theme(plot.title = element_text(hjust = 0.5)) 
                 ggtitle("Cell Type")
                 
p2 <- plot_cells(cds, 
                 color_cells_by = "cell.type.main") +
                 theme(plot.title = element_text(hjust = 0.5)) +
                 ggtitle("Cell Type With Trajectory")

p3 <- plot_cells(cds, 
                 color_cells_by = "partition", 
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE) +
                 theme(plot.title = element_text(hjust = 0.5)) +
                 ggtitle("Partition With Trajectory")
```

```{r}
# Pseudotime Trajectory Plots
cds <- order_cells(cds, root_pr_nodes='Y_4')
  
ordered_pseudotime_values <- cds$pseudotime_order
pseudotime_values <- pseudotime(cds)

p4 <- plot_cells(cds,
                 color_cells_by = "pseudotime",
                 label_cell_groups=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 graph_label_size=1.5) +
                 theme(plot.title = element_text(hjust = 0.5)) +
                 ggtitle("Pseudotime Trajectories")

p4
```

## 3D Trajectory

```{r}
# 3D Trajectory Plot (cluster)
# 3D models
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
```

```{r}
cds_3d <- order_cells(cds_3d, root_pr_nodes='Y_4')

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="cluster")

cds_3d_plot_obj 
```

```{r}
imm$pseudotime <- pseudotime(cds)
FeaturePlot(imm, features = "pseudotime")
```

```{r}
# Create plot(s) of gene expression for marker genes, plotting expression along pseudotime with a line going through average expression

# Specify LOESS regression (non-parametric regression method) to estimate the relationship between pseudotime and gene expression by fitting mutiple regression lines locally. 

qplot(NK_Myeloid$pseudotime,
      as.numeric(NK_Myeloid@assays$RNA@data["Klrb1c",]), 
      xlab="pseudotime", 
      ylab="Expression", 
      main="Klrb1c") +
      geom_smooth(se = FALSE, method = "loess") + # Specify LOESS regression 
      theme_bw()

qplot(NK_Myeloid$pseudotime,
      as.numeric(NK_Myeloid@assays$RNA@data["Adgre1",]), 
      xlab="pseudotime", 
      ylab="Expression", 
      main="Adgre1") +
      geom_smooth(se = FALSE, method = "loess") + # Specify LOESS regression 
      theme_bw()


```

```{r}
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Klrb1c", "Adgre1")))
cds_subset <- cds[my_genes,]
alt_plot <- plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )
```

```{r}
my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("Csf1r", "Adgre1" )))
cds_subset <- cds[my_genes,]
alt_plot <- plot_genes_in_pseudotime(cds_subset, color_cells_by = "monocle3_pseudotime" )
```

## Alternative Find Markers

```{r}
# Find markers between marker genes
nk_markers_list <- c("Gzmk", "Klrb1c", "Nkg2d")
myeloid_markers <- c("Csf1r", "Adgre1", "Il6r")

Myeloid_list <- c("Monocytes", "Neutrophils", "Eosinophils", "DC", "Mast cells", 
                  "Macrophages") 

Idents(imm) <- imm@meta.data$cell.type.main

fm_nk <- FindMarkers(imm, ident.1 = "NK cells", 
                     min.pct = 0.25, 
                     assay = "RNA") %>%
                     filter(p_val_adj < 0.05) %>%
                     tibble::rownames_to_column("Gene")

fm_my <- FindMarkers(imm, ident.1 = Myeloid_list, 
                     min.pct = 0.25, 
                     assay = "RNA") %>%
                     filter(p_val_adj < 0.05) %>%
                     tibble::rownames_to_column("Gene")

fm_nk$Cell.Type <- "NK cell"
fm_my$Cell.Type <- "Myeloid"

co_expressed <- inner_join(fm_nk, fm_my, by = "Gene")

# Note: This method simplies finds the gene markers for NK cells, or Myeloids seperately after annotation
# Used as a method to explore if there are any overlaps in genes that contain both gene markers in both NK and Myeloid
```

```{r}
# Alternatively, we can use find markers between the idents of NK and Myeloid cells together.
# Used as a method to explore markers in comparison
test_markers <- FindMarkers(imm, 
                            ident.1 = "NK cells", 
                            ident.2 = Myeloid_list, 
                            min.pct = 0.25,
                            assay = "RNA") %>%
                            filter(p_val_adj < 0.05) %>%
                            tibble::rownames_to_column("Gene")
```

# Timer

```{r}
log.txt <- tic.log(format = TRUE)
log.txt
```
