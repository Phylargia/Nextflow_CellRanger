before <- gc(reset = TRUE)
peakRAM::peakRAM(library(Seurat))

## Data Path
# Define pathway to input file
data_path <- "C:/Users/Kenod/Pictures/count/filtered_feature_bc_matrix"

## Read in Data
# Read in cellranger count aggregated matrix
count_matrix <- Read10X(data.dir = data_path)

# Remove excess variable
rm(data_path)
peakRAM(pbmc <- CreateSeuratObject(
    counts = count_matrix, 
    min.cells = 3, 
    min.genes = 200, 
    project = "pbmc_k23"))

peakRAM(pbmc <- NormalizeData(
    object = pbmc, 
    normalization.method = "LogNormalize",
    scale.factor = 10000))

peakRAM(pbmc <- FindVariableFeatures(
    object = pbmc, 
    mean.function = ExpMean,
    dispersion.function = LogVMR))

peakRAM(pbmc <- ScaleData(pbmc, features = all.genes))

peakRAM(pbmc <- RunPCA(
    object = pbmc,
    pc.genes = pbmc@var.genes, 
    do.print = TRUE, 
    pcs.print = 1:5, 
    genes.print = 5))

peakRAM(pbmc <- FindNeighbors(
    object = pbmc, 
    dims = 1:10))

peakRAM(pbmc <- FindClusters(
    object = pbmc, 
    resolution = c(0.3)))

peakRAM(pbmc <- RunTSNE(
    object = pbmc, 
    dims.use = 1:10, 
    do.fast = TRUE))

peakRAM(pbmc <- RunUMAP(
    object = pbmc, 
    dims = 1:10))

after <- gc(reset = TRUE)
