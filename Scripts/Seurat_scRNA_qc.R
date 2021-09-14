library("stringr")
library("Seurath")
library("dplyr")
library("Matrix")

meta_info = read.table("~/Dek")
scrna_raw = load

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "~/Downloads/filtered_gene_bc_matrices/hg19/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size 

sparse.size <- object.size(x = pbmc.data) 
sparse.size

dense.size/sparse.size

pbmc <- CreateSeuratObject(
    raw.data = pbmc.data,
    min.cells = 3,
    min.genes = 200,
    project = "Tosti"
)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)

percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")

pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
