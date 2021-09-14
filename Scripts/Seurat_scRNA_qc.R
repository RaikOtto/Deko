library("stringr")
library("Seurat")
library("dplyr")
library("Matrix")

meta_info = read.table("~/Deko_Projekt/Misc/Tosti_Metadaten.tsv",sep ="\t", header = T)
rownames(meta_info) = meta_info$Cell
scrna_raw = readRDS("~/Dropbox/Tosti.scRNA.S112563.RDS")

# Examine the memory savings between regular and sparse matrices
dense.size = object.size(x = as.matrix(x = scrna_raw))
dense.size 

sparse.size <- object.size(x = scrna_raw) 
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

pbmc = FilterCells(
    object = pbmc,
    subset.names = c("nGene", "percent.mito"),
    low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05)
)

pbmc = NormalizeData(
    object = pbmc,
    normalization.method = "LogNormalize",
    scale.factor = 10000
)

pbmc = FindVariableGenes(
    object = pbmc,
    mean.function = ExpMean,
    dispersion.function = LogVMR,
    x.low.cutoff = 0.0125,
    x.high.cutoff = 3,
    y.cutoff = 0.5
)

length(x = pbmc@var.genes)

pbmc = ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))