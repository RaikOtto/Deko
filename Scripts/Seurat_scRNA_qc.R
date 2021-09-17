library("stringr")
library("Seurat")
library("dplyr")
library("Matrix")

meta_info = read.table("~/Deko_Projekt/Misc/Tosti_Metadaten.tsv",sep ="\t", header = T)
rownames(meta_info) = meta_info$Cell

scrna_raw = readRDS("~/Downloads/Tosti.scRNA.S112563.RDS")
meta_data = meta_info[colnames(scrna_raw),]
age_exclusion_vec = which(meta_data$age<5)
scrna_raw = scrna_raw[,-age_exclusion_vec]
dim(scrna_raw)
meta_data = meta_info[colnames(scrna_raw),]

pbmc = CreateSeuratObject(
    counts = scrna_raw,
    min.cells = 3,
    min.genes = 200,
    project = "Tosti"
)

mito.genes = grep(pattern = "^MT-", x = hgnc_list, value = FALSE)

nUMI_vec = meta_data$UMI_count
names(nUMI_vec) = rownames(meta_data)
nGene = meta_data$Gene_count
names(nGene) = rownames(meta_data)
donor_id = meta_data$patient_ID
names(donor_id) = rownames(meta_data)
cell_type = meta_data$Cluster
names(cell_type) = meta_data$Cell

percent.mito = Matrix::colSums(scrna_raw[mito.genes, ])/Matrix::colSums(scrna_raw)
pbmc = AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pbmc = AddMetaData(object = pbmc, metadata = nUMI_vec, col.name = "nUMI")
pbmc = AddMetaData(object = pbmc, metadata = nGene, col.name = "nGene")
pbmc = AddMetaData(object = pbmc, metadata = donor_id, col.name = "donor_id")
pbmc = AddMetaData(object = pbmc, metadata = cell_type, col.name = "cell_type")

#VlnPlot(object = pbmc, features = c("nGene", "nUMI", "percent.mito","donor_id"), log = TRUE, split.plot = TRUE, group.by = donor_id)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")

pbmc <- subset(
    pbmc,
    subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito < 0.05
)

pbmc = NormalizeData(
    pbmc,
    normalization.method = "LogNormalize",
    scale.factor = 1e4)

pbmc = FindVariableFeatures(
    pbmc,
    selection.method = 'vst',
    nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes = rownames(pbmc)
length(all.genes)
length(unique(all.genes))
pbmc = ScaleData(pbmc, features = all.genes)

pbmc = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc[['pca']], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = 'pca')
DimPlot(pbmc, reduction = 'pca')

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

pbmc = RunUMAP(pbmc, dims = 1:14)
table(meta_data$Cluster)
DimPlot(pbmc, reduction = "umap")

beta_markers = FindMarkers(pbmc, ident.1 = "Beta", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)