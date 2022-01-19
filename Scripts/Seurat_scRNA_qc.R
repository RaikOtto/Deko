library("stringr")
library("Seurat")
library("dplyr")
library("Matrix")

meta_info = read.table("~/Downloads/drive-download-20220114T151242Z-001/barcodes.tsv",sep =",", header = T)
rownames(meta_info) = meta_info$cell.name

#scrna_raw = readRDS("~/Downloads/Tosti.Seurat.normalized.S78048.RDS")
scrna_raw = read.table("~/Downloads/count-matrix.txt", sep =" ", header = TRUE)
scrna_raw[1:5,1:5]
meta_data = meta_info[match(colnames(scrna_raw), meta_info$cell.name,nomatch = 0),]
table(meta_data$patient)
table(meta_data$is_edge)

n_finder = str_detect(colnames(scrna_raw), pattern = "N")
colnames(scrna_raw)[n_finder]
n_finder = str_detect(meta_info$cell.name, pattern = "N")
meta_info[n_finder,]
"N1_AAACCTGCAACCGCCA" %in% colnames(scrna_raw)
"N1_AAACCTGCAACCGCCA" %in% meta_info$cell.name
colnames(scrna_raw)[str_detect(colnames(scrna_raw), pattern ="N")]

candidate_samples = which(colnames(scrna_raw) %in% meta_info$cell.name)
scrna.S12608 = scrna_raw[,candidate_samples]
meta_data = meta_info[match(colnames(scrna.S12608), meta_info$cell.name,nomatch = 0),]
table(meta_data$patient)
table(meta_data$is_edge)
saveRDS(scrna.S12608, "~/Downloads/Peng.S12608.RDS")

meta_data = meta_info[match(colnames(scrna.S12608), meta_info$cell.name,nomatch = 0),]
candidate_not_neoplastic = str_detect(meta_data$patient, pattern = "N")
scrna = scrna.S12608[,candidate_not_neoplastic]
dim(scrna)
saveRDS(scrna.S12608, "~/Downloads/Peng.S12608.RDS")
dim(scrna)

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
DimPlot(pbmc, reduction = "umap")

beta_markers = FindMarkers(pbmc, ident.1 = "Beta", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

pbmc.markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)

normalized_data = pbmc
normalized_data[1:5,1:5]

write.table(
    as.matrix(GetAssayData(object = pbmc, slot = "data")), 
    '~/Downloads/Tosti.Seurat.normalized.S78048.tsv', 
    sep = ',', row.names = T, col.names = T, quote = F)

#write.table(meta_data,"~/Downloads/Meta_info.Edge.tsv",sep ="\t",quote =F)


####

#expr_raw = readRDS("~/Downloads/Tosti.Seurat.normalized.S78048.RDS")
#rownames = as.character(sapply(rownames(expr_raw), FUN = function(vec){return(as.character(unlist(str_split(vec,pattern = "\\|")))[1])}))
rownames = as.character(sapply(rownames(expr_raw), FUN = function(vec){return(as.character(unlist(str_split(vec,pattern = "\\.")))[1])}))
rownames = str_to_upper(rownames(expr_raw))
hgnc_list = rownames

#### downsampling

meta_info = read.table("~/Deko_Projekt/Misc/Tosti_Metadaten.tsv",sep ="\t", header = T)
rownames(meta_info) = meta_info$Cell

expr_raw = readRDS("~/Downloads/Tosti.Seurat.normalized.S78048.RDS")
meta_data = meta_info[colnames(expr_raw),]

candidates_cells = which((meta_data$cluster == "Acinar cell") & (meta_data$is_edge == "1"))
meta_data[candidates_cells,"cluster"] = "Acinar_edge_cell"

cell_type_vec_overall = str_to_lower(meta_data$cluster)
table(cell_type_vec_overall)

#

#cell_type_vec_overall[cell_type_vec_overall %in% c("acinar-i","acinar-s","acinar-reg+")] = "acinar"
#cell_type_vec_overall[cell_type_vec_overall %in% c("muc5b+ ductal","ductal")] = "ductal"

#

#cell_type_vec = cell_type_vec_overall[!(cell_type_vec_overall %in% c("schwann","endothelial","activated stellate","quiescent stellate","macrophage"))]
cell_type_vec = cell_type_vec_overall#[!(cell_type_vec_overall %in% c("alpha","beta","gamma","delta","schwann","endothelial","activated stellate","quiescent stellate","macrophage"))]
table(cell_type_vec)
selected_samples = c()

for ( cell_type in unique(cell_type_vec)){
    coords = which(cell_type_vec_overall == cell_type )
    
    if (length(coords) >= 200)
        coords = sample(coords, size = 200)
    
    selected_samples = c(selected_samples, coords)
}
length(selected_samples)

expr = expr_raw[,selected_samples]
model_name = "Gopalan_S784"

library("devtools")
load_all("~/artdeco")
source("~/Deko_Projekt/CIBERSORT_package/CIBERSORT.R")
library("bseqsc")

subtypes = cell_type_vec_overall[selected_samples]
table(subtypes)

dim(expr)
length(subtypes)

add_deconvolution_training_model_bseqsc(
    transcriptome_data = expr,
    model_name = model_name,
    subtype_vector =  str_to_lower(subtypes),
    training_p_value_threshold = 0.05,
    training_nr_permutations = 0,
    training_nr_marker_genes = 200
)

