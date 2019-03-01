library(ggplot2)

#m_data = read.table("~/Deko/Data/GSE87375.TPM.HGNC.tsv", sep ="\t",header = T, stringsAsFactors = F)
m_data = read.table("~/Deko/Data/ProcessedData.tsv", sep ="\t",header = T, stringsAsFactors = F, row.names = 1)
colnames(m_data) = stringr::str_replace_all(colnames(m_data), pattern ="^X","")
m_data[1:5,1:5]

row_var = apply( m_data, FUN = var, MARGIN = 1)

ggplot2::qplot(sort((row_var),decreasing = T))
#sample_index = sample(1:ncol(expr_raw), 500)
#gene_index = order(row_var,decreasing = T)[1:300]

#sample_data = expr_raw[ gene_index , sample_index]
#rownames(sample_data) = rownames(expr_raw)[gene_index]

umap_plot = umap::umap(t(m_data))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")

dist_mat = dist((vis_data))
k_clust = hclust(dist_mat)
clusterCut <- cutree(k_clust, 3)
clusterCut[clusterCut == 1] = "darkgreen"
clusterCut[clusterCut == 2] = "red"
clusterCut[clusterCut == 3] = "brown"

sort((row_var),decreasing = T)
hgnc_gene_name = "NEUROG3"
neurog3 = (as.double(m_data[rownames(m_data) == hgnc_gene_name,]))
size_indicator = neurog3
p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, size = size_indicator, color = clusterCut)
p = p + ggtitle(hgnc_gene_name) +  xlab("umap x") + ylab("umap y") + geom_text(aes(label=names(clusterCut)),hjust=0, vjust=0)
p

meta_data = data.frame(
    "Cluster" = as.character(clusterCut),
    "Name" = as.character(names(clusterCut)),
    stringsAsFactors = F
)
rownames(meta_data) = names(clusterCut)
pheatmap::pheatmap(
  cor(m_data),
  show_rownames = F,
  show_colnames = F,
  annotation_col = meta_data["Cluster"]
)

#ggplot2::ggplot( data = vis_data, aes( x, y ), size = neurog3[sample_index])

#count_data = read.table("~/Deko/Data/Count_data.Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)
#count_data[1:5,1:5]

s_match = match( colnames(count_data), seg_meta$Extract.Name, nomatch = 0)
subtypes = as.character(seg_meta$Characteristics.cell.type.)[s_match]
subtypes = str_replace_all(subtypes, pattern = " cell", "")

count_data = count_data[ ,  subtypes %in% c("Alpha","Beta", "Gamma", "Delta") ]
subtypes = subtypes[ subtypes %in% c("Alpha","Beta", "Gamma", "Delta") ]
dim(count_data)
length(subtypes)

intersecting_genes = intersect(rownames(count_data), rownames(m_data))
new_seger = count_data[rownames(count_data) %in% intersecting_genes,]
new_botton = m_data[rownames(m_data) %in% intersecting_genes,]
new_seger = new_seger[order(rownames(new_seger)),]
new_botton = new_botton[order(rownames(new_botton)),]
dim(new_seger)
dim(new_botton)
new_seger[1:5,1:5]
new_botton[1:5,1:5]

new_merge = cbind(new_seger,new_botton)
#write.table(new_merge,"~/Deko/Data/Merge_Seger_Botton.tsv",sep ="\t", col.names = T, row.names = T, quote = F)
#new_merge = read.table("~/Deko/Data/Count_data.Segerstolpe.tsv",sep ="\t", header = T, stringsAsFactors = F)
dim(new_merge)
length(cutree(k_clust, 3))
subtypes = c(subtypes, cutree(k_clust, 3) )
subtypes[subtypes == 1] = "Botton_1"
subtypes[subtypes == 2] = "Botton_2"
subtypes[subtypes == 3] = "Botton_3"

dim(new_merge)
length(subtypes)

#subtypes = read.table("~/Deko/Data/Merge_Subtypes.tsv", header = F, sep = "\t", stringsAsFactors = F)[,1]
#write.table(subtypes, "~/Deko/Data/Merge_Subtypes.tsv", row.names = F, col.names = F, quote = F, sep ="\t")

aggregate(as.double(m_data[rownames(m_data) == "NGN3_nas",]), by = list(as.character(clusterCut)), FUN = mean)
