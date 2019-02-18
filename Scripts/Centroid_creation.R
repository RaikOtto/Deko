library("org.Hs.eg.db")
library("pamr")
library("stringr")

meta_info = read.table("~/Deko/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
rownames(meta_info) = meta_info$Name
colnames(meta_info) = str_replace(colnames(meta_info),pattern = "\\.","_")

count_data_centroid = read.table("~/Deko/Data/Human_differentiated_pancreatic_islet_cells_scRNA/Segerstolpe.tsv",sep="\t", header = T)
colnames(count_data_centroid) = str_replace(colnames(count_data_centroid), pattern = "\\.", "_")
meta_data = meta_info[colnames(count_data_centroid),]
subtype_vec = as.character(meta_data$Subtype)
table(subtype_vec)

count_data_centroid = count_data_centroid[ ,subtype_vec %in% c("Alpha","Beta","Gamma","Delta")]
meta_data = meta_info[colnames(count_data_centroid),]
subtype_vec = as.character(meta_data$Subtype)
table(subtype_vec)

count_data_centroid = t(apply( count_data_centroid, FUN = function(vec){return((as.double(unlist(vec))))} , MARGIN = 1 ))
count_data_centroid[1:5,1:5]


my_data = list(
  x = count_data_centroid,
  y = subtype_vec,
  genenames = rownames( count_data_centroid ),
  geneid = as.character(1:dim( count_data_centroid )[1])
)

fit = pamr.train(
  my_data,
  n.threshold = 30
)
min_threshold = tail(which( fit$errors == min(fit$errors)),1)

centroids = as.data.frame( fit$centroids,ncol=5)
centroids$means = rowMeans( centroids )
centroids = centroids[order( centroids$means, decreasing = T),]
centroids = centroids[,-ncol(centroids)]

# optics

pamr::pamr.plotcen(fit = fit, data = my_data,threshold = min_threshold)

final_centroid_list = pamr::pamr.listgenes(
  fit = fit,
  data = my_data,
  threshold = min_threshold,
  genenames = T
)
colnames(final_centroid_list) = str_replace_all(colnames(final_centroid_list) , pattern = "-score", "")

avg_importance = as.double(rowSums(
  abs(
    matrix(
      as.double(
        final_centroid_list[
          ,colnames(final_centroid_list) %in% c("Alpha","Beta","Gamma","Delta")
        ]
      ),ncol = 4
    )
  )
))
final_centroid_list = cbind( final_centroid_list, as.double(avg_importance))
final_centroid_list = final_centroid_list[order(abs(as.double(final_centroid_list[,7])),decreasing = T),]
#write.table(final_centroid_list,"~/Deko/Models/Centroid_Alpha_Beta_Gamma_Delta_Segerstolpe.tsv",sep="\t",quote =F,row.names=F)

