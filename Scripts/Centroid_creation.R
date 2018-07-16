library("pamr")

classifier = meta_data$Deco_type
count_data_centroid = expr_raw

my_data = list(
  x = as.matrix(expr_raw),
  y = classifier,
  genenames = rownames( expr_raw ),
  geneid = as.character(1:dim( count_data_centroid )[1])
)

fit = pamr.train(
  my_data,
  n.threshold = 30#,
  #ngroup.survival = 5
)
min_threshold = tail(which( fit$errors == min(fit$errors)),1)

centroids = as.data.frame( fit$centroids,ncol=5)
centroids$means = rowMeans( centroids )
centroids = centroids[order( centroids$means, decreasing = T),]
centroids = centroids[,-ncol(centroids)]

final_centroid_list[,1] = rownames()
avg_importance = rowSums( abs((matrix(as.double(final_centroid_list[,3:6]), ncol = 4)) ))
final_centroid_list = cbind( final_centroid_list, avg_importance)
final_centroid_list = final_centroid_list[order(final_centroid_list[,ncol(final_centroid_list)],decreasing = T),]
write.table(final_centroid_list,"~/MAPTor_Net_RNA_data/Results/8685/Centroids_IT9.tsv",sep="\t",quote =F,row.names=F)

