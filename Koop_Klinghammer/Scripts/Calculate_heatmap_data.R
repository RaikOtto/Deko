library("stringr")
library("ggplot2")

# pure data pre-processing

#pure_data = read.table("~/Koop_Klinghammer/Data/35S.14_03_2018.normalized.tsv", sep ="\t", header = T, row.names = 1, stringsAsFactors = F)
colnames(pure_data) = str_replace_all(colnames(pure_data), pattern = "^X", "")
pure_data = pure_data[, !( colnames(pure_data) %in% c("25") ) ]
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t",header =T, stringsAsFactors = F)
meta_data = meta_info[match(colnames(pure_data), meta_info$Name),]

meta_data$Subtype[! meta_data$Sig] = "Not_sig"
meta_data$P_value = 1-meta_data$P_value
meta_data$Subtype_2[! meta_data$Sig_2] = "Not_sig"
meta_data$P_value_2 = 1-meta_data$P_value_2

###

aka3 = list(
  Subtype   = c(
      BA = "red",
      CL = "darkgreen",
      MS = "darkblue",
      Not_sig = "black"
  ),
  Subtype_2   = c(
    BA = "red",
    CL = "darkgreen",
    MS = "darkblue",
    Not_sig = "black"
  ),
  AREG = c(low = "white", middle= "gray", high = "black"),
  IDO1 = c(low = "white", middle= "gray", high = "black"),
  VEGF = c(low = "white", middle= "gray", high = "black")
)

### vis

cor_mat = cor(pure_data)
pca = prcomp(t(cor_mat));

ggbiplot::ggbiplot(
  pca,
  obs.scale = 1,
  var.scale = 1, 
  labels.size = 4,
  labels = colnames(pure_data),
  alpha = 1,
  groups = meta_data$Subtype,
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)+ theme(legend.position = "top")

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data["Subtype"],
  show_rownames = F,
  show_colnames = F,
  clustering_method = "average",
  annotation_colors = aka3
)

# tSNE

colors =as.character( meta_data$Subtype)
colors[ colors == "MS"] = "blue"
colors[ colors == "BA"] = "red"
colors[ colors == "CL"] = "green"
colors[ colors == "Not_sig"] = "black"
tsne_out <- Rtsne::Rtsne( t(pure_data), perplexity = 3) # Run TSNE
plot( tsne_out$Y, col = colors)

### GOI

GOI = c("AKR1C1","AKR1C3","AREG","CDKN2A","CXCL10","E2F2","EGFR","HIF1A","IDO1","IFNG","IL17A","KRT9","LAG3","STAT1","VEGF")

dd = as.data.frame( t( pure_data[rownames(pure_data) %in% GOI,] ) )
dd$Subtype = meta_data$Subtype
dd2 = reshape2::melt( dd, value.name = "Subtype" )
colnames( dd2 )  = c( "Subtype", "Gene", "Exp" )

ggplot( dd2, aes( Gene, Exp, fill = Subtype ) ) + geom_boxplot( position = "dodge" ) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top")

