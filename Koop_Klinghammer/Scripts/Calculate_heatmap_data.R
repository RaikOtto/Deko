library("stringr")
library("ggplot2")

s_match = match( as.character( colnames(pure_data)), as.character( meta_info$Name), nomatch = 0)
meta_data = meta_info[s_match,]
rownames(meta_data) = meta_data$Name

threshold = 0.05
meta_data$Sig[ which( meta_data$P_value < threshold ) ] = "Sig"
meta_data$Sig[ which( meta_data$P_value >= threshold ) ] = "Not_Sig"
meta_data$Subtype[ which( meta_data$P_value >= threshold ) ] = "Not_sig"

#meta_data$Sig[ (as.character(meta_data$Erhaltungstherapie) != "2" ) ] = "Not_sig"

meta_data = meta_data[meta_data$Sig == "Sig",]
pure_data = pure_data[,colnames(pure_data) %in% meta_data$Name]

meta_data$OS = str_replace_all(meta_data$OS, pattern = ",", ".")
meta_data$OS = as.double(meta_data$OS)
meta_data$OS = as.double(meta_data$OS)
meta_data$OS = log(meta_data$OS)

meta_data$PFS = str_replace_all(meta_data$PFS, pattern = ",", ".")
meta_data$PFS = as.double(meta_data$PFS)
meta_data$PFS = as.double(meta_data$PFS)
meta_data$PFS = log(meta_data$PFS)

#meta_data$Subtype_2[! meta_data$Sig_2] = "Not_sig"
#meta_data$P_value_2 = 1-meta_data$P_value_2

#hc.cut <- factoextra::hcut(
#  scale(t(pure_data)),
#  k = 3,
#  hc_method = "ward.D"
#)
#factoextra::fviz_dend(hc.cut, show_labels = FALSE, rect = TRUE)
#factoextra::fviz_cluster(hc.cut, ellipse.type = "convex")
#meta_data$Subtype = hc.cut$cluster

###

aka3 = list(
  Subtype   = c(
      BA = "red",
      CL = "darkgreen",
      MS = "darkblue",
      Not_sig = "black"
  ),
  OS = c(
      low = "white",
      middle = "gray",
      high = "black"
  ),
  PFS = c(
    low = "white",
    middle = "gray",
    high = "black"
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
  #groups = as.character( meta_data$Erhaltungstherapie ),
  groups = as.character( meta_data$Subtype ),
  ellipse = TRUE, 
  circle = TRUE,
  var.axes = F
)+ theme(legend.position = "top")

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Loc_Primaertumor","PFS","OS","Best_response","Erhaltungstherapie","Subtype")],
  show_rownames = F,
  show_colnames = F,
  clustering_method = "complete",
  annotation_colors = aka3
)

OS = as.double(as.character(meta_data$OS))
Subtype = meta_data$Subtype[!(is.na(OS))]
response_vec = meta_data$Best_response[!(is.na(OS))]
pfs_vec = meta_data$OS[!(is.na(OS))]

OS = OS[!(is.na(OS))]
OS = OS[ (response_vec != "NE") & (response_vec != "")]
pfs_vec = pfs_vec[ (response_vec != "NE") & (response_vec != "")]
pfs_vec = as.double(pfs_vec)
Subtype = Subtype[ (response_vec != "NE") & (response_vec != "")]

response_vec = response_vec[ (response_vec != "NE") & (response_vec != "")]
aggregate(pfs_vec, by = list(Subtype), FUN = mean)
library(survival)
fit <- survfit(
  Surv(OS) ~ Subtype
)
survminer::ggsurvplot(fit, data = meta_data, risk.table = FALSE)

tab_mat = data.frame(cbind(os_vec, subtype_vec))
table(tab_mat)
fisher.test(table(tab_mat))

### GOI

GOI = c("SPRR3","KRT13","KRT1","CCL19","CXCL9","CCR7","VNN2")

dd = as.data.frame( t( pure_data[rownames(pure_data) %in% GOI,] ) )
dd$Subtype = meta_data$Subtype
dd2 = reshape2::melt( dd, value.name = "Subtype" )
colnames( dd2 )  = c( "Subtype", "Gene", "Exp" )

ggplot( dd2, aes( Gene, Exp, fill = Subtype ) ) + geom_boxplot( position = "dodge" ) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top")

###

data = table(meta_data[,c("Subtype","Best_response")])
fisher.test(data)

###

data = aggregate(meta_data$Chemozyklen, by = list(meta_data$Subtype), FUN = c)
t.test( 
  as.double( as.character(unlist(data[ data$Group.1 == "BA",2] ) ) ),
  as.double( as.character(unlist(data[ data$Group.1 == "MS",2] ) ) )
)

vis_mat = reshape2::melt(meta_data[,c("Subtype","Best_response")])
ggplot( vis_mat, aes( Subtype, Best_response ) ) + geom_bar(aes( x = Subtype, y = Best_response, fill = Best_response),stat="identity" ) + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position = "top")

sub_vec = meta_data$Subtype[s_match]
r_mat[1:5,1:5]

exp_mat = apply(
  r_mat[ , !( colnames(r_mat) %in% c("Sample", "OS") )], 
  FUN = function(vec){
     return(
        aggregate(vec, by = list(sub_vec) , c)
     )
  },
  MARGIN = 1
)
