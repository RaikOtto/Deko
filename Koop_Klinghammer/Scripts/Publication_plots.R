library("stringr")
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

### Prep

library("grid")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep = "\t",header = T,stringsAsFactors = F)
meta_data = meta_info[ match( colnames(pure_data), meta_info$Name), ]
rownames(meta_data) = meta_data$Name

###

aka3 = list(
  Subtype = c(BA = "blue", CL ="darkgreen", MS = "red" ),
  Location = c(Primary = "white", Metastasis = "black"),
  OS = c(high = "White", medium = "gray", low = "black"),
  Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
  Included = c(Yes = "green", No = "red"))

###

cor_mat = cor(expr);pcr = prcomp(t(pure_data))

## Figure 1

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Subtype")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7#,
  #clustering_method = "single"
)

## Same patient only

find_vec = meta_data$pID
match_vec = match(find_vec, unique(find_vec),nomatch = 0)
multi_match = which( table(match_vec) > 1  )
multi_match = multi_match[multi_match != 12]
true_match = meta_data$pID[which( match_vec %in% multi_match)]

multi_data = pure_data[,which( meta_data$pID %in% true_match )]
multi_cor = cor(multi_data)
meta_data$pID = as.character(meta_data$pID)

pheatmap::pheatmap(
  multi_cor,
  annotation_col = meta_data[c("Subtype","OS","pID")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

### pca

ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data$Subtype),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = meta_data$Name
) 

## Selected patient only
# pID = c(2, 362, 365, 397)

multi_meta_data = as.data.frame(meta_data[which( meta_data$pID %in%  365),])
multi_data = pure_data[, multi_meta_data$Name]
multi_cor = cor(multi_data)
pcr = prcomp(t(multi_cor))

pheatmap::pheatmap(
  multi_cor,
  annotation_col = meta_data[c("Subtype","OS","pID")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

ggbiplot::ggbiplot(
  pcr,
  groups = as.character(multi_meta_data$Subtype),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = multi_meta_data$Name
) 
