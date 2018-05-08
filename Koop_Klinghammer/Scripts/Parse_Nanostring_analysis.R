require('NanoStringNorm')
library("stringr")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t",header =T, stringsAsFactors = F)
raw_data  = read.markup.RCC( rcc.path = "~/Koop_Klinghammer//Data/Raw_data/", rcc.pattern = "*.RCC")

### tmp
colnames(raw_data$x)[-seq(3)] = 1:length(colnames(raw_data$x)[-seq(3)])

sample_names = names(raw_data$header)
sample_names = str_replace(sample_names, pattern = "^X","")
sample_names = str_replace_all(sample_names, pattern = "\\.","_")
mis_samp =  sum( !( sample_names %in% meta_info$Raw_Name) )

if (  mis_samp != 0 ){
  new_mat = matrix( rep("", mis_samp * ncol(meta_info) ), ncol = ncol(meta_info), nrow = mis_samp)
  colnames(new_mat) = colnames(meta_info)
  meta_info = rbind( meta_info, new_mat )
  meta_info$Name = as.character(1:nrow(meta_info))
  meta_info$Raw_Name = sample_names
  write.table( meta_info, "~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t", quote = F, row.names = F)
}
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t",header =T, stringsAsFactors = F)
meta_info = subset( meta_info, Included == TRUE)

### normalization

#norm.comp.results.test = norm.comp(raw_data, verbose = T)
eset = NanoStringNorm::NanoStringNorm( 
  raw_data,
  CodeCount.methods = "sum",
  Background.methods = "mean.2sd",
  SampleContent.methods = "none",
  OtherNorm.methods = "vsn",
  take.log = T,
  round.values = T,
  return.matrix.of.endogenous.probes = F,
  #traits = meta_data,
  verbose = T
)

m    = matrix( as.character(unlist( eset$normalized.data)), nrow=  dim(eset$normalized.data)[1], ncol = dim(eset$normalized.data)[2])
info = m[,seq(3)]
data = matrix( as.double(m[,-seq(3)]), nrow=  dim(eset$normalized.data)[1], ncol = dim(eset$normalized.data)[2]-3)
data = round(data,1)

rownames(data) = rownames(eset$normalized.data)
col_labels = str_replace( colnames(eset$normalized.data)[-seq(3)], pattern = "^X", "") 

colnames(data) = col_labels
res  = cbind( info,data )
res = cbind(rownames(eset$normalized.data), res)

pure_data = as.character(res)[-seq(dim(eset$normalized.data)[1]*4)]
pure_data = matrix( as.double( pure_data ), nrow = dim(eset$normalized.data)[1] )
rownames( pure_data ) = info[ ,2 ]
pure_data = pure_data[ info[,1] == "Endogenous"  ,]
colnames(pure_data) = 1:ncol(pure_data)

s_match = match( colnames(pure_data), meta_info$Name, nomatch = 0)
meta_data = meta_info[s_match,]

### optional normalization

design <- model.matrix(~0 + meta_data$Subtype)
colnames(design) = c("BA","CL", "MS")

DGE = edgeR::DGEList(pure_data)
DGE = edgeR::calcNormFactors(DGE,method =c("TMM"))
v = limma::voom(DGE,design)

pure_data = DGE$counts
#pure_data = v$E
boxplot(pure_data)



# export
#meta_info$Subtype[!meta_info$Sig] = "Not_sig"
#write.table(meta_info,"~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t",row.names =F, quote =F)

meta_data$AKR1C1 = as.character( pure_data[rownames(pure_data) == "AKR1C1",] )
meta_data$AKR1C3 = as.character( pure_data[rownames(pure_data) == "AKR1C3",] )
meta_data$AREG = as.character( pure_data[rownames(pure_data) == "AREG",] )
meta_data$CDKN2A = as.character( pure_data[rownames(pure_data) == "CDKN2A",] )
meta_data$CXCL10 = as.character( pure_data[rownames(pure_data) == "CXCL10",] )
meta_data$E2F2 = as.character( pure_data[rownames(pure_data) == "E2F2",] )
meta_data$EGFR = as.character( pure_data[rownames(pure_data) == "EGFR",] )
meta_data$HIF1A = as.character( pure_data[rownames(pure_data) == "HIF1A",] )
meta_data$IDO1 = as.character( pure_data[rownames(pure_data) == "IDO1",] )
meta_data$IFNG = as.character( pure_data[rownames(pure_data) == "IFNG",] )
meta_data$IL17A = as.character( pure_data[rownames(pure_data) == "IL17A",] )
meta_data$KRT9 = as.character( pure_data[rownames(pure_data) == "KRT9",] )
meta_data$LAG3 = as.character( pure_data[rownames(pure_data) == "LAG3",] )
meta_data$STAT1 = as.character( pure_data[rownames(pure_data) == "STAT1",] )
meta_data$VEGF = as.character( pure_data[rownames(pure_data) == "VEGF",] )

meta_match = match(meta_data$Name, meta_info$Name, nomatch = 0)

meta_info$AKR1C1[meta_match] = as.character( meta_data$AKR1C1 )
meta_info$AKR1C3[meta_match] = as.character( meta_data$AKR1C3 )
meta_info$AREG[meta_match] = as.character( meta_data$AREG )
meta_info$CDKN2A[meta_match] = as.character( meta_data$CDKN2A )
meta_info$CXCL10[meta_match] = as.character( meta_data$CXCL10 )
meta_info$E2F2[meta_match] = as.character( meta_data$E2F2 )
meta_info$EGFR[meta_match] = as.character( meta_data$EGFR )
meta_info$HIF1A[meta_match] = as.character( meta_data$HIF1A )
meta_info$IDO1[meta_match] = as.character( meta_data$IDO1 )
meta_info$IFNG[meta_match] = as.character( meta_data$IFNG )
meta_info$IL17A[meta_match] = as.character( meta_data$IL17A )
meta_info$KRT9[meta_match] = as.character( meta_data$KRT9 )
meta_info$LAG3[meta_match] = as.character( meta_data$LAG3 )
meta_info$STAT1[meta_match] = as.character( meta_data$STAT1 )
meta_info$VEGF[meta_match] = as.character( meta_data$VEGF )

#write.table(meta_info, "~/Koop_Klinghammer/Misc/Meta_information.tsv", sep ="\t", quote = F, col.names = T)

