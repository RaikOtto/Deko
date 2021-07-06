meta_info_map = read.table("~/MAPTor_NET/Misc/Meta_information.tsv",sep ="\t",header = T)
rownames(meta_info_map) = meta_info_map$Sample

meta_info_deko = read.table("~/Deko_Projekt/Misc/Meta_information.tsv",sep ="\t",header = T)
rownames(meta_info_deko) = meta_info_deko$Sample

matcher = match(meta_info_deko$Sample,meta_info_map$Sample,nomatch = 0)
meta_info_deko[matcher != 0,"NEC_NET"] = meta_info_map[matcher,"NEC_NET_PCA"]

meta_info_deko[matcher != 0,"NEC_NET"]

#write.table(meta_info_deko,"~/Deko_Projekt/Misc/Meta_information.tsv",sep ="\t", quote=F,row.names = F)
