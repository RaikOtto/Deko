# TPM count transformation
#gene_length_t = read.table("~/Deko/Misc/gene_length.tsv",sep ="\t", header = T, stringsAsFactors = F)
#l_match = match(rownames(count_data), gene_length_t$hgnc_symbol, nomatch = 0)
#count_data = count_data[ l_match != 0,]
#l_match = match( gene_length_t$hgnc_symbol, rownames(count_data), nomatch = 0)
#gene_length = gene_length_t$transcript_length[l_match]
#x = count_data / gene_length
#count_data = t(x) * 1e6 / colSums(x)
#count_data = t(count_data)

#seg_meta = read.table("~/Deko/Misc/Segerstolpe_Meta_info.tsv", sep ="\t", header = T)

#count_data = t(count_data)