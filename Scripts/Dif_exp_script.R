library("stringr")
library("limma")

meta_data = meta_info[colnames(count_data),]
meta_data[1:5,1:5]

groups = meta_data$Subtype
groups[groups == "HSC"] = "CASE"
groups[groups != "CASE"] = "CTRL"
design <- model.matrix(~0 + groups)
colnames(design) = c("HSC", "Not_HSC")

vfit <- lmFit(count_data,design)
contr.matrix = makeContrasts( contrast = HSC - Not_HSC,  levels = design )
vfit <- contrasts.fit( vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))

result_t = topTable( efit, coef = "contrast", number  = nrow(count_data), adjust  ="none", p.value = 1, lfc = 0)
result_t$hgnc_symbol = rownames(result_t)
colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")

result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
result_t = result_t[order(result_t$P_value, decreasing = F),]
result_t$Log_FC = round(result_t$Log_FC, 1)
result_t$Average_Expr = round(result_t$Average_Expr, 1)
result_t = result_t[order(result_t$Log_FC,decreasing = T),]

#result_t$ABS_FC = abs(result_t$Log_FC)
result_t = result_t[order(result_t$Log_FC, decreasing = T),]

write.table("~/Deko//Results/Dif_Exp/Dif_exp_HSC_vs_Not_HSC_Merge_Baron_Stanescu_Yan_scRNA.tsv", x = result_t, sep = "\t", quote = F, row.names = F)