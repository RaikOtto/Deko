library("stringr")
library("survminer")
library("survival")

data_t = read.csv("~/Deko/Misc/Auswahl_Samples.tsv",sep="\t", row.names = 1,stringsAsFactors = F)

survival = as.double(str_replace(data_t$OS_Tissue, pattern = ",", "."))
censored = data_t$Zensur

target_vector = data_t$AADAT_LKyn_sum / data_t$TPHsum
target_vector = data_t$M_LKyn_c / data_t$M_Serotonin
target_vector = as.double(target_vector)

target_vector[target_vector < mean(target_vector)] = "low"
target_vector[target_vector != "low"] = "high"

data_t$survival = survival
data_t$target_vector = target_vector

fit = survfit(
    Surv( survival, data_t$Zensur ) ~ target_vector
)

ggsurvplot(
    fit,
    data = data_t,
    risk.table = T,
    pval = T,
    censor.size = 10
)

### boxplot

data = data_t[,c("target_vector","NEC_NET")]
vis_t = reshape2::melt(data)
colnames(vis_t) = c("NEC_NET","X","AADAT_TPH")
vis_t$Kyn_Sero = log(vis_t$AADAT_TPH)

p = ggplot( data = vis_t,aes( x = NEC_NET, y = AADAT_TPH, fill =NEC_NET))#, min = ROC-SD, max = ROC+SD) )
p = p + geom_boxplot()
p


t.test(
    subset(vis_t$AADAT_TPH,vis_t$NEC_NET == "NET"),
    subset(vis_t$AADAT_TPH,vis_t$NEC_NET == "NEC")
)
#p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
#p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
#p = p + annotate("text", x=1:57,y = 5.5,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
#p = p + xlab("") + ylab("MEN1 expression in log TPM and MEN1 mutation allele frequency") + theme(legend.position = "top")
p + geom_errorbar(aes(),  position = "dodge")
