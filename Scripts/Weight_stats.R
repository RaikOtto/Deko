library("ggplot2")
library("stringr")

stats_t = read.table("~/Deko/Misc/Weight_stats.tsv",sep ="\t", header = T, stringsAsFactors = F)
stats_t = stats_t[!is.na(stats_t$Krafthase) | !is.na(stats_t$Beaver),]

vis_mat = reshape::melt(stats_t[,c("Nr","Krafthase","Beaver")], id.vars = "Nr")
colnames(vis_mat) = c("Nr","Person","Weight")

# Warnings of mismatched lengths are suppressed
vis_mat_krafthase = subset(vis_mat, Person == "Krafthase")
p <- ggplot(vis_mat_krafthase, aes(x = Nr, y = Weight))
p <- p + geom_line(aes(y = Weight, colour = "Krafthase"))
p = p + geom_smooth(method = "lm", colour = "Red") + ylab("Krafthase") + xlab("Day")

matcher = 1.32
anti_matcher = 1 / matcher
p = p + scale_y_continuous(sec.axis = sec_axis(~.*matcher, name = "Beaver"))

# adding the relative humidity data, transformed to match roughly the range of the temperature
vis_mat_beaver = subset(vis_mat, Person == "Beaver")
vis_mat_beaver$Weight = vis_mat_beaver$Weight * anti_matcher
p = p + geom_line(data = vis_mat_beaver,aes(y = Weight, colour = "Beaver"))
p <- p + scale_color_manual( values = c("darkblue","Red") )
p <- p + theme(legend.position = "top")
p = p + geom_smooth(data= vis_mat_beaver,method = "lm", colour = "darkblue")

lm_beaver = lm( data = stats_t, Nr ~ Beaver)
summary(lm_beaver)
daily_change_beaver = round((stats_t$Beaver[1] - stats_t$Beaver[length(stats_t$Beaver)]) / length(stats_t$Beaver),2)
days_remaining_beaver = round((stats_t$Beaver[length(stats_t$Beaver)] - 78) / daily_change_beaver,1)
p = p + annotate("text", x = 12, y = 65.5, label = 
    paste(collapse = "", c("Daily change Beaver: ",daily_change_beaver, " Kg"))
)
p = p + annotate("text", x = 12, y = 65.0, label = 
    paste(collapse = "", c("Days until aspired Beaver weight of 78Kg: ",days_remaining_beaver))
)

lm_krafthase = lm( data = stats_t, Nr ~ Krafthase)
summary(lm_krafthase)
daily_change_krafthase = round((stats_t$Krafthase[1] - stats_t$Krafthase[length(stats_t$Krafthase)]) / length(stats_t$Krafthase),2)
days_remaining_krafthase = round((stats_t$Krafthase[length(stats_t$Krafthase)] - 60) / daily_change_krafthase,1)
p = p + annotate("text", x = 7, y = 62.0, label = 
    paste(collapse = "", c("Daily change Krafthase: ",daily_change_krafthase, " Kg"))
)
p = p + annotate("text", x = 8, y = 61.5, label = 
    paste(collapse = "", c("Days until aspired Krafthase weight of 60Kg: ",days_remaining_krafthase))
)
p