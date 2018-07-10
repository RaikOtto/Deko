library("NMF")
library("stringr")

#nmf_res = nmf(t(count_data), 4,nrun=5, .opt='vP8')
nmf_res = readRDS("~/Deko/Data/NMF.06_07_2018.RDS")

W = as.data.frame(basis(nmf_res))
H = as.data.frame(coef(nmf_res))

colnames(W)[which.max(W[rownames(W) == "INS",])] = "Beta"
colnames(W)[which.max(W[rownames(W) == "GCG",])] = "Alpha"
colnames(W)[which.max(W[rownames(W) == "PPY",])] = "Gamma"
colnames(W)[which.max(W[rownames(W) == "SST",])] = "Delta"
colnames(W)[which.max(W[rownames(W) == "TTR",])] = "Alpha_five"

W[1:5,1:5]
H[1:5,1:5]

pred_labels = apply( H, MARGIN = 2, FUN = function(vec){
    return( colnames(W)[which.max(vec)] )})
pred_labels[1:5]

comparison_t = rbind( subtypes, pred_labels )
table(subtypes,pred_labels)

# prediction 
# W %*% H = V
# H = V %/% W

expr_raw[1:5,1:5]

W_match = match( rownames(W), rownames(count_data), nomatch = 0)
V_new = count_data[ W_match,]

V_match = match( rownames(W), rownames(V_new), nomatch = 0)
W_new   =  W[W_match != 0 ,]

W_new[1:5,1:5]
V_new[1:5,1:5]

dim(W_new)
dim(V_new)

H_new = MASS::ginv(as.matrix(W_new)) %*% as.matrix(V_new)
H_new[1:5,1:5]

pred_labels = apply( H_new, MARGIN = 2, FUN = function(vec){
  return( colnames(W_new)[which.max(vec)] )})

table(pred_labels)
