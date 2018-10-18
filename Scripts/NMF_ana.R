library("NMF")
library("stringr")

#type_vec = as.factor( c(1,2,3,4,5) )

#rand_m = matrix(1:16,ncol = 5)
nmf_res = nmf(count_data, 5,nrun=8, .opt='vP8')

basis(nmf_res)
coef(nmf_res)

round(basis(nmf_res) %*% coef(nmf_res),1)

summary( nmf_res, target= type_vec)

#init = nmfModel(4,W = .5)

