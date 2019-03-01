library("powsimR")

expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.ENSEMBL.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = stringr::str_replace(colnames(expr_raw), pattern = "^X", "")

test_dat = expr_raw[ rownames(expr_raw) %in%as.character(unlist(pancreasMarkers))  ,]

## simulating bulk RNA-seq experiment

new_mat <<- matrix(as.double(),nrow = nrow(test_dat), ncol = ncol(test_dat))

for (subtype in unique(subtypes) ){
    expr_sub_mat = count_data[ rownames(count_data) %in% as.character(unlist(pancreasMarkers)) , subtypes ==  subtype]
    expr_sub_mat[1:5,1:5]
    
    ngenes = nrow(expr_sub_mat)
    nsamples = ncol(expr_sub_mat)    
    
    row_mean = rowMeans(log(expr_sub_mat + 1))
    row_sd = apply( log(expr_sub_mat + 1), MARGIN = 1, FUN = sd)
    
    true.means <- 2^rnorm(ngenes, mean = row_mean, sd = row_sd)
    true.dispersions = 3/true.means + 0.1
    
    sf.values = rnorm(nsamples, mean = row_mean, sd= row_sd)
    sf.means  = outer(true.means, sf.values, '*')
    
    cnts <- matrix(
        rnbinom(
            ngenes*nsamples,
            mu   = sf.means,
            size = 1/true.dispersions),
        ncol = nsamples
    )
    rownames(cnts) = rownames(test_dat)
    colnames(cnts) = colnames(test_dat)
    ## estimating negative binomial parameters
    estparam <- estimateParam(
        countData = cnts,
        Distribution = 'NB',
        RNAseq = "bulk",
        normalisation = 'MR',
        sigma = 1.96
        
    )
    plotParam(estparam, annot=F)
    
}



## simulate 2 groups of samples
p.med <- function(x) rnorm(x, mean=0, sd=0.5) # medium DE differences
simcounts.med <- simulateCounts(
    n=c(6,6),
    ngenes=10000,
    p.DE=0.05,
    pLFC = p.med,
    p.B=0,
    bLFC=NULL,
    bPattern="uncorrelated",
    p.M=NULL, mLFC=NULL,
    params=estparam,
    size.factors="equal",
    spike=NULL,
    spikeIns=FALSE,
    downsample=FALSE,
    geneset=FALSE,
    sim.seed=34628,
    verbose=TRUE
)

p.minor <- function(x) rnorm(x, mean=0, sd=0.1) # minor differences
simcounts.minor <- simulateCounts(
    n=c(6,10), ngenes=10000,
    p.DE=0.01, pLFC = p.minor,
    p.B=0, bLFC=NULL, bPattern="uncorrelated",
    p.M=NULL, mLFC=NULL,
    params=estparam,
    size.factors="equal",
    spike=NULL, spikeIns=FALSE,
    downsample=FALSE, geneset=FALSE,
    sim.seed=34628, verbose=TRUE
)

plotCounts(simCounts = simcounts.med, Distance = 'euclidean', Scale = T, DimReduce = "PCA", verbose = T)

plotCounts(simCounts = simcounts.minor, Distance = 'euclidean', Scale = T, DimReduce = "PCA", verbose = T)
