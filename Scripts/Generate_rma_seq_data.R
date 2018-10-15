library(powsimR)

expr_raw = read.table("~/MAPTor_NET/BAMs/Kallisto_three_groups/Groetzinger_Scarpa.TPM.filtered.ENSEMBL.Voom.TMM.normalized.tsv",sep="\t", stringsAsFactors =  F, header = T)
colnames(expr_raw) = stringr::str_replace(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
test_dat = expr_raw[1:100,1:10]

## simulating bulk RNA-seq experiment
ngenes = nrow(test_dat)
nsamples = ncol(test_dat)
true.means <- 2^rnorm(ngenes, mean=8, sd=2)
true.dispersions <- 3/true.means + 0.1
sf.values <- rnorm(nsamples, mean=1, sd=0.1)
sf.means <- outer(true.means, sf.values, '*')
cnts <- matrix(
  rnbinom(
    ngenes*nsamples,
    mu=sf.means,
    size=1/true.dispersions),
    ncol=nsamples
)
rownames(cnts) = rownames(test_dat)
colnames(cnts) = colnames(test_dat)
## estimating negative binomial parameters
estparam <- estimateParam(
  countData=cnts,
  Distribution='NB',
  RNAseq="bulk",
  normalisation='MR',
  sigma=1.96
  
)
plotParam(estparam, annot=F)


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
