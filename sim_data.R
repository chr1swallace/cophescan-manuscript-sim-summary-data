#!/home/cew54/localc/bin/Rscript
library(randomFunctions)
library(magrittr)
library(abind)
## args=getArgs(defaults=list(chr=2, begin=43082452, end=44086664), numeric=c("chr","begin","end"))
source("dirs.rb")
args=getArgs(defaults=list(block="chr2_block0",d="data/simdata_v2"))

blockfile=list.files(paste0(DIR,"/reference/byblock"),pattern=".vcf.gz") %>% sample(., 1)
args$block=sub(".vcf.gz","",blockfile)

################################################################################
## get reference haplotype data to swap alleles - use 1000GP hg38
message("reading haplotypes for block ",args$block)
tmp.imp=tempfile(tmpdir="data/tmp")
paste("./make_ref_block.rb ",args$block,tmp.imp) %>% system()

library(data.table)
leg=fread(paste0(tmp.imp, ".impute.legend"))
message("snps found: ",nrow(leg))
haps=scan(paste0(tmp.imp, ".impute.hap"),what=0) %>% matrix(., ncol=nrow(leg))
colnames(haps)=leg$ID

unlink(paste0(tmp.imp,".*"))

if(ncol(haps)>1000)
  haps=haps[,1:1000]
if(interactive() & ncol(haps)>100)
  haps=haps[,1:100] # for testing

################################################################################

library(simGWAS)
freq=(haps + 1) %>% as.data.frame()
freq$Probability <- 1/nrow(freq)
sum(freq$Probability)

### store MAF, LD
snps <- colnames(haps)
LD=cor(haps)
MAF=colMeans(haps) %>% pmin(., 1-.)

## functions

## inner simulation function
simall=function(cv,null=FALSE) {
  if(null==TRUE) {
    g=0; nrep=18
  } else {
    g <- rnorm(length(cv), 0.04); nrep=1
  }
  z <- f.z(cv,g,nrep)
  ## abs(z) %>% apply(., 1, max) %>% summary()
  ## (apply(abs(z), 1, max) > 4.89) %>% mean
  vbeta <- f.vb(cv,g,nrep)
  beta <- z * sqrt(vbeta)
  colnames(beta)=colnames(vbeta)=colnames(haps)
  list(beta=beta,vbeta=vbeta,g=g,cv=cv)
}

## selects CVs, calls simall to do sim, and tidies results
repsim=function(cv) {
  CV=sample(snps,4)
  Hc=simall(CV[1])
  Ha=simall(CV[2])
  Hc2=simall(CV[c(1,3)])
  Ha2=simall(CV[c(2,4)])
  Hn=simall(CV[1], null=TRUE)
  H=list(Hc,Ha,Hc2,Ha2,Hn)
  result=list(beta=lapply(H, "[[", "beta") %>% do.call("rbind",.),
              vbeta=lapply(H, "[[", "vbeta") %>% do.call("rbind",.),
              beta=lapply(H, "[[", "beta") %>% do.call("rbind",.),
              g=list(Hc=Hc$g,Ha=Ha$g,Hc2=Hc2$g,Ha2=Ha2$g),
              cv=list(Hc=Hc$cv,Ha=Ha$cv,Hc2=Hc2$cv,Ha2=Ha2$cv))
  rownames(result$beta)=rownames(result$vbeta)=c("Hc","Ha","Hc2","Ha2",paste0("Hn_",1:18))
  result
}

f.z=function(cv,g,nrep=1) {
            simulated_z_score(N0=10000, # number of controls
                            N1=10000, # number of cases
                            snps=snps, # column names in freq of SNPs for which Z scores should be generated
                            W=cv, # causal variants, subset of snps
                            gamma.W=g, # log odds ratios
                            freq=freq, # reference haplotypes
                            nrep=nrep)
}
f.vb=function(cv,g,nrep=1) {
            simulated_vbeta(N0=10000, # number of controls
                            N1=10000, # number of cases
                            snps=snps, # column names in freq of SNPs for which Z scores should be generated
                            W=cv, # causal variants, subset of snps
                            gamma.W=g, # log odds ratios
                            freq=freq, # reference haplotypes
                            nrep=nrep)
}

datasets=replicate(20, repsim(), simplify=FALSE)

save(datasets,MAF,LD,
     file=tempfile(pattern=paste0(args$block,"_"),tmpdir=args$d,fileext=".RData"))
