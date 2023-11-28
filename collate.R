## remotes::install_github("ichcha-m/cophescan")
library(cophescan)
library(magrittr)
library(randomFunctions)
args=getArgs(default=list(d="data/simdata_v2",outfile="data/simdata_v2/collated.RData"))

files=list.files(args$d,full=TRUE,pattern="block")
message(length(files))

if(interactive())
  files=sample(files,100)

DATA=vector("list",length(files))
for(i in seq_along(files)) {
  (load(files[i]))
  ## str(datasets)
  allsnps=colnames(datasets[[1]]$beta)
  for(j in seq_along(datasets)) { # reorder so first snp (column) is causal for Hc in each set
    snps=c(datasets[[j]]$cv$Hc[1], setdiff(allsnps,datasets[[j]]$cv$Hc[1]))
    datasets[[j]]$beta=datasets[[j]]$beta[,snps]
    datasets[[j]]$vbeta=datasets[[j]]$vbeta[,snps]
  }
  beta=lapply(datasets, "[[", "beta") %>% do.call("rbind", .)
  vbeta=lapply(datasets, "[[", "vbeta") %>% do.call("rbind", .)
  truth=sub("_.*","",rownames(beta))
  z=beta/sqrt(vbeta)
  sd.prior=0.2
  r <- sd.prior^2/(sd.prior^2 + vbeta)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  lABF.Hc=lABF[,1]
  lABF.Ha=apply(lABF[,-1],1,logsumexp)
  DATA[[i]]=list(lABF.Hc=lABF.Hc, lABF.Ha=lABF.Ha, truth=sub("_.*","",rownames(beta)),
                 nsnps=rep(ncol(lABF),nrow(beta)))
}

################################################################################
save(DATA,file=args$outfile)
