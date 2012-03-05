
############ New summaries and diagnostics for SpatioTemporal models
############ 	Started by Assaf, Winter 2012

#### quick t-statistic calcluation for parameter estimates

tstat=function(pardat,alphas=TRUE) {
outv=rep(NA,length(pardat$par.all))
tmp=solve(-pardat$hessian.all)
if(!is.na(tmp[1])) outv=pardat$par.all/sqrt(diag(tmp))
if (!alphas) return(outv)

return(data.frame(Estimate=pardat$par.all,Std.Err=pardat$par.all/outv,t.stat=outv))
}

#### For an 'estimateCV' list, pool together all point estimates, t-statistics and Hessian eigenvalues 
#### to a single structure.

CVbasics=function(cvout) {

parpar=cvout$res.all[[1]]$par.all

outmt=matrix(0,nrow=length(parpar),ncol=length(cvout$res.all))
outmei=outmt
rownames(outmt)=names(parpar)

outma=outmt


for (a in 1:length(cvout$res.all)) {

    outmt[,a]=tstat(cvout$res.all[[a]],alphas=FALSE)

    outma[,a]=cvout$res.all[[a]]$par.all

    outmei[,a]=eigen(-cvout$res.all[[a]]$hessian.all)$value
}
list(alphas=outma,tees=outmt,eigens=outmei)
}
