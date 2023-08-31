source('basicfunctions.R')
sl=ls()
n0_xy=floor(prop_gwas_overlap_Xs_and_Y*min(c(sample_size_Y,sample_size_Xs)))
nall=sample_size_Y+sample_size_Xs-n0_xy
indY=1:sample_size_Y
indX=(nall-sample_size_Xs+1):nall
G=rbinom(number_of_causal_SNPs*nall,2,mafs_of_causal_SNPs) # assuming independence for now
G=matrix(G,nr=nall,nc=number_of_causal_SNPs)
G=apply(G,2,std)
# order will go UHP, CHP, valid, weak
### CHP IVs
gammaC=rep(0,number_of_causal_SNPs)
if(number_of_CHP_causal_SNPs>0) {
  chpix=1:number_of_CHP_causal_SNPs
  gammaC_=runif(number_of_CHP_causal_SNPs,-1/2,1/2);
  adj=U_variance_explained_by_CHP/sum(gammaC_^2)
  gammaC_=sqrt(adj)*gammaC_; gammaC[chpix]=-gammaC_
} else {
  chpix=c()
}
### setting U
eU=rnorm(nall,0,sqrt(1-U_variance_explained_by_CHP))
U=G%*%gammaC+eU
### setting X
CorrXX=parthcorr(phenotypic_correlation_Xs,n=number_of_exposures)
GenCorrXX=parthcorr(genetic_correlation_Xs,n=number_of_exposures)
LD=parthcorr(LD_causal_SNPs,number_of_causal_SNPs)
K=kronecker(GenCorrXX,LD)
Thsq=chol(solve(LD))
B=rmvn(1,rep(0,dim(K)[1]),K) # can effectively add LD the G by adding LD to B
B=matrix(B,nr=number_of_causal_SNPs,nc=number_of_exposures)
th=chol(solve(GenCorrXX))
cop=pnorm(B%*%th)
#cor(cop)
cop=cop%*%chol(GenCorrXX)
list(gencor=round(GenCorrXX,2),cop=round(cor(cop),2),corB=round(cor(B),2))
# rescale to match heritability
adj=Xs_variance_explained_by_g/colSums(B^2)
for(i in 1:ncol(B)) B[,i]=sqrt(adj[i])*B[,i]
# effect of confounder
pix=sqrt(Xs_variance_explained_by_U)
sdeX=diag(1-Xs_variance_explained_by_g-pix^2,number_of_exposures)
sdeX=sqrt(sdeX)
SigmaEX=sdeX%*%CorrXX%*%sdeX
eX=rmvn(nall,rep(0,number_of_exposures),SigmaEX)
X=G%*%B+matrix(pix*U,nr=nall,nc=number_of_exposures)+eX
### model for Y
vXY=Y_variance_explained_by_Xs
theta=vXY*rep(1,number_of_exposures)*(length(vXY)==1)+vXY*(length(vXY)>1)
theta=theta*signs_of_causal_effects
adj=vXY/sum(theta^2)
theta=sqrt(adj)*theta
piy=-1*sqrt(Y_variance_explained_by_U) # always negative
gammaU=rep(0,number_of_causal_SNPs)
if(number_of_UHP_causal_SNPs>0 & length(chpix)>0) {
  uhpix=(max(chpix)+c(1:number_of_UHP_causal_SNPs))
} else if(number_of_UHP_causal_SNPs>0 & length(chpix)==0){
  uhpix=1:number_of_UHP_causal_SNPs
} else {
  uhpix=c()
}
if(number_of_UHP_causal_SNPs>0) {
  gammaU_=runif(number_of_UHP_causal_SNPs,-1,1);
  adj=Y_variance_explained_by_UHP/sum(gammaU_^2)
  gammaU_=sqrt(adj)*gammaU_; gammaU[uhpix]=gammaU_
}
vUHPY=Y_variance_explained_by_UHP
eY=rnorm(nall,0,sqrt(1-vXY-vUHPY-piy^2))
Y=X%*%theta+piy*U+G%*%gammaU+eY
RhoXY=cor(cbind(Y,X))
Sinv=solve(RhoXY[2:ncol(RhoXY),2:ncol(RhoXY)])
p=ncol(B);nY=length(indY); nX=sample_size_Xs
nn=sqrt(c(nY,rep(nX,p)));nn=nn%*%t(nn) # for RhoME (correlations between measurement errors)
no=nn*0+n0_xy;no[1,1]=nY;no[2:(p+1),2:(p+1)]=nX
RhoME=no/nn*RhoXY
### GWAS
gwas_y=biggwas(Y[indY],G[indY,]) # indY created near the top
by=gwas_y$est; byse=gwas_y$std
bx=bxse=matrix(nr=number_of_causal_SNPs,nc=number_of_exposures)
for(j in 1:ncol(bx)) {
  fit=biggwas(X[indX,j],G[indX,])
  bx[,j]=fit$est
  bxse[,j]=fit$std
}
### uncomment if you want to see plots
# par(mfrow=c(1,2))
# pfun(B%*%theta,B%*%theta+piy*gammaC+gammaU,chpix,uhpix,xlab='bx*theta',ylab='by')
# pfun(bx%*%theta,by,chpix,uhpix,xlab='bxhat*theta',ylab='byhat')
### standardize?
bxunstd=bx; bxseunstd=bxse
data=parthstd(bx,by,bxse,byse,mafs_of_causal_SNPs,sample_size_Xs,sample_size_Y,MR_standardization_type)
bx=data$bx;bxse=data$bxse;by=data$by;byse=data$byse; m=nrow(bx)
### instrument selection based on P-value
######
### any P-value pruning?
st=tolower(substr(simtype,start=1,stop=3))
if(st=='win') {
  keep=pjs=c()
  Ruu=RhoME[-1,-1] # take off outcome-relevant term
  # if user wants to use joint testing for IV selection
  if(grepl('joint',tolower(MVMR_IV_selection_type),fixed=TRUE)) {
    for(j in 1:m) {
      Th=diag(bxse[j,])%*%Ruu%*%diag(bxse[j,])
      Th=solve(Th)
      v=t(bx[j,])%*%Th%*%bx[j,]
      pj=1-pchisq(v,p); pjs[j]=pj
      if(pj<IV_Pvalue_threshold) keep=c(keep,j)
    }
    # if user wants the union set of significant SNPs as IVs in MVMR
  } else {
    z=bx/bxse
    keep=which(apply(z^2,1,function(h) any(h>qchisq(1-IV_Pvalue_threshold,1))))
  }
  ix=pruning(pjs[keep],LD[keep,keep],LD_pruning_r2)
  ix=keep[ix]
  ### or pruning to achieve certain variance explained by IVs?
} else if(st=='wea') {
  ff=setf(bxunstd,nX,fix_Fstatistic_at)
  ix=ff$ix
  fs=ff$fs
} else {
  stop('simtype must be one of -weak- or -winners-')
}
### return data
mIVs=length(ix)
LD=LD[ix,ix]
bx=bx[ix,];by=by[ix];bxse=bxse[ix,];byse=byse[ix]
bxunstd=bxunstd[ix,]; bxseunstd=bxseunstd[ix,]
IVtype=classIVs(ix,uhpix,chpix)
el=ls()
# vetify weak IV set
# plot(fs);abline(v=which.min(abs(fix_Fstatistic_at-fs)));abline(h=fix_Fstatistic_at)
# h2=colSums(bxunstd^2)
# (meff=nrow(bxunstd))
# (nX-meff-1)/meff*h2/(1-h2)
# abline(h=mean((nX-meff-1)/meff*h2/(1-h2)),col='red')
# keep original objects plus the important ones I made
nn=c('Outcome',paste0('Exposure',1:p))
rownames(RhoME)=colnames(RhoME)=nn
ko=c(sl,'bx','bxse','by','byse','mStart','mIVs','RhoME','LD','theta','IVtype','bxunstd','bxseunstd') 
rm(list=el[el %!in% ko])
