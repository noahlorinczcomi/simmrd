source('basicfunctions.R')
sl=ls()
n0_xy=floor(proportion_overlapping_in_Xs_and_Y_GWAS*min(c(sample_size_Y,sample_size_Xs)))
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
  gammaC_=runif(number_of_CHP_causal_SNPs,1,5);
  adj=variance_in_U_explained_by_CHP_causal_SNPs/sum(gammaC_^2)
  gammaC_=sqrt(adj)*gammaC_; gammaC[chpix]=-gammaC_
} else {
  chpix=c()
}
### setting U
eU=rnorm(nall,0,sqrt(1-variance_in_U_explained_by_CHP_causal_SNPs))
U=G%*%gammaC+eU
### setting X
CorrXX=parthcorr(phenotypic_correlation_Xs,n=number_of_exposures)
GenCorrXX=parthcorr(genetic_correlation_Xs,n=number_of_exposures)
LD=parthcorr(LD_causal_SNPs,number_of_causal_SNPs)
K=kronecker(GenCorrXX,LD)
B=rmvn(1,rep(5,dim(K)[1]),K) # can effectively add LD the G by adding LD to B
B=matrix(B,nr=number_of_causal_SNPs,nc=number_of_exposures)
# need to consider some IVs as weak at this stage
# first k will be strong, remaining m-k will not be
# heritability will be partitioned accordingly
if(number_of_weak_causal_SNPs==0) variance_in_Xs_explained_by_weak_causal_SNPs=0
strongix=1:(number_of_causal_SNPs-number_of_weak_causal_SNPs)
strongh2=variance_in_Xs_explained_by_all_causal_SNPs-variance_in_Xs_explained_by_weak_causal_SNPs
weakh2=variance_in_Xs_explained_by_all_causal_SNPs-strongh2
strongadj=strongh2/colSums(B[strongix,]^2)
weakadj=weakh2/colSums(B[-strongix,]^2)
for(j in 1:ncol(B)) {
  B[strongix,j]=sqrt(strongadj[j])*B[strongix,j]
  B[-strongix,j]=sqrt(weakadj[j])*B[-strongix,j]
}
# effect of confounder
pix=sqrt(variance_in_Xs_explained_by_U)
sdeX=diag(1-variance_in_Xs_explained_by_all_causal_SNPs-pix^2,number_of_exposures)
sdeX=sqrt(sdeX)
SigmaEX=sdeX%*%CorrXX%*%sdeX
eX=rmvn(nall,rep(0,number_of_exposures),SigmaEX)
X=G%*%B+matrix(pix*U,nr=nall,nc=number_of_exposures)+eX
### model for Y
vXY=variance_in_Y_explained_by_Xs
theta=vXY*rep(1,number_of_exposures)*(length(vXY)==1)+vXY*(length(vXY)>1)
adj=vXY/sum(theta^2)
theta=sqrt(adj)*theta
piy=-sign(theta[1])*sqrt(variance_in_Y_explained_by_U)
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
  adj=variance_in_Y_explained_by_UHP_causal_SNPs/sum(gammaU_^2)
  gammaU_=sqrt(adj)*gammaU_; gammaU[uhpix]=gammaU_
}
vUHPY=variance_in_Y_explained_by_UHP_causal_SNPs
eY=rnorm(nall,0,sqrt(1-vXY-vUHPY-piy^2))
Y=X%*%theta+piy*U+G%*%gammaU+eY
RhoXY=cor(cbind(Y,X))
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
### uncomment if you want to see lots
# par(mfrow=c(1,2))
# pfun(B%*%theta,B%*%theta+piy*gammaC+gammaU,chpix,uhpix,xlab='bx*theta',ylab='by')
# pfun(bx%*%theta,by,chpix,uhpix,xlab='bxhat*theta',ylab='byhat')
### standardize?
data=parthstd(bx,by,bxse,byse,mafs_of_causal_SNPs,sample_size_Xs,
              sample_size_Y,MR_standardization_type,outcome_type,exposure_types)
### instrument selection based on P-value
bx=data$bx;bxse=data$bxse;by=data$by;byse=data$byse
keep=pjs=c()
m=number_of_causal_SNPs; mStart=m; Ruu=RhoME[-1,-1] # take off outcome-relevant term
for(j in 1:m) {
  Th=(bxse[j,]%*%t(bxse[j,]))*Ruu
  Th=solve(Th)
  v=t(bx[j,])%*%Th%*%bx[j,]
  pj=1-pchisq(v,p); pjs[j]=pj
  if(pj<instrument_selection_Pvalue_threshold) keep=c(keep,j)
}
if(length(keep)==0) stop(cat('ERROR: The IV selection P-value threshold does not keep any SNPs. \n This can happen if you simulated a very large number of causal \n SNPs because it implies they each have a very small contribution to heritability'))
mSelected=length(keep)
### any pruning?
ix=pruning(pjs[keep],LD[keep,keep],instrument_selection_LD_pruning_r2)
ix=keep[ix]
mNotPruned=length(ix)
LD=LD[ix,ix]
### return data
el=ls()
# keep original objects plus the important ones I made
nn=c('Outcome',paste0('Exposure',1:p))
rownames(RhoME)=colnames(RhoME)=nn
ko=c(sl,'bx','bxse','by','byse','mStart','mSelected','mNotPruned','RhoME','LD','theta') 
rm(list=el[el %!in% ko])
