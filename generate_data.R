source('basicfunctions.R')
sl=ls()
n0_xy=floor(proportion_overlapping_in_Xs_and_Y_GWAS*min(c(sample_size_Y,sample_size_Xs)))
nall=sample_size_Y+sample_size_Xs-n0_xy
indY=1:sample_size_Y
indX=(nall-sample_size_Xs+1):nall
G=rbinom(number_of_instruments*nall,2,mafs_of_instruments) # assuming independence for now
G=matrix(G,nr=nall,nc=number_of_instruments)
G=apply(G,2,std)
# order will go UHP, CHP, valid, weak
### CHP IVs
gammaC=rep(0,number_of_instruments)
chpix=1:number_of_CHP_instruments
if(number_of_CHP_instruments>0) {
  gammaC_=runif(number_of_CHP_instruments,0,0.5);
  adj=variance_in_U_explained_by_CHP_instruments/sum(gammaC_^2)
  gammaC_=sqrt(adj)*gammaC_; gammaC[chpix]=-gammaC_
}
### setting U
eU=rnorm(nall,0,sqrt(1-variance_in_U_explained_by_CHP_instruments))
U=G%*%gammaC+eU
### setting X
CorrXX=parthcorr(phenotypic_correlation_Xs,n=number_of_exposures)
GenCorrXX=parthcorr(genetic_correlation_Xs,n=number_of_exposures)
LD=parthcorr(LD_instruments,number_of_instruments)
K=kronecker(GenCorrXX,LD)
B=rmvn(1,rep(0,dim(K)[1]),K) # can effectively add LD the G by adding LD to B
B=matrix(B,nr=number_of_instruments,nc=number_of_exposures)
# need to consider some IVs as weak at this stage
# first k will be strong, remaining m-k will not be
# heritability will be partitioned accordingly
strongix=1:(number_of_instruments-number_of_weak_instruments)
strongh2=variance_in_Xs_explained_by_all_instruments-variance_in_Xs_explained_by_weak_instruments
weakh2=variance_in_Xs_explained_by_all_instruments-strongh2
strongadj=strongh2/colSums(B[strongix,]^2)
weakadj=weakh2/colSums(B[-strongix,]^2)
for(j in 1:ncol(B)) {
  B[strongix,j]=sqrt(strongadj[j])*B[strongix,j]
  B[-strongix,j]=sqrt(weakadj[j])*B[-strongix,j]
}
# effect of confounder
pix=sqrt(variance_in_Xs_explained_by_U)
sdeX=diag(1-variance_in_Xs_explained_by_all_instruments-pix^2,number_of_exposures)
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
gammaU=rep(0,number_of_instruments)
uhpix=(max(chpix)+c(1:number_of_UHP_instruments))
if(number_of_UHP_instruments>0) {
  gammaU_=runif(number_of_CHP_instruments,-1,1);
  adj=variance_in_Y_explained_by_UHP_instruments/sum(gammaU_^2)
  gammaU_=sqrt(adj)*gammaU_; gammaU[uhpix]=gammaU_
}
# INSIDE ASSUMPTION ... ?
vUHPY=variance_in_Y_explained_by_UHP_instruments
vCHPY=variance_in_Y_explained_by_CHP_instruments
eY=rnorm(nall,0,sqrt(1-vXY-vUHPY-vCHPY))
Y=X%*%theta+piy*U+G%*%gammaU+eY
### GWAS
gwas_y=biggwas(Y[indY],G[indY,]) # indY created near the top
by=gwas_y$est; byse=gwas_y$std
bx=bxse=matrix(nr=number_of_instruments,nc=number_of_exposures)
for(j in 1:ncol(bx)) {
  fit=biggwas(X[indX,j],G[indX,])
  bx[,j]=fit$est
  bxse[,j]=fit$std
}
### plot with true values and estimated using true theta
# par(mfrow=c(1,2))
# pfun(B%*%theta,B%*%theta+piy*gammaC+gammaU,chpix,uhpix)
# pfun(bx%*%theta,by,chpix,uhpix)
### standardize?
data=parthstd(bx,by,bxse,byse,mafs_of_instruments,sample_size_Xs,
              sample_size_Y,MR_standardization_type,outcome_type,exposure_types)
### return data
el=ls()
ko=c(sl,'bx','bxse','by','byse') # keep original objects and bx, bxse, by, byse only
rm(list=el[el %!in% ko])
