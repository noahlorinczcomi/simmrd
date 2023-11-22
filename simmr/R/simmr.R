#' A Data-Generating Function
#'
#' This function generates simulated data given a list of parameters
#' @param params List of parameters used to generate simulated data
#' @keywords data
#' @export 
#' @import mvnfast
#' @examples
#' generate()
generate=function(params){
  # assign values in params to local environment
    for(i in 1:length(params)) assign(names(params)[i],params[[i]])
  # make sure user specified reasonable values
    exposure_lengths=c(number_of_exposures,length(Y_variance_explained_by_Xs),length(signs_of_causal_effects))
    if(exposure_lengths[1]<exposure_lengths[2] | exposure_lengths[1]<exposure_lengths[3]) {
      stop("Error in generate(params) :\n",
           "please make sure the parameters for the following values each have the same length:\n",
           "\t- `sample_size_Xs`\n",
           "\t- `Y_variance_explained_by_Xs`\n",
           "\t- `signs_of_causal_effects`\n")
    }
    snp_lengths=c(number_of_causal_SNPs,number_of_UHP_causal_SNPs,number_of_CHP_causal_SNPs)
    if(snp_lengths[1]<(snp_lengths[2]+snp_lengths[3])) {
      stop("Error in generate(params) :\n",
           "please make sure the parameters for the following values are possible:\n",
           "\t- `number_of_causal_SNPs`\n",
           "\t- `number_of_UHP_causal_SNPs`\n",
           "\t- `number_of_CHP_causal_SNPs`\n")
    }
  # generate data
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
    LD=makeBlocks(LD_causal_SNPs,number_of_causal_SNPs,number_of_LD_blocks)
    # LD=parthcorr(LD_causal_SNPs,number_of_causal_SNPs)
    K=kronecker(GenCorrXX,LD)
    Thsq=chol(solve(LD))
    B=mvnfast::rmvn(1,rep(0,dim(K)[1]),K) # can effectively add LD the G by adding LD to B
    B=matrix(B,nr=number_of_causal_SNPs,nc=number_of_exposures)
    if(length(chpix)>0) B[nrow(B):(nrow(B)-length(chpix)),]=0
    th=chol(solve(GenCorrXX))
    cop=pnorm(B%*%th)
    #cor(cop)
    cop=cop%*%chol(GenCorrXX)
    #list(gencor=round(GenCorrXX,2),cop=round(cor(cop),2),corB=round(cor(B),2))
    # rescale to match heritability
    adj=Xs_variance_explained_by_g/colSums(B^2)
    for(i in 1:ncol(B)) B[,i]=sqrt(adj[i])*B[,i]
    # effect of confounder
    pix=sqrt(Xs_variance_explained_by_U)
    sdeX=diag(1-Xs_variance_explained_by_g-pix^2,number_of_exposures)
    sdeX=sqrt(sdeX)
    SigmaEX=sdeX%*%CorrXX%*%sdeX
    eX=mvnfast::rmvn(nall,rep(0,number_of_exposures),SigmaEX)
    X=G%*%B+matrix(pix*U,nr=nall,nc=number_of_exposures)+eX
    ### model for Y
    vXY=Y_variance_explained_by_Xs
    #theta=vXY*signs_of_causal_effects
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
    etaX=X%*%theta
    #if(length(chpix)>0) etaX[which(gammaC!=0)]=0
    Y=etaX+piy*U+G%*%gammaU+eY
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
    ### standardize?
    bx_unstd=bx; bxse_unstd=bxse; by_unstd=by; byse_unstd=byse
    data=parthstd(bx,by,bxse,byse,mafs_of_causal_SNPs,sample_size_Xs,sample_size_Y,MR_standardization)
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
          D_=bxse[j,]
          Th=diag(D_,length(D_))%*%as.matrix(Ruu)%*%diag(D_,length(D_))
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
      ff=setf(bx_unstd,nX,fix_Fstatistic_at)
      ix=ff$ix
      fs=ff$fs
    } else {
      stop('simtype must be one of -weak- or -winners-')
    }
    ### return data
    mIVs=length(ix)
    LD=LD[ix,ix]
    bx=bx[ix,];by=by[ix];bxse=bxse[ix,];byse=byse[ix]
    bx_unstd=bx_unstd[ix,]; bxse_unstd=bxse_unstd[ix,]
    IVtype=classIVs(ix,uhpix,chpix)
    ### did user want to imprecisely estimate LD among the IVs?
    if(N_of_LD_ref<Inf) {LD=rWishart(1,N_of_LD_ref,LD)[,,1]; LD=cov2cor(LD)}
    # verify weak IV set
    # plot(fs);abline(v=which.min(abs(fix_Fstatistic_at-fs)));abline(h=fix_Fstatistic_at)
    # h2=colSums(bx_unstd^2)
    # (meff=nrow(bx_unstd))
    # (nX-meff-1)/meff*h2/(1-h2)
    # abline(h=mean((nX-meff-1)/meff*h2/(1-h2)),col='red')
    nn=c('Outcome',paste0('Exposure',1:p))
    rownames(RhoME)=colnames(RhoME)=nn
    out=list(bx=bx,bxse=bxse,by=by,byse=byse,RhoME=RhoME,LD=LD,theta=theta,IVtype=IVtype,bx_unstd=bx_unstd,bxse_unstd=bxse_unstd,by_unstd=by_unstd,byse_unstd=byse_unstd)
    return(out)
}

#' Helper function
#'
#' Helper function
#' @param x phenotype vector
#' @param G genotype matrix 
#' @keywords 
#' @export
#' @examples
#' biggwas()
biggwas=function(x,G){
  x=as.vector(x)
  ux=mean(x)
  vx=var(x);vx=as.numeric(vx)
  ug=colMeans(G)
  G=t(t(G)-ug)
  vg=colSums(G^2)
  b=(t(G)%*%(x-ux))/vg
  sdb=(vx-b^2*vg/length(x))/length(x)
  A=list()
  A$est=as.vector(b)
  A$std=as.vector(sqrt(sdb))
  return(A)
}

#' Helper function
#'
#' Helper function
#' @param n 
#' @param rho
#' @keywords 
#' @export
#' @examples
#' ar1()
ar1=function(n,rho=0.5) {x=matrix(rho,n,n)^(toeplitz(1:n)-1);diag(x)=1;x}

#' Helper function
#'
#' Helper function
#' @param x
#' @keywords 
#' @export
#' @examples
#' std()
std=function(x) (x-mean(x))/sd(x)

#' Helper function
#'
#' Helper function
#' @param x 
#' @param n
#' @keywords 
#' @export
#' @examples
#' parthcorr()
parthcorr=function(x,n) {
  keys=c('ar','toeplitz','cs')
  boo=sapply(keys,function(h) grepl(h,tolower(x),fixed=TRUE))
  if(boo[1]) return(ar1(n,as.numeric(substr(x,start=5,stop=nchar(x)-1))))
  if(boo[2]) return(1/toeplitz(1:n))
  if(boo[3]) {a_=as.numeric(substr(x,start=4,stop=nchar(x)-1)); return(diag(n)+a_-diag(a_,n))}
  return(diag(n))
}

#' Helper function
#'
#' Helper function
#' @param x 
#' @param y
#' @param chpix
#' @param uhpix
#' @keywords 
#' @export
#' @examples
#' pfun()
pfun=function(x,y,chpix,uhpix,...) {
  plot(x,y,pch=16,col='gray80',...)
  lc=length(chpix)>0
  lu=length(uhpix)>0
  if(lc) points(x[chpix],y[chpix],col='royalblue',pch=16)
  if(lu) points(x[uhpix],y[uhpix],col='indianred',pch=16)
  if(lc & !lu) legend('bottomright','CHP',pch=16,col='royalblue')
  if(!lc & lu) legend('bottomright','UHP',pch=16,col='indianred')
  if(lc & lu) legend('bottomright',c('CHP','UHP'),pch=c(16,16),col=c('royalblue','indianred'))
}

#' Helper function
#'
#' Helper function
#' @param bx 
#' @param by
#' @param bxse
#' @param byse
#' @param maf
#' @param nx
#' @param ny
#' @param MR_standardization_type
#' @keywords 
#' @export
#' @examples
#' parthstd()
parthstd=function(bx,by,bxse,byse,maf,nx,ny,MR_standardization_type) {
  mst=tolower(MR_standardization_type)
  if(mst=='none') return(list(bx=bx,bxse=bxse,by=by,byse=byse))
  if(mst=='z') return(list(bx=bx/bxse,bxse=bx/bx,by=by/byse,byse=by/by))
  # will then do qi & chatterjee if top two didn't return
  by=by/byse/sqrt(ny); byse=1/sqrt(ny)
  bx=bx/bxse/sqrt(nx); bxse=1/sqrt(nx)
  return(list(bx=bx,bxse=bxse,by=by,byse=byse))
}

#' Pruning SNPs
#'
#' Helper function
#' @param jointPs joint p-degree of freedom chi-square tests for IVs
#' @param R LD correlation matrix for SNPs
#' @param r2 upper squared LD r2 threshold for pruning
#' @keywords 
#' @export
#' @examples
#' pruning()
pruning=function(jointPs,R,r2) {
  n=nrow(R);
  ps=cbind(1:n,jointPs); ps=ps[order(ps[,2]),]
  ord=ps[,1]; ps=ps[,2]
  R0=R[ord,ord]; R0[lower.tri(R0,diag=TRUE)]=0
  drop=c()
  for(i in 1:n) {
    if(i %in% drop) next
    vR=R0[i,]
    w=which((vR)^2>r2)
    if(length(w)==0) next else drop=c(drop,w)
  }
  if(sum(drop)>0) keep=c(1:n)[-ord[unique(drop)]] else keep=1:n
  return(keep)
}

#' Helper function
#'
#' Helper function
#' @param ix
#' @param uhpix
#' @param chpix
#' @keywords 
#' @export
#' @examples
#' classIVs()
classIVs=function(ix,uhpix,chpix) {
  keys=c('UHP','CHP')
  ll=list(uhpix,chpix)
  boo=sapply(1:2,function(h) ix %in% ll[[h]])
  boo=matrix(boo,nr=length(ix))
  cl=c();for(i in 1:nrow(boo)) {toa=keys[which(boo[i,])];cl[i]=ifelse(length(toa)==0,'valid',toa)}
  return(cl)
}

#' Helper function
#'
#' Helper function
#' @param bxunstd
#' @param nX
#' @param fix_Fstatistic_at
#' @keywords 
#' @export
#' @examples
#' setf()
setf=function(bxunstd,nX,fix_Fstatistic_at) {
  # currently agnostic to LD structure
  m=nrow(B);p=ncol(bxunstd)
  v=diag(bxunstd%*%t(bxunstd)); 
  ord=cbind(1:m,v);
  ord=ord[order(v,decreasing=FALSE),]
  h2s=a=b=c(); 
  for(i in 1:m) {
    ixi=ord[1:i,1]
    h2s[i]=mean(colSums(matrix(bxunstd[ixi,]^2,nc=p)))
    a[i]=(nX-i-1)/i
    b[i]=h2s[i]/(1-h2s[i])
  }
  fs=a*b
  wf=which.min(abs(fs-fix_Fstatistic_at))
  ix=ord[1:wf,1]
  return(list(ix=ix,fs=fs,h2s=h2s))
}

#' Helper function
#'
#' Helper function
#' @param data direct output from generate()
#' @keywords 
#' @import ggplot2 
#' @export
#' @examples
#' plot_simdata()
plot_simdata=function(data) {
  for(i in 1:length(data)) assign(names(data)[i],data[[i]])
  nX=sample_size_Xs
  cm=colMeans(as.matrix(bxse)); Ruu=diag(cm)%*%RhoME[-1,-1]%*%diag(cm)
  lpse=t(theta)%*%Ruu%*%theta
  lpse=sqrt(lpse)
  pdf=data.frame(lp=c(bx%*%theta),y=c(by),lpse=lpse,yse=c(byse),label=IVtype)
  slp=t(pdf$y)%*%pdf$lp/sum(pdf$lp^2)
  #bxt=chol(solve(LD))%*%as.matrix(bx)%*%chol(solve(Ruu))
  Fs=mean(colSums(bx_unstd^2)); 
  mx=length(c(by))
  Fs=(nX-mx-1)/(mx)*Fs/(1-Fs)
  Fs=paste0('F=',round(Fs))
  annopos=c(min(c(bx)), max(c(by)))*0.75
  if(ncol(bx)==1) xl='SNP-exposure association' else xl='linear predictor for SNP-exposures associations'
  p=ggplot(pdf,aes(lp,y,fill=label)) +
    geom_errorbar(aes(ymin=y-2*yse,ymax=y+2*yse),width=0,color='gray60') +
    geom_errorbarh(aes(xmin=lp-2*lpse,xmax=lp+2*lpse),height=0,color='gray60') +
    geom_vline(xintercept=0,linetype='dashed') +
    geom_hline(yintercept=0,linetype='dashed') +
    geom_abline(intercept=0,slope=1,color='gray80') +
    geom_abline(intercept=0,slope=slp) +
    theme_bw() +
    annotate('text',x=annopos[1],y=annopos[2],label=Fs) +
    labs(x=xl,y='SNP-outcome associations') +
    theme(legend.position='bottom')+
    guides(fill=guide_legend(title='IV type',override.aes=list(size=3))) +
    geom_point(pch=21,size=2)
  p
}

#' A function to make LD blocks
#'
#' Helper function
#' @param LD_causal_SNPs the LD structure of the causal SNPs
#' @param number_of_causal_SNPs the total number of causal SNPs
#' @param nblocks the number of independent LD blocks 
#' @keywords 
#' @import ggplot2 
#' @export
#' @examples
#' makeBlocks()
makeBlocks=function(LD_causal_SNPs,number_of_causal_SNPs,nblocks=1) {
  m=number_of_causal_SNPs
  mat=matrix(0,m,m)
  p=ceiling(m/nblocks)
  blockCor=parthcorr(LD_causal_SNPs,p)
  ixs=ceiling(seq(0,p*nblocks,length.out=nblocks+1))
  for(i in 1:(length(ixs)-1)) {
    ix=(ixs[i]+1):min(ixs[i+1],m)
    mat[ix,ix]=blockCor[1:length(ix),1:length(ix)]
  }
  return(mat)
}
