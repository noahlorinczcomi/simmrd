#' An Individual-Level Data-Generating Function
#'
#' This function generates simulated individual-level and data given a list of parameters
#' @param params List of parameters used to generate simulated data
#' @keywords data
#' @export 
#' @import mvnfast
#' @examples
#' \dontrun{
#' individual_params <-
#'   list(
#'     sample_size_Xs = 5e4, # exposure GWAS sample sizes
#'     sample_size_Y = 5e4, # outcome GWAS sample size
#'     prop_gwas_overlap_Xs_and_Y = 0.5, # proportion of exposure and outcome GWAS overlap
#'     number_of_exposures = 2, # 4 number of exposures
#'     phenotypic_correlation_Xs = 0.2, # phenotypic correlation between exposures
#'     genetic_correlation_Xs = 0, # genetic correlation between exposures
#'     Xs_variance_explained_by_U = 1/4 - 0.12, # exposures variance explained by confounder
#'     Y_variance_explained_by_Xs = c(0, 0.5), # outcome variance explained by exposures
#'     signs_of_causal_effects = c(1, 1), # signs of causal effects
#'     Y_variance_explained_by_U = 0.1, # outcome variance explained by confounder
#'     number_of_causal_SNPs = 200, # number of SNPs causing exposures
#'     mafs_of_causal_SNPs = stats::runif(100, 0.1, 0.5), # minor allele frequency of causal SNPs
#'     Xs_variance_explained_by_g = 0.12, # exposures variance explained by SNPs 
#'     number_of_UHP_causal_SNPs = 30, # number of UHP exposure SNPs
#'     number_of_CHP_causal_SNPs = 10, # number of CHP exposure SNPs
#'     Y_variance_explained_by_UHP = 0.05, # outcome variance explained by UHP SNPs
#'     U_variance_explained_by_CHP = 0.05, # outcome variance explained by CHP SNPs
#'     LD_causal_SNPs = 'I', # independent causal exposure SNPs
#'     number_of_LD_blocks = 1, # number of independent LD blocks
#'     MR_standardization = 'Z', # standardization of GWAS summary statistics 
#'     simtype = 'weak', # simulation performed using weak instruments
#'     MVMR_IV_selection_type = 'joint', # P-values for IV selection based on joint test for exposures
#'     IV_Pvalue_threshold = 1, # P-value threshold for candidate IVs
#'     LD_pruning_r2 = 1, # upper boundary of squared LD correlation
#'     N_of_LD_ref = Inf, # size of the LD reference panel
#'     fix_Fstatistic_at = 10 # average across exposures using full MVMR IV set
#'   )
#' gwas_data <- generate_individual(individual_params)
#' }
generate_individual=function(params){
  # assign values in params to local environment
    for(i in 1:length(params)) assign(names(params)[i],params[[i]])
  # make sure user specified reasonable values
    exposure_lengths=c(number_of_exposures,length(Y_variance_explained_by_Xs),length(signs_of_causal_effects))
    if(exposure_lengths[1]<exposure_lengths[2] | exposure_lengths[1]<exposure_lengths[3]) {
      stop("please make sure the parameters for the following values each have the same length:\n",
           "\t- `sample_size_Xs`\n",
           "\t- `Y_variance_explained_by_Xs`\n",
           "\t- `signs_of_causal_effects`\n")
    }
    snp_lengths=c(number_of_causal_SNPs,number_of_UHP_causal_SNPs,number_of_CHP_causal_SNPs)
    if(snp_lengths[1]<(snp_lengths[2]+snp_lengths[3])) {
      stop("please make sure the parameters for the following values are possible:\n",
           "\t- `number_of_causal_SNPs`\n",
           "\t- `number_of_UHP_causal_SNPs`\n",
           "\t- `number_of_CHP_causal_SNPs`\n")
    }
    if(length(Y_variance_explained_by_Xs)==1) Y_variance_explained_by_Xs=rep(Y_variance_explained_by_Xs,number_of_exposures)
    if(length(prop_gwas_overlap_Xs_and_Y)==1) prop_gwas_overlap_Xs_and_Y=rep(prop_gwas_overlap_Xs_and_Y,number_of_exposures)
  # generate data
    n0_xy=floor(prop_gwas_overlap_Xs_and_Y*min(c(sample_size_Y,sample_size_Xs)))
    nall=sample_size_Y+sample_size_Xs-n0_xy
    indY=1:sample_size_Y
    indX=(nall-sample_size_Xs+1):nall
    mafs_of_causal_SNPs=0.3 
    G=stats::rbinom(number_of_causal_SNPs*nall,2,mafs_of_causal_SNPs) # assuming independence for now
    G=matrix(G,nrow=nall,ncol=number_of_causal_SNPs)
    G=apply(G,2,std)
    # order will go UHP, CHP, valid, weak
    ### CHP IVs
    gammaC=rep(0,number_of_causal_SNPs)
    if(number_of_CHP_causal_SNPs>0) {
    chpix=1:number_of_CHP_causal_SNPs
    gammaC_=stats::runif(number_of_CHP_causal_SNPs,-1/2,1/2);
    adj=U_variance_explained_by_CHP/sum(gammaC_^2)
    gammaC_=sqrt(adj)*gammaC_; gammaC[chpix]=-gammaC_
    } else {
    chpix=c()
    }
    ### setting U
    eU=stats::rnorm(nall,0,sqrt(1-U_variance_explained_by_CHP))
    U=G%*%gammaC+eU
    ### setting X
    CorrXX=parthcorr(phenotypic_correlation_Xs,n=number_of_exposures)
    GenCorrXX=parthcorr(genetic_correlation_Xs,n=number_of_exposures)
    LD=makeBlocks(LD_causal_SNPs,number_of_causal_SNPs,number_of_LD_blocks)
    # LD=parthcorr(LD_causal_SNPs,number_of_causal_SNPs)
    K=kronecker(GenCorrXX,LD)
    Thsq=chol(solve(LD))
    B=mvnfast::rmvn(1,rep(0,dim(K)[1]),K) # can effectively add LD the G by adding LD to B
    B=matrix(B,nrow=number_of_causal_SNPs,ncol=number_of_exposures)
    if(length(chpix)>0) B[nrow(B):(nrow(B)-length(chpix)),]=0
    th=chol(solve(GenCorrXX))
    cop=stats::pnorm(B%*%th)
    #stats::cor(cop)
    cop=cop%*%chol(GenCorrXX)
    #list(gencor=round(GenCorrXX,2),cop=round(stats::cor(cop),2),corB=round(stats::cor(B),2))
    # rescale to match heritability
    adj=Xs_variance_explained_by_g/colSums(B^2)
    for(i in 1:ncol(B)) B[,i]=sqrt(adj[i])*B[,i]
    # effect of confounder
    pix=sqrt(Xs_variance_explained_by_U)
    sdeX=diag(1-Xs_variance_explained_by_g-pix^2,number_of_exposures)
    sdeX=sqrt(sdeX)
    SigmaEX=sdeX%*%CorrXX%*%sdeX
    eX=mvnfast::rmvn(nall,rep(0,number_of_exposures),SigmaEX)
    X=G%*%B+matrix(pix*U,nrow=nall,ncol=number_of_exposures)+eX
    ### model for Y
    vXY=Y_variance_explained_by_Xs
    #theta=vXY*signs_of_causal_effects
    theta=vXY*rep(1,number_of_exposures)*(length(vXY)==1)+vXY*(length(vXY)>1)
    theta=theta*signs_of_causal_effects
    adj=vXY/sum(theta^2)
    adj[is.nan(adj)]=0
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
    gammaU_=stats::runif(number_of_UHP_causal_SNPs,-1,1);
    adj=Y_variance_explained_by_UHP/sum(gammaU_^2)
    gammaU_=sqrt(adj)*gammaU_; gammaU[uhpix]=gammaU_
    }
    IVtype=rep('valid',number_of_causal_SNPs)
    if(length(uhpix)>0) IVtype[uhpix]='UHP'
    if(length(chpix)>0) IVtype[chpix]='CHP'
    vUHPY=Y_variance_explained_by_UHP
    eY=stats::rnorm(nall,0,sqrt(1-vXY-vUHPY-piy^2))
    etaX=X%*%theta
    #if(length(chpix)>0) etaX[which(gammaC!=0)]=0
    Y=etaX+piy*U+G%*%gammaU+eY
    cyx=c(stats::cov(X,Y))
    rho2=t(cyx)%*%solve(stats::cov(X))%*%cyx
    RhoXY=stats::cor(cbind(Y,X))
    Sinv=solve(RhoXY[2:ncol(RhoXY),2:ncol(RhoXY)])
    p=ncol(B);nY=length(indY); nX=sample_size_Xs
    nn=sqrt(c(nY,rep(nX,p)));nn=nn%*%t(nn) # for RhoME (correlations between measurement errors)
    no=nn*0+n0_xy;no[1,1]=nY;no[2:(p+1),2:(p+1)]=nX
    RhoME=no/nn*RhoXY
    ### GWAS
    gwas_y=biggwas(Y[indY],G[indY,]) # indY created near the top
    by=gwas_y$est; byse=gwas_y$std
    bx=bxse=matrix(nrow=number_of_causal_SNPs,ncol=number_of_exposures)
    for(j in 1:ncol(bx)) {
    fit=biggwas(X[indX,j],G[indX,])
    bx[,j]=fit$est
    bxse[,j]=fit$std
    }
    ### standardize?
    bx_unstd=bx; bxse_unstd=bxse; by_unstd=by; byse_unstd=byse
    data=parthstd(bx,by,bxse,byse,mafs_of_causal_SNPs,sample_size_Xs,sample_size_Y,MR_standardization)
    bx=data$bx;bxse=data$bxse;by=data$by;byse=data$byse; m=nrow(bx)
    ### saving GWAS estimates from all SNPs
    bxunstd_all=bx
    bxseunstd_all=bxse
    byunstd_all=by
    byseunstd_all=byse
    ### did user want to imprecisely estimate LD among the IVs?
    R=LD
    R0=R
    if(N_of_LD_ref<Inf) R=stats::rWishart(1,N_of_LD_ref,R)[,,1]
    if(length(c(R))>1) R=stats::cov2cor(R) else R=R/N_of_LD_ref
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
          pj=1-stats::pchisq(v,p); pjs[j]=pj
          if(pj<IV_Pvalue_threshold) keep=c(keep,j)
        }
        # if user wants the union set of significant SNPs as IVs in MVMR
    } else {
        z=bx/bxse
        keep=which(apply(z^2,1,function(h) any(h>stats::qchisq(1-IV_Pvalue_threshold,1))))
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
    R=as.matrix(R[ix,ix])
    R0=as.matrix(R0[ix,ix])
    IVtype=IVtype[ix]
    bx=as.matrix(bx[ix,]);bxse=as.matrix(bxse[ix,])
    by=by[ix];byse=byse[ix]
    bx_unstd=as.matrix(bx_unstd[ix,]); bxse_unstd=as.matrix(bxse_unstd[ix,])
    by_unstd=by_unstd[ix]; byse_unstd=byse_unstd[ix]
    # if only have one IV selected
    if(length(by)==1) {bx=t(bx);bxse=t(bxse);bx_unstd=t(bx_unstd);bxse_unstd=t(bxse_unstd)} 
    # verify weak IV set
    # plot(fs);abline(v=which.min(abs(fix_Fstatistic_at-fs)));abline(h=fix_Fstatistic_at)
    # h2=colSums(bx_unstd^2)
    # (meff=nrow(bx_unstd))
    # (nX-meff-1)/meff*h2/(1-h2)
    # abline(h=mean((nX-meff-1)/meff*h2/(1-h2)),col='red')
    nn=c('Outcome',paste0('Exposure',1:p))
    rownames(RhoME)=colnames(RhoME)=nn
    out=list(bx=bx,bxse=bxse,by=by,byse=byse,RhoME=RhoME,LDMatrix=R0,LDhatMatrix=R,theta=theta,IVtype=IVtype,bx_unstd=bx_unstd,bxse_unstd=bxse_unstd,by_unstd=by_unstd,byse_unstd=byse_unstd)
    return(out)
}

#' An Summary Statistic Data-Generating Function
#'
#' This function generates simulated summary-level and data given a list of parameters
#' @param params List of parameters used to generate simulated data
#' @keywords data
#' @export 
#' @import mvnfast
#' @examples
#' \dontrun{
#' summary_params <-
#'   list(
#'     sample_size_Xs = 30000, # exposure GWAS sample sizes
#'     sample_size_Y = 30000, # outcome GWAS sample size
#'     prop_gwas_overlap_Xs_and_Y = 1, # proportion of exposures' and outcome GWAS overlap
#'     number_of_exposures = 3, # number of exposures
#'     number_of_causal_SNPs = 100, # number of SNPs causing each exposure
#'     number_of_UHP_causal_SNPs = 0, # number of UHP causal SNPs
#'     number_of_CHP_causal_SNPs = 20, # number of CHP causal SNPs
#'     ratio_of_UHP_variance = 0.15, # ratio of UHP variance to valid IV variance
#'     ratio_of_CHP_variance = 0.25, # ratio of CHP variance to valid IV variance
#'     CHP_correlation = -0.5, # correlation between CHP and valid IV effect sizes
#'     simtype = 'winners', # performs IV selection based on P-value
#'     fix_Fstatistic_at = 10, # ignored because simtype='winners'
#'     prop_gwas_overlap_Xs = 1, # overlap of exposures' GWAS
#'     phenotypic_correlation_Xs = 0.3, # phenotypic correlations between exposures
#'     genetic_correlation_Xs = 0.15,  # genetic correlation between exposures
#'     phenotypic_correlations_Xs_and_Y = 0.3, # phenotypic correlations b/w exposures and outcome
#'     true_causal_effects = 0.3, # true causal effect sizes
#'     Xs_variance_explained_by_g = 0.10, # exposure variance explained by SNPs
#'     LD_causal_SNPs = 'ar1(0.5)', # LD between causal exposure SNPs
#'     number_of_LD_blocks = 3, # number of independent LD blocks
#'     MR_standardization = 'none', # does not standardize GWAS estimates
#'     MVMR_IV_selection_type = 'union', # SNPs associated with >0 exposures are candidate IVs
#'     IV_Pvalue_threshold = 5e-8, # only SNPs with P<this threshold are candidate IVs
#'     LD_pruning_r2 = 1, # the upper LD r2 pruning threshold for candidate IVs
#'     N_of_LD_ref = Inf # the sample size of the LD reference panel
#'   )
#' gwas_data <- generate_summary(summary_params)
#' }
generate_summary=function(params) {
  # assign values in params to local environment
  for(i in 1:length(params)) assign(names(params)[i],params[[i]])
  ### checks to make sure input makes sense
  snp_lengths=c(number_of_causal_SNPs,number_of_UHP_causal_SNPs,number_of_CHP_causal_SNPs)
  if(snp_lengths[1]<(snp_lengths[2]+snp_lengths[3])) {
    stop("please make sure the parameters for the following values are possible:\n",
         "\t- `number_of_causal_SNPs`\n",
         "\t- `number_of_UHP_causal_SNPs`\n",
         "\t- `number_of_CHP_causal_SNPs`\n")
  }
  p=number_of_exposures
  SigmaBB=parthcorr(genetic_correlation_Xs,p)
  if(length(Xs_variance_explained_by_g)==1) D=sqrt(diag(rep(Xs_variance_explained_by_g,p),p))
  SigmaBB=D%*%SigmaBB%*%D
  R=parthcorr(LD_causal_SNPs,number_of_causal_SNPs)
  Rxx=parthcorr(phenotypic_correlation_Xs,p)
  if(length(phenotypic_correlations_Xs_and_Y)==1) rxy=rep(phenotypic_correlations_Xs_and_Y,p) else rxy=phenotypic_correlations_Xs_and_Y
  if(length(prop_gwas_overlap_Xs_and_Y)==1) Propxy=rep(prop_gwas_overlap_Xs_and_Y,p) else Propxy=prop_gwas_overlap_Xs_and_Y
  Propxx=parthcorr(prop_gwas_overlap_Xs,p)
  if(length(sample_size_Xs)==1) Nxx=rep(sample_size_Xs,p) else Nxx=sample_size_Xs
  Prop=diag(p+1)
  Prop[2:(p+1),1]=Prop[1,2:(p+1)]=Propxy
  Prop[2:(p+1),2:(p+1)]=Propxx
  ns=c(sample_size_Y,Nxx)
  N=ns%*%t(ns) # denominator
  Ni=matrix(0,p+1,p+1) # numerator
  for(i in 1:(p+1)) {
    for(j in 1:(p+1)) {
      Ni[i,j]=Prop[i,j]*min(ns[i],ns[j])
    }
  }
  Ryx=diag(1,p+1) # constants
  Ryx[2:(p+1),1]=Ryx[1,2:(p+1)]=rxy
  Ryx[2:(p+1),2:(p+1)]=Rxx
  Pwu=Ryx*Ni/N # covariance between GWAS estimation errors for (outcome, exposures)
  wu=rmvn(number_of_causal_SNPs,rep(0,ncol(Pwu)),Pwu) # N(0,I,Pwu)
  if(toupper(LD_causal_SNPs)=='I') { # if iid SNPs
    B=rmvn(number_of_causal_SNPs,rep(0,p),SigmaBB/number_of_causal_SNPs)
  } else { # if SNPs are correlated
    Th=kronecker(SigmaBB/number_of_causal_SNPs,R)
    B=rmvn(number_of_causal_SNPs*p,rep(0,number_of_causal_SNPs*p),Th)
    B=matrix(B,number_of_causal_SNPs,p)
  }
  ## do not use copula method on B
  # linear predictor and true variance explained
  if(length(true_causal_effects)<p) true_causal_effects=rep(true_causal_effects,p)
  a=B%*%true_causal_effects
  D=diag(1-Xs_variance_explained_by_g,p)
  rho2=t(true_causal_effects)%*%(t(B)%*%R%*%B+D)%*%true_causal_effects
  # CHP
  IVtype=rep('valid',number_of_causal_SNPs)
  gammaU=gammaC=rep(0,number_of_causal_SNPs)
  lp_var=stats::var(B%*%c(true_causal_effects)+as.matrix(wu[,-1])%*%c(true_causal_effects)+wu[,1])
  if(ratio_of_CHP_variance>0 & number_of_CHP_causal_SNPs>0) {
    ix=1:number_of_CHP_causal_SNPs
    IVtype[ix]='CHP'
    gammaC[ix]=as.matrix(B[ix,])%*%c(true_causal_effects)*(-1+CHP_correlation)+stats::rnorm(length(ix),0,sqrt(ratio_of_CHP_variance*lp_var))
  }
  # UHP
  if(ratio_of_UHP_variance>0 & number_of_UHP_causal_SNPs>0) {
    ix0=match('valid',IVtype)+1
    ix=ix0:(ix0+number_of_UHP_causal_SNPs-1)
    IVtype[ix]='UHP'
    gammaU[ix]=stats::rnorm(length(ix),0,sqrt(ratio_of_UHP_variance*lp_var))
  }
  # Bhat, ahat
  valid_ix=which(IVtype=='valid')
  a=a+gammaU+gammaC
  by=a+wu[,1]
  bx=B+wu[,-1]
  bxse=bx*0; for(i in 1:p) bxse[,i]=stats::rchisq(number_of_causal_SNPs,Nxx[i]-1)/(Nxx[i]-1)/Nxx[i] # sigma2=1/Nxx[i]
  bxse=sqrt(bxse)
  byse=stats::rchisq(number_of_causal_SNPs,sample_size_Y-1)/(sample_size_Y-1)/sample_size_Y
  byse=sqrt(byse)
  # cols=rep('black',number_of_causal_SNPs)
  # cols[IVtype=='UHP']='red'
  # cols[IVtype=='CHP']='blue'
  # plot(bhat%*%true_causal_effects,ahat,pch=3,col=cols)
  ### IV selection
  ### standardize?
  bx_unstd=bx; bxse_unstd=bxse; by_unstd=by; byse_unstd=byse
  mafs_of_causal_SNPs=0.3
  data=parthstd(bx,by,bxse,byse,mafs_of_causal_SNPs,sample_size_Xs,sample_size_Y,MR_standardization)
  bx=data$bx;bxse=data$bxse;by=data$by;byse=data$byse; m=nrow(bx)
  ### saving GWAS estimates from all SNPs
  bxunstd_all=bx
  bxseunstd_all=bxse
  byunstd_all=by
  byseunstd_all=byse
  ### measurement error correlation matrix
  RhoME=stats::cov2cor(Pwu)
  ### did user want to imprecisely estimate LD among the IVs?
  R0=R
  if(N_of_LD_ref<Inf) R=stats::rWishart(1,N_of_LD_ref,R)[,,1]
  if(length(c(R))>1) R=stats::cov2cor(R) else R=R/N_of_LD_ref
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
        pj=stats::pchisq(v,p,lower.tail=FALSE); pjs[j]=pj
        if(pj<IV_Pvalue_threshold) keep=c(keep,j)
      }
      # if user wants the union set of significant SNPs as IVs in MVMR
    } else {
      z=bx/bxse
      keep=which(apply(z^2,1,function(h) any(h>stats::qchisq(1-IV_Pvalue_threshold,1))))
      pjs=apply(z^2,1,function(h)stats::pchisq(h,1,lower.tail=FALSE))
    }
    if(length(keep)==0) stop('no IVs were selected given your selection criteria')
    ix=pruning(pjs[keep],R[keep,keep],LD_pruning_r2)
    ix=keep[ix]
    ### or pruning to achieve certain variance explained by IVs?
  } else if(st=='wea') {
    ff=setf(bx_unstd,Nxx,fix_Fstatistic_at)
    ix=ff$ix
    fs=ff$fs
  } else {
    stop('simtype must be one of -weak- or -winners-')
  }
  ### return data
  mIVs=length(ix)
  R=as.matrix(R[ix,ix])
  R0=as.matrix(R0[ix,ix])
  IVtype=IVtype[ix]
  bx=as.matrix(bx[ix,]);bxse=as.matrix(bxse[ix,])
  by=by[ix];byse=byse[ix]
  bx_unstd=as.matrix(bx_unstd[ix,]); bxse_unstd=as.matrix(bxse_unstd[ix,])
  by_unstd=by_unstd[ix]; byse_unstd=byse_unstd[ix]
  # if only one IV selected
  if(mIVs==1) {bx=t(bx);bxse=t(bxse);bx_unstd=t(bx_unstd);bxse_unstd=t(bxse_unstd)} 
  # verify weak IV set
  # plot(fs);abline(v=which.min(abs(fix_Fstatistic_at-fs)));abline(h=fix_Fstatistic_at)
  # h2=colSums(bx_unstd^2)
  # (meff=nrow(bx_unstd))
  # (nX-meff-1)/meff*h2/(1-h2)
  # abline(h=mean((nX-meff-1)/meff*h2/(1-h2)),col='red')
  nn=c('Outcome',paste0('Exposure',1:p))
  rownames(RhoME)=colnames(RhoME)=nn
  out=list(bx=bx,bxse=bxse,by=by,byse=byse,RhoME=RhoME,LDMatrix=R0,LDhatMatrix=R,theta=true_causal_effects,IVtype=IVtype,bx_unstd=bx_unstd,bxse_unstd=bxse_unstd,by_unstd=by_unstd,byse_unstd=byse_unstd)
  return(out)
}

#' Helper function
#'
#' Helper function
#' @param exposure_overlap_proportions scalar or matrix of overlap proportions between exposures GWAS
#' @param prop_gwas_overlap_Xs_and_Y scalar or vector of overlap proportions between exposures and outcome GWAS
#' @param number_of_exposures number of exposures
#' @export
#' @examples
#' adj_overlap(
#'   exposure_overlap_proportions = 0.2,
#'   prop_gwas_overlap_Xs_and_Y = 0.1,
#'   number_of_exposures = 3
#' )
adj_overlap=function(exposure_overlap_proportions,prop_gwas_overlap_Xs_and_Y,number_of_exposures){
  # limit exposure_overlap_proportions (some values will not be possible given other overlap parameters)
  # this function will very likely only be used with generate_summary() since it may be near impossible
  # to generate an arbitrary overlap structure of individual level data
  ep=exposure_overlap_proportions
  ep=prop_gwas_overlap_Xs_and_Y*exposure_overlap_proportions/number_of_exposures
  p=number_of_exposures
  Prop=diag(p+1)
  Prop[2:(p+1),1]=Prop[1,2:(p+1)]=prop_gwas_overlap_Xs_and_Y
  Prop[2:(p+1),2:(p+1)]=ep
  diag(Prop)=1
  cn=c('Outcome',paste0('Exposure',1:p))
  colnames(Prop)=rownames(Prop)=cn
  return(Prop)
}


#' Helper function
#'
#' Helper function
#' @param x phenotype vector
#' @param G genotype matrix 
#' @export
#' @examples
#' \dontrun{
#' biggwas()
#' }
biggwas=function(x,G){
  x=as.vector(x)
  ux=mean(x)
  vx=stats::var(x);vx=as.numeric(vx)
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
#' @param n The number of rows (and columns) of the matrix
#' @param rho rho
#' @export
#' @examples
#' ar1(2)
ar1=function(n,rho=0.5) rho^stats::toeplitz(0:(n-1))

#' Helper function
#'
#' Helper function
#' @param x Vector to standardise
#' @export
#' @examples
#' std(0:10)
std=function(x) (x-mean(x))/stats::sd(x)

#' Helper function
#'
#' Helper function
#' @param x Matrix or numeric value
#' @param n Number of rows and columns of the matrix
#' @export
#' @examples
#' \dontrun{
#' parthcorr()
#' }
parthcorr=function(x,n) {
  if(is.matrix(x)) return(x)
  if(is.numeric(x)) {M=matrix(x,n,n);diag(M)=1;return(M)}
  keys=c('ar','toeplitz','cs')
  boo=sapply(keys,function(h) grepl(h,tolower(x),fixed=TRUE))
  if(boo[1]) return(ar1(n,as.numeric(substr(x,start=5,stop=nchar(x)-1))))
  if(boo[2]) return(1/stats::toeplitz(1:n))
  if(boo[3]) {a_=as.numeric(substr(x,start=4,stop=nchar(x)-1)); return(diag(n)+a_-diag(a_,n))}
  return(diag(n))
}

#' Helper function
#'
#' Helper function
#' @param x x
#' @param y y
#' @param chpix chipx
#' @param uhpix uhpix
#' @param ... Additional arguments passed to \code{plot()}
#' @export
#' @examples
#' \dontrun{
#' pfun()
#' }
pfun=function(x,y,chpix,uhpix,...) {
  plot(x,y,pch=16,col='gray80',...)
  lc=length(chpix)>0
  lu=length(uhpix)>0
  if(lc) graphics::points(x[chpix],y[chpix],col='royalblue',pch=16)
  if(lu) graphics::points(x[uhpix],y[uhpix],col='indianred',pch=16)
  if(lc & !lu) graphics::legend('bottomright','CHP',pch=16,col='royalblue')
  if(!lc & lu) graphics::legend('bottomright','UHP',pch=16,col='indianred')
  if(lc & lu) graphics::legend('bottomright',c('CHP','UHP'),pch=c(16,16),col=c('royalblue','indianred'))
}

#' Helper function
#'
#' Helper function
#' @param bx bx
#' @param by by
#' @param bxse bxse
#' @param byse byse
#' @param maf Minor allele frequency
#' @param nx nx
#' @param ny ny
#' @param MR_standardization_type Standardization type
#' @export
#' @examples
#' \dontrun{
#' parthstd()
#' }
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
#' Pruning SNPs
#' @param jointPs joint p-degree of freedom chi-square tests for IVs
#' @param R LD correlation matrix for SNPs
#' @param r2 upper squared LD r2 threshold for pruning
#' @export
#' @examples
#' \dontrun{
#' pruning()
#' }
pruning=function(jointPs,R,r2) {
  R=as.matrix(R)
  n=nrow(R);
  ps=cbind(1:n,jointPs)
  ps=ps[order(ps[,2]),]
  ps=matrix(ps,nrow=n,ncol=2)
  ord=ps[,1]; ps=ps[,2]
  R0=R[ord,ord]; 
  R0=as.matrix(R0)
  R0[lower.tri(R0,diag=TRUE)]=0
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
#' @param ix ix
#' @param uhpix uhpix
#' @param chpix chpix
#' @export
#' @examples
#' \dontrun{
#' classIVs()
#' }
classIVs=function(ix,uhpix,chpix) {
  keys=c('UHP','CHP')
  ll=list(uhpix,chpix)
  boo=sapply(1:2,function(h) ix %in% ll[[h]])
  boo=matrix(boo,nrow=length(ix))
  cl=c();for(i in 1:nrow(boo)) {toa=keys[which(boo[i,])];cl[i]=ifelse(length(toa)==0,'valid',toa)}
  return(cl)
}

#' Helper function
#'
#' Helper function
#' @param bxunstd bxunstd
#' @param nX nX
#' @param fix_Fstatistic_at Value to fix the F-statistic at
#' @export
#' @examples
#' # setf()
setf=function(bxunstd,nX,fix_Fstatistic_at) {
  # currently agnostic to LD structure
  # keeps only the weakest IVs
  bxunstd=as.matrix(bxunstd)
  m=nrow(bxunstd);p=ncol(bxunstd)
  v=diag(bxunstd%*%t(bxunstd)); 
  ord=cbind(1:m,v);
  ord=ord[order(v,decreasing=FALSE),]
  h2s=a=b=c(); 
  for(i in 1:m) {
    ixi=ord[1:i,1]
    h2s[i]=mean(colSums(matrix(bxunstd[ixi,]^2,ncol=p)))
    a[i]=(stats::median(nX)-i-1)/i
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
#' @import ggplot2 
#' @export
#' @examples
#' # plot_simdata_lower()
plot_simdata_lower=function(data,params=params,showFstat=TRUE) {
  for(i in 1:length(data)) assign(names(data)[i],data[[i]])
  nX=params$sample_size_Xs
  bx=as.matrix(bx)
  p=ncol(as.matrix(bx))
  cm=colMeans(as.matrix(bxse)); Ruu=diag(cm,p)%*%RhoME[-1,-1]%*%diag(cm,p)
  lpse=t(theta)%*%Ruu%*%theta
  lpse=sqrt(lpse)
  pdf=data.frame(lp=c(as.matrix(bx)%*%c(theta)),y=c(by),lpse=c(lpse),yse=c(byse),label=IVtype)
  # slope using valid IVs
  valid_ix=which(data$IVtype=='valid')
  validslp=t(pdf$y[valid_ix])%*%pdf$lp[valid_ix]/sum(pdf$lp[valid_ix]^2)
  # slope using all IVS
  invalidslp=t(pdf$y)%*%pdf$lp/sum(pdf$lp^2)
  #bx_unstd=chol(solve(LD))%*%as.matrix(bx_unstd)%*%chol(solve(Ruu))
  Fs=mean(colSums(as.matrix(bx_unstd)^2)); 
  mx=length(c(by))
  Fs=(nX-mx-1)/(mx)*Fs/(1-Fs)
  Fs=paste0('F=',round(Fs))
  if(p>1) Fs=paste0('Mean MVMR ',Fs) else Fs=paste0('UVMR ',Fs)
  annopos=c(min(c(bx)), max(c(by)))*0.75
  if(ncol(bx)==1) xl='SNP-exposure association' else xl='linear predictor for \n SNP-exposures associations'
  p=ggplot(pdf,aes(lp,y,fill=label)) +
    geom_errorbar(aes(ymin=y-2*yse,ymax=y+2*yse,color=label),width=0,color='gray60') +
    geom_errorbarh(aes(xmin=lp-2*lpse,xmax=lp+2*lpse,color=label),height=0,color='gray60') +
    geom_vline(xintercept=0,linetype='dashed') +
    geom_hline(yintercept=0,linetype='dashed') +
    geom_abline(intercept=0,slope=invalidslp,color='gray70',lwd=1) +
    geom_abline(intercept=0,slope=validslp,color='royalblue',lwd=1) +
    theme_bw() +
    labs(x=xl,y='SNP-outcome associations') +
    theme(legend.position='bottom')+
    guides(fill=guide_legend(title='IV type',override.aes=list(size=3))) +
    geom_point(pch=21,size=2,alpha=0.75)
  if(showFstat) p=p + labs(subtitle=Fs)
  p
}

#' Helper function
#'
#' Helper function
#' @param data direct output from generate()
#' @import ggplot2 
#' @export
#' @examples
#' # plot_simdata()
plot_simdata=function(data,params=params,exposure_specific_plot='total',verbose=TRUE) {
  p=ncol(as.matrix(data$bx)) # number of exposures
  # if p>1, plot linear predictor
  p0=plot_simdata_lower(data,params)
  if(p>1) p0=p0+ggtitle('All exposures') # linear predictor
  ll=list()
  ll[[1]]=p0
  if(p>1) {
    # do I only want to plot their total causal effects?
    counter=1 # counting the number of plots I am making (and not skipping)
    for(i in 1:p) {
      data0=data
      doskip=FALSE
      if(exposure_specific_plot=='total' | !(exposure_specific_plot%in%c('joint','conditional'))) {
        total=TRUE
        bxi=data0$bx[,i]
        bxsei=data0$bxse[,i]
        ixi=which((bxi/bxsei)^2>stats::qchisq(1-params$IV_Pvalue_threshold,1))
        if(length(ixi)==0) {
          if(verbose) cat('Exposure ', i, ' has no significant IVs\n',sep='')
          doskip=TRUE
        }
      } else {
        total=FALSE
        ixi=1:nrow(as.matrix(data0$bx))
      }
      if(doskip) next
      counter=counter+1
      data0$bx=as.matrix(data$bx[ixi,i])
      data0$bxse=as.matrix(data$bxse[ixi,i])
      data0$RhoME=as.matrix(data$RhoME[c(1,i+1),c(1,i+1)])
      data0$theta=ifelse(total,1,data$theta[i]) # if doing total causal effect, ignore theta
      data0$by=data0$by[ixi]
      data0$byse=data0$byse[ixi]
      data0$IVtype=data0$IVtype[ixi]
      obj=plot_simdata_lower(data0,params,showFstat=TRUE)+ggplot2::ggtitle(paste0('Exposure ',i))
      ll[[i+1]]=obj
    }
    ix1=floor(counter/2)
    ix2=ceiling(counter/ix1)
    ll=Filter(function(x) !is.null(x), ll) # remove NULL values if they are present (happens if an exposure has no IVs)
    allp=ggpubr::ggarrange(plotlist=ll,nrow=ix1,ncol=ix2,common.legend=TRUE)
    return(allp)
  } else {
    # if p==1, just plot
    return(p0)
  }
}

#' A function to make LD blocks
#'
#' Helper function
#' @param LD_causal_SNPs the LD structure of the causal SNPs
#' @param number_of_causal_SNPs the total number of causal SNPs
#' @param nblocks the number of independent LD blocks 
#' @import ggplot2 
#' @export
#' @examples
#' # makeBlocks()
makeBlocks=function(LD_causal_SNPs,number_of_causal_SNPs,nblocks=1) {
  m=number_of_causal_SNPs
  mat=matrix(0,m,m)
  p=ceiling(m/nblocks)
  blockCor=parthcorr(LD_causal_SNPs,p)
  if(dim(blockCor)[1]==m) return(blockCor)
  ixs=ceiling(seq(0,p*nblocks,length.out=nblocks+1))
  for(i in 1:(length(ixs)-1)) {
    ix=(ixs[i]+1):min(ixs[i+1],m)
    mat[ix,ix]=blockCor[1:length(ix),1:length(ix)]
  }
  return(mat)
}
