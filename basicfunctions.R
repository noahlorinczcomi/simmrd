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
ar1=function(n,rho=0.5) {x=matrix(rho,n,n)^(toeplitz(1:n)-1);diag(x)=1;x}
std=function(x) (x-mean(x))/sd(x)
parthcorr=function(x,n) {
  keys=c('ar','toeplitz','cs')
  boo=sapply(keys,function(h) grepl(h,tolower(x),fixed=TRUE))
  if(boo[1]) return(ar1(n,as.numeric(substr(x,start=5,stop=nchar(x)-1))))
  if(boo[2]) return(1/toeplitz(1:n))
  if(boo[3]) {a_=as.numeric(substr(x,start=4,stop=nchar(x)-1)); return(diag(n)+a_-diag(a_,n))}
  return(diag(n))
}
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
parthstd=function(bx,by,bxse,byse,maf,nx,ny,MR_standardization_type) {
  mst=tolower(MR_standardization_type)
  if(mst=='none') return(list(bx=bx,bxse=bxse,by=by,byse=byse))
  if(mst=='z') return(list(bx=bx/bxse,bxse=bx/bx,by=by/byse,byse=by/by))
  # will then do qi & chatterjee if top two didn't return
  by=by/byse/sqrt(ny); byse=1/sqrt(ny)
  bx=bx/bxse/sqrt(nx); bxse=1/sqrt(nx)
  return(list(bx=bx,bxse=bxse,by=by,byse=byse))
}
`%!in%`=Negate(`%in%`)
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
classIVs=function(ix,uhpix,chpix) {
  keys=c('UHP','CHP')
  ll=list(uhpix,chpix)
  boo=sapply(1:2,function(h) ix %in% ll[[h]])
  boo=matrix(boo,nr=length(ix))
  cl=c();for(i in 1:nrow(boo)) {toa=keys[which(boo[i,])];cl[i]=ifelse(length(toa)==0,'valid',toa)}
  return(cl)
}
setf=function(bxunstd,nX,fix_Fstatistic_at) {
  # currently agnostic to LD structure
  m=nrow(B);p=ncol(bxunstd)
  v=diag(bxunstd%*%t(bxunstd)); 
  ord=cbind(1:m,v);
  ord=ord[order(v,decreasing=FALSE),]
  #bxunstd=bxunstd[ord[,1],]
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
  return(list(ix=ix,fs=fs))
}
plot_simdata=function() {
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
