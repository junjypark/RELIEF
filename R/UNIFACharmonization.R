UNIFACharmonize=function(data.mat, batches=NULL, sources=NULL, covariates=NULL,
                          scale.features=T, eps=1e-3, max.iter=1000, verbose=T){
  if (verbose) {
    if (!is.null(covariates)){ cat(paste0("[UNIFAC] Performing UNIFAC harmonization with ", ncol(covariates), " covariates\n")) }
    else{ cat(paste0("[UNIFAC] Performing UNIFAC harmonization without covariates\n")) }
  }
  if (is.null(batches)){ stop("batch information must be provided\n") }
  p=nrow(data.mat); n=ncol(data.mat)

  batches.f=as.factor(batches); batches = as.numeric(batches.f)
  batch.id=unique(batches); n.batches=length(batch.id)
  if (verbose) cat(paste0("[UNIFAC] ",n.batches," batches identified\n"))

  sources=rep(1,p)
  
  if (!is.null(covariates)){
    XZ=model.matrix(~covariates+batches.f)
    X=XZ[,which(attr(XZ,"assign")<2)]
    
    regfit=lm(t(data.mat)~XZ-1)
    
    beta=regfit$coefficients[1:ncol(X),]
    Xbeta=t(X%*%beta)
    data.mat=t(regfit$residuals)
    # sigma.features=   tcrossprod(sigma(regfit),rep(1,n))
    # data.mat=t(data.mat)/sigma.features
  } else{
    regfit=lm(t(data.mat)~batches.f)
    data.mat=t(regfit$residuals)
    # sigma.features= tcrossprod(sigma(regfit),rep(1,n))
    # data.mat=t(data.mat)/sigma.features
    Xbeta=matrix(0,p,n)
  }

  # mn=matrix(0,p,n)
  # for (i in 1:n){
  #   mn[,i]=mean(data.mat[,i])
  #   data.mat[,i]=data.mat[,i]-mn[,i]
  # }
  # 
  # for (b in 1:n.batches){
  #   dd=mean(data.mat[,batches==b])
  #   mn[,batches==b]=mn[,batches==b]+dd
  #   data.mat[,batches==b]=data.mat[,batches==b]+dd
  # }
  
  sub.batch = unlist(lapply(c(1,n.batches), combn, x = batch.id, simplify = FALSE), recursive = FALSE)

  nvec=rep(NA, n.batches); 
  sigma.mat.batch=Matrix(1, p, n)

  sig=sigma.rmt(data.mat)
  data.mat2=data.mat
  for (b in 1:n.batches){
    order.temp.batch=which(batches==batch.id[b])
    nvec[b]=length(order.temp.batch)
  
    s=sigma.rmt(data.mat[,order.temp.batch])
    sigma.mat.batch[, order.temp.batch]=sigma.mat.batch[, order.temp.batch]*s
    data.mat[, order.temp.batch]=data.mat[, order.temp.batch]/s
  }

  lambda.set=matrix(NA,1,length(sub.batch))
  for (b in 1:length(sub.batch)){
    lambda.set[b]=sqrt(p)+sqrt(sum(nvec[sub.batch[[b]]]))
  }

  index.set.batch=foreach(b=1:length(sub.batch))%do%{ which(batches%in%sub.batch[[b]]) }

  estim=matrix(list(), 1, length(sub.batch))
  for (b in 1:length(sub.batch)){
      estim[[1,b]]=Matrix(0, p, n, sparse=T)
  }

  bool=TRUE
  count=1; crit0=0
  foo=matrix(1:length(sub.batch),1)

  if (verbose) {
    cat(paste0("[UNIFAC] Start optimizing...\n"))
    pb = txtProgressBar(min = 0, max=max.iter, initial=0, char="-", style = 3)
  }

  while (bool){
    if (verbose){  setTxtProgressBar(pb, count)  }
    crit0.old = crit0

    nuc.temp=matrix(NA,1, length(sub.batch))
    for (b in length(sub.batch):1){
        temp=softSVD( (data.mat-Reduce("+", estim[-foo[1,b]]))[,index.set.batch[[b]]],lambda.set[1,b])
        estim[[1,b]][,index.set.batch[[b]]]=temp$out
        nuc.temp[1,b]=temp$nuc
    }

    crit0 = 1/2*frob(data.mat-Reduce("+", estim))+sum(lambda.set*nuc.temp,na.rm=T)
    if (abs(crit0.old-crit0)<eps){ bool=FALSE }
    else if (count==max.iter){ bool=FALSE}
    else{ count = count+1 }
  }


  if (verbose & count<max.iter){
    setTxtProgressBar(pb, max.iter)
    cat(paste0("\n[UNIFAC] Convergence reached. Finish harmonizing.\n"))
  }
  if (verbose & count==max.iter){
    cat(paste0("\n[UNIFAC] Convergence not reached. Increase max.iter.\n"))
  }
  
  E0=(data.mat-Reduce("+", estim))
  G0=Reduce("+",estim[foo[1, length(index.set.batch)]])
  I0=Reduce("+", estim[1,-length(index.set.batch)])
  
  e=as.matrix(data.mat2-sigma.mat.batch*(G0+I0))
  s2=numeric(n.batches)
  for (b in 1:n.batches){
    s2[b]=var(as.numeric(e[,batches==b]))
  }

  # s2=unique(as.numeric(sigma.mat.batch))^2
  sb=sqrt(sum(nvec*s2)/sum(nvec))
  
  # sig=sd(data.mat2-sigma.mat.batch*(G0+I0))
  # E=sb*sigma.features*E0
  # G=sb*sigma.features*G0
  # I=sb*sigma.features*I0

  E=sb*E0
  G=sigma.mat.batch*G0
  I=sb*I0
  
  harmonized=Xbeta+G+E

  return(list(dat.unifac=harmonized,Xbeta=Xbeta,batches=batches.f,
              # mn=mn,
              # sigma.features=sigma.features, sig=sig,
              sigma.batch=unique(as.numeric(sigma.mat.batch)),
              G=G, E=E,I=I))
}
