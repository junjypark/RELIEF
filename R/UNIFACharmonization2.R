UNIFACharmonize2=function(data.mat, batches=NULL, sources=NULL, covariates=NULL,
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
  
  if (is.null(sources)){ sources=rep(1,p) }
  sources.f=as.factor(sources); sources = as.numeric(sources.f)
  source.id=unique(sources); n.sources=length(source.id)
  if (verbose) cat(paste0("[UNIFAC] ",n.sources," data source identified\n"))
  
  mn=apply(data.mat,1,mean)
  data.mat=data.mat-mn
  regfit=lm(t(data.mat)~batches.f)
  data.mat=(regfit$residuals)
  sigma.features=  tcrossprod(sigma(regfit),rep(1,nrow(data.mat)))
  data.mat=t(data.mat)/sigma.features
  
      
  # XZ=model.matrix(~covariates+batches.f)
  # X=XZ[,which(attr(XZ,"assign")<2)]
  # 
  # beta=regfit$coefficients[1:ncol(X),]
  # Xbeta=t(X%*%beta)

  sub.batch = unlist(lapply(c(1,n.batches), combn, x = batch.id, simplify = FALSE), recursive = FALSE)
  sub.source = unlist(lapply(c(1,n.sources), combn, x = source.id, simplify = FALSE), recursive = FALSE)
  sub.source=unique(sub.source)
  
  nvec=rep(NA, n.batches); mvec=rep(NA, n.sources)
  sigma.mat.batch=Matrix(1, p, n)
  
  sig=sigma.rmt(data.mat)
  data.mat2=data.mat
  for (b in 1:n.batches){
    order.temp.batch=which(batches==batch.id[b])
    nvec[b]=length(order.temp.batch)
    for (s in 1:n.sources){
      order.temp.source=which(sources==source.id[s])
      mvec[s]=length(order.temp.source)
      s=sigma.rmt(data.mat[order.temp.source, order.temp.batch])
      sigma.mat.batch[order.temp.source, order.temp.batch]=sigma.mat.batch[order.temp.source, order.temp.batch]*s
      data.mat[order.temp.source, order.temp.batch]=data.mat[order.temp.source, order.temp.batch]/s
    }
  }
  
  lambda.set=matrix(NA, length(sub.source),length(sub.batch))
  for (b in 1:length(sub.batch)){
    for (s in 1:length(sub.source)){
      lambda.set[s,b]=sqrt(sum(mvec[sub.source[[s]]]))+sqrt(sum(nvec[sub.batch[[b]]]))
    }
  }
  
  index.set.batch=foreach(b=1:length(sub.batch))%do%{ which(batches%in%sub.batch[[b]]) }
  index.set.source=foreach(s=1:length(sub.source))%do%{ which(sources%in%sub.source[[s]]) }
  
  estim=matrix(list(), length(sub.source), length(sub.batch))
  for (b in 1:length(sub.batch)){
    for (s in 1:length(sub.source)){
      estim[[s,b]]=Matrix(0, p, n, sparse=T)
    }
  }
  
  bool=TRUE
  count=1; crit0=0
  foo=matrix(1:(length(sub.source)*length(sub.batch)),length(sub.source))
  
  if (verbose) {
    cat(paste0("[UNIFAC] Start optimizing...\n"))
    pb = txtProgressBar(min = 0, max=max.iter, initial=0, char="-", style = 3)
  }
  
  Xbeta=0
  count.vec=c()
  for (bool2 in 1:5){
    #Update R,I,E
    while (bool){
      if (verbose){  setTxtProgressBar(pb, count)  }
      crit0.old = crit0
      
      nuc.temp=matrix(NA,length(sub.source), length(sub.batch))
      for (b in length(sub.batch):1){
        for (s in length(sub.source):1){
          temp=softSVD( (data.mat-Xbeta-Reduce("+", estim[-foo[s,b]]))[index.set.source[[s]],index.set.batch[[b]]],lambda.set[s,b])
          estim[[s,b]][index.set.source[[s]],index.set.batch[[b]]]=temp$out
          nuc.temp[s,b]=temp$nuc
        }
      }
      
      crit0 = 1/2*frob(data.mat-Xbeta-Reduce("+", estim))+sum(lambda.set*nuc.temp,na.rm=T)
      if (abs(crit0.old-crit0)<eps){ bool=FALSE }
      else if (count==max.iter){ bool=FALSE}
      else{ count = count+1 }
    }
    
    #Update Beta    
    E=data.mat-Reduce("+", estim)
    fit.reg=lm(as.matrix(t(E))~covariates)
    Xbeta=t(fit.reg$fitted.values)
    count.vec=c(count.vec, crit0)
  }
  
  
  
  if (verbose & count<max.iter){
    setTxtProgressBar(pb, max.iter)
    cat(paste0("\n[UNIFAC] Convergence reached. Finish harmonizing.\n"))
  }
  if (verbose & count==max.iter){
    cat(paste0("\n[UNIFAC] Convergence not reached. Increase max.iter.\n"))
  }
  
  E0=(data.mat-Reduce("+", estim))
  G0=Reduce("+",estim[foo[1:length(index.set.source), length(index.set.batch)]])
  I0=sigma.mat.batch*Reduce("+", estim[1:length(index.set.source),-length(index.set.batch)])
  
  sig=sd(data.mat2-sigma.mat.batch*(G0+I0))
  
  E=sig*sigma.features*E0
  G=sigma.mat.batch*sigma.features*G0
  I=sigma.mat.batch*sigma.features*I0
  
  harmonized=Xbeta+G+E
  
  return(list(dat.unifac=harmonized,Xbeta=Xbeta,batches=batches.f,
              sigma.features=sigma.features, sig=sig,
              sigma.batch=unique(as.numeric(sigma.mat.batch)),
              G=G, I=I, E=E))
}
