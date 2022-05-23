rearrange.unifac=function(fit){
  batches=fit$batches
  batches.id=unique(batches)
  n.batches=length(batches.id)
  p=nrow(fit$dat.unifac); n=ncol(fit$dat.unifac)
  Xbeta.ordered=G.ordered=I.ordered=E.ordered=matrix(NA, p, n)
  orderlist=NULL
  tempv=hclust(dist(fit$I))$order
  for (b in 1:n.batches){
    batch=batches.id[b]
    index=which(batches==batch)
    temp=hclust(dist(t(fit$I[,index])))$order
    index2=index[temp]
    Xbeta.ordered[,index]=as.matrix(fit$Xbeta)[tempv,index2]
    orderlist=c(orderlist, index2)
    G.ordered[,index]=as.matrix(fit$G)[tempv,index2]
    I.ordered[,index]=as.matrix(fit$I)[tempv,index2]
    E.ordered[,index]=as.matrix(fit$E)[tempv,index2]
  }

  for (j in 1:p){
    Xbeta.ordered[j,]=Xbeta.ordered[j,]-mean(Xbeta.ordered[j,])
  }

  return(list(
    Xbeta=Matrix(Xbeta.ordered),
    G=Matrix(G.ordered),
    I=Matrix(I.ordered),
    E=Matrix(E.ordered),
    orderlist=orderlist,
    tempv=tempv
  ))
}
