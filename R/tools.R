library(denoiseR) 
frob=function(X){ sum(X^2,na.rm=T) }

sigma.rmt=function(X){ estim_sigma(X,method="MAD") }

softSVD=function(X, lambda){
  svdX=svd(X)
  nuc=pmax(svdX$d-lambda,0)
  out=tcrossprod(svdX$u, tcrossprod( svdX$v,diag(nuc) ))
  return(list(out=out, nuc=sum(nuc)))
}
