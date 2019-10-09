#' A random graph data generation function 
#'
#' This function generates random observations from a DAG graph
#' @param n The sample size  
#' @param amat The adjacency matrix of a DAG
#' @param bmat The intervention matrix
#' @param Sig The error variance of each node  
#' @return Gaussian data with the given sample size
#' @importFrom stats rnorm
#' @export
#' @examples
#' amat <- matrix(c(0,1,0,0),2,2)
#' bmat <- matrix(c(1,1,0,1),2,2)
#' rmvDAG_int(50, amat, bmat)

rmvDAG_int <-function(n, amat, bmat, Sig = NULL){
  stopifnot(is.numeric(n), is.matrix(amat), is.matrix(bmat)) 
  p = dim(amat)[1]
  w = dim(bmat)[2]
  if(is.null(Sig)) Sig <- rep(1, p) 
  sd.Sig = sqrt(Sig)
  X        = matrix(rnorm(n * w), n, w)
  Y        = matrix(0, n, p)
  Y[, 1]   = rnorm(n, mean=0, sd=sd.Sig[1]) + bmat[1, ] %*% t(X)
  for (j in 2:p) {
    ij     = 1:(j - 1)
    Y[, j] = Y[, ij, drop = FALSE] %*% amat[j, ij] + as.numeric((bmat[j, ])%*%t(X))+ rnorm(n, mean=0, sd.Sig[j])
  }
  return(cbind(Y, X)) 
}