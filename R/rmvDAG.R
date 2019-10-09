#' A random graph data generation function 
#'
#' This function generates random observations from a DAG graph
#' @param n The sample size  
#' @param amat The adjacency matrix of a DAG
#' @param Sig The error variance of each node  
#' @return Gaussian data with the given sample size
#' @importFrom stats rnorm
#' @export
#' @examples
#' amat=matrix(c(0,1,0,0),2,2)
#' rmvDAG_obs(50,amat)

rmvDAG_obs <- function(n, amat, Sig = NULL){
    stopifnot(is.numeric(n), is.matrix(amat)) 
    p = dim(amat)[1]
    if(is.null(Sig)) Sig <- rep(1,p)
    errMat <- matrix(rnorm(n * p), nrow = n)
    X <- matrix(0, n, p)
    X[, 1] <- errMat[, 1] * sqrt(Sig[1])
    for (j in 2:p) {
        ij <- 1:(j - 1)
        X[, j] <- X[, ij, drop = FALSE] %*% amat[j,ij] + errMat[, j] * sqrt(Sig[j])
    }
    X
}
