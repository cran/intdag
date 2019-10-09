#' A DAG function for observational data
#'
#' This function allows you to learn the DAG structure from observational data
#' @param X The n by p data matrix  
#' @param lambda tuning parameter for the first penalty of the adjacency matrix
#' @param tau tuning parameter of the TLP function
#' @param rho the ADMM penalty parameter, default is 1
#' @param A_NZ0 An p by p matrix indicating nonzero elements as initial values
#' @param A0 An p by p matrix as initial values for A
#' @param opts.tol Tolerance for convergence 
#' @param maxIter maximum number of iterations in ADMM loop
#' @return Estimated adjacency matrix
#' @export
#' @useDynLib intdag, .registration = TRUE
#' @examples
#' p <- 10
#' amat <- matrix(0, p, p)
#' amat[2:p, 1] <- 1
#' Sig <- seq(1, 0.5, length.out=p)
#' X <- rmvDAG_obs(100, amat, Sig)
#' out <- obsDAG(X, 5, 0.01)
obsDAG <- function(X, lambda, tau, rho = 1, A_NZ0 = NULL, A0 = NULL, 
  opts.tol = 1e-4, maxIter = 1000){
  p <- ncol(X)
  n <- nrow(X)

  # default values for A_NZ0
  if(is.null(A_NZ0)){
    A_NZ0 <- matrix(0, p, p)
    nonzero <- 0
  } else{
    if(all.equal(dim(A_NZ0), c(p, p)))
      nonzero <- sum(A_NZ0 == 1)
    else stop("Invalid input A_NZ0!")
  }

  # default values for A0
  if(is.null(A0)){
    A0 <- matrix(0, p, p)
  }else{
    if(!all.equal(dim(A0), c(p, p))) {
      stop("Invalid input A0!")
    }
  }

  # fix the variance as 1
  sigma <- rep(1, p)

  # The inverse of XTX to update A, each dimension corresponds to each row in A
  XTX <- t(X)%*%X
  XTX_inv <- array(0, c(p-1, p-1, p))
  for(i in 1:p){
    XTX_inv[,,i] <- solve((XTX[-i,-i] + rho*diag(p-1)))
  }

  out <- .C("DAG_obs", X=as.double(X), A=as.double(A0), m=as.integer(p), 
    n=as.integer(n), lambda=as.double(lambda), tau=as.double(tau), 
    A_NZ=as.integer(A_NZ0), nonzero=as.integer(nonzero), sigma=as.double(sigma), 
    tol=as.double(opts.tol), obj=as.double(0), XTX=as.double(XTX),
    XTX_inv=as.double(XTX_inv), rho=as.double(rho), maxIter=as.integer(maxIter))

  NZout   <- matrix(out[[7]], p, p)
  Aout    <- matrix(out$A, p, p)
  Aout    <- Aout * NZout
  return(list(A = Aout))
}
