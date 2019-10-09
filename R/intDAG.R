#' A DAG function for interventional data
#'
#' This function allows you to learn the DAG structure from interventional data
#' @param X The n by M data matrix
#' @param p The dimension of the adjacency matrix 
#' @param lambda tuning parameter for the first penalty of the adjacency matrix
#' @param lambda2 tuning parameter for the sparsity penalty of the intervention matrix
#' @param tau tuning parameter of the TLP function, default is 0.05
#' @param rho the ADMM penalty parameter, default is 1
#' @param A_NZ0 An p by M matrix indicating nonzero elements as initial values
#' @param A0 An p by M matrix as initial values for (A, B)
#' @param sigma the parameter in the variance constraint, not needed when variance.constraint is set to FALSE
#' @param Sig vector of length p, the error variances of each node, not needed when variance.constraint is set to FALSE
#' @param variance.constraint a flag indicating if the variance constraint is included, default is TRUE
#' @param opts.tol Tolerance for convergence 
#' @param maxIter maximum number of iterations in ADMM loop
#' @return A list with components
#'   \item{A}{Estimated adjacency matrix}
#'   \item{B}{Estimated intervention matrix}
#'   \item{Sig}{Estimated vector of error variances of each node}
#'   \item{sigma}{Estimated paramter in the variance constraint}
#' @export
#' @useDynLib intdag, .registration = TRUE
#' @examples
#' \donttest{p <- w <- 10
#' s0 <- p # number of edges
#' lower <- rep(0, (p*(p-1)/2)) # num of possible edges
#' nz_set <- sample(1:(p*(p-1)/2), s0) # sample a non-zero edge set
#' lower[nz_set] <- 0.5
#' amat <- matrix(0, p, p)
#' amat[lower.tri(amat)] <- lower
#' bmat <- diag(sqrt(seq(1, 1.5, length=p)))
#' Sig <- seq(1.5, 1, length = p)
#' X <- rmvDAG_int(100, amat, bmat, Sig)
#' Sig0 <- rep(1, p)
#' sigma0 <- 3
#' out <- intDAG(X, p, 2, 2, 0.05, rho=10, sigma=sigma0, Sig=Sig0)}
intDAG <- function(X, p, lambda, lambda2, tau = 0.05, rho = 1, A_NZ0 = NULL, 
  A0 = NULL, sigma = NULL, Sig = NULL, variance.constraint = TRUE, 
  opts.tol = 1e-3, maxIter = 1000){

  M <- ncol(X)
  n <- nrow(X)

  # default values for nonzero0
  if(is.null(A_NZ0)){
    A_NZ0 <- matrix(0, p, M)
    nonzero <- 0
  } else{
    if(all.equal(dim(A_NZ0), c(p, M)))
      nonzero <- sum(A_NZ0 == 1)
    else stop("Invalid input A_NZ0!")
  }

  # default values for A0
  if(is.null(A0)){
    A0 <- matrix(0, p, M)
  }else{
    if(!all.equal(dim(A0),c(p, M))) {
      stop("Invalid input A0!")
    }
  }

  B0  <- A0
  Sig <- rep(1, p)

  ## cache the inverse of XTX + rho I
  XTX <- t(X)%*%X
  XTX_inv <- array(0,c(M-1, M-1, p))
  for(i in 1:p){
    XTX_inv[,,i] <- solve((XTX[-i,-i] + Sig[i]*rho*diag(M-1)))
  }

  # col sums of X^2
  X2 = rep(0, p)
  for(j in 1:p){
    X2[j] = sum(X[, j]^2)
  }

  out <- .C("DAG_int", as.double(X), A=as.double(A0), B=as.double(B0), 
    as.integer(p), as.integer(M), as.integer(n), as.double(lambda), 
    as.double(lambda2), as.double(tau), A_NZ=as.integer(A_NZ0), 
    as.integer(nonzero), sigma=as.double(1), sigmaX=as.double(1), 
    Sig=as.double(Sig), as.double(opts.tol), obj=as.double(0), as.double(XTX), 
    as.double(XTX_inv), as.double(X2), as.double(rho), as.double(rho), 
    as.integer(maxIter))
  
  NZout   <- matrix(out$A_NZ, p, M)
  nonzero <- sum(NZout == 1)

  Aout  <- matrix(out$A, p, M)
  Bout  <- matrix(out$B, p, M)
  NZout <- matrix(out$A_NZ, p, M)
  sigma <- out$sigma
  Sig   <- out$Sig
  # transform into original parameteriation
  Aout  <- diag(1/Sig) %*% Aout
  Bout  <- diag(1/Sig) %*% Bout
  Sig   <- 1/Sig^2

  if(variance.constraint){
    C0  <- matrix(0, p, M-p)
    U0  <- matrix(0, p, M)
    
    out <- .C("DAG_int_var", as.double(X), A=as.double(A0), B=as.double(B0), 
      C=as.double(C0), U=as.double(U0), as.integer(p), as.integer(M), 
      as.integer(n), as.double(lambda), as.double(lambda2), as.double(tau), 
      A_NZ=as.integer(A_NZ0), as.integer(nonzero), sigma=as.double(sigma), 
      sigmaX=as.double(1), Sig=as.double(1/Sig), as.double(opts.tol), 
      obj=as.double(0), as.double(XTX), as.double(XTX_inv), as.double(rho), 
      as.double(rho), as.double(rho), as.integer(maxIter))
    
    NZout   <- matrix(out$A_NZ, p, M)
    nonzero <- sum(NZout == 1)
  
    Aout  <- matrix(out$A, p, M)
    Bout  <- matrix(out$B, p, M)
    sigma <- out$sigma
    Sig   <- out$Sig
    # transform into original parameteriation
    Aout  <- diag(1/Sig) %*% Aout
    Bout  <- diag(1/Sig) %*% Bout
    Sig   <- 1/Sig^2
  }
  return(list(A = Aout, B = Bout, Sig = Sig, sigma = sigma))
}

