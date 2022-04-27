#' Find Covariance Matrix for Rotation
#'
#' @param X Matrix of vector observations
#' @param Y Matrix of vector observations
#'
#' @return Covariance matrix
#' @export
#'
#' @examples
get_covariance <- function(X,Y,units = "euler", sd = NA) {

  if (is.na(sd)) {
    stop("Must provide standard error for vector observations.")
  }

  n <- nrow(X)
  # attiude profile matrix
  B <- (t(X)  %*% Y)
  S <- B + t(B)
  sigma <- sum(diag((B)))
  z <- c(B[2,3] - B[3,2],B[3,1] - B[1,3], B[1,2] - B[2,1])
  # davenport matrix K
  K <- matrix(NA, nrow = 4, ncol = 4)
  K[1:3,1:3] <- S - sigma*diag(3)
  K[4,4] <- sigma
  K[4,1:3] <- z
  K[1:3,4] <- z
  eigenK <- eigen(K)
  # eigen vector associated with largest eigenvalue is the optimal rotation.
  q <- eigenK$vectors[,1]

  # get predicted Yhat vectors to find sigma2 assuuming Y \sim N(0,\sigma^2 I_3)
  R <- quaterion_to_rotation(q)
  Yhat <- t(apply(X,1,function(x) R %*% x))
  # sig2 <- sum((Y - Yhat)^2)/(3*n-2)

  CovCart <- sd^2*solve(matrix(rowMeans(apply(X,1, function(x) diag(3) -( (R %*% x) %*% t((R %*% x))))), nrow = 3, ncol = 3,byrow= T))

  if (units == "rotation") {
    out <- CovCart
    rownames(out) <- c("x", "y", "z")
    colnames(out) <- c("x", "y", "z")
  }

  if (units == "quaterion") {
    out <- CovCart/4
    rownames(out) <- c("i", "j", "k")
    colnames(out) <- c("i", "j", "k")
  }

  if (units == "euler") {
    H <- cartesian_to_euler_covariance(q)
    CovEuler <- H %*% CovCart %*% t(H)
    out <- CovEuler
    rownames(out) <- c("pitch", "roll", "yaw")
    colnames(out) <- c("pitch", "roll", "yaw")
  }

  return(out)
}



#' Convert Cartesian Covariance Matrix to Euler Angle Covariance Matrix
#'
#' @param q Unit quaterion.
#'
#' @return H matrix (3x3) to properly transform cartesian covariance matrix to euler covariance matrix.
#' @export
#'
#' @examples
cartesian_to_euler_covariance <- function(q) {

  euler <- quaterion_to_euler(q)
  p <- euler$pitch
  r <- euler$roll
  y <- euler$yaw

  R11 <- function(p,r,y) {cos(y)*cos(p)}
  R12 <- function(p,r,y) {cos(y)*sin(p)*sin(y)-sin(y)*cos(r)}
  R13 <- function(p,r,y) {cos(y)*sin(p)*cos(r)+sin(y)*sin(r)}

  R21 <- function(p,r,y) {sin(y)*cos(p)}
  R22 <- function(p,r,y) {sin(y)*sin(p)*sin(r)+cos(y)*cos(r)}
  R23 <- function(p,r,y) {sin(y)*sin(p)*cos(r)-cos(y)*sin(r)}

  R31 <- function(p,r,y) {-sin(p)}
  R32 <- function(p,r,y) {cos(p)*sin(r)}
  R33 <- function(p,r,y) {cos(p)*cos(r)}


  R1 <- function(p,r,y) {c(R11(p,r,y), R21(p,r,y), R31(p,r,y))}
  R2 <- function(p,r,y) {c(R12(p,r,y), R22(p,r,y), R32(p,r,y))}
  R3 <- function(p,r,y) {c(R13(p,r,y), R23(p,r,y), R33(p,r,y))}

  R1p <- Deriv(R1,"p")
  R1r <- Deriv(R1,"r")
  R1y <- Deriv(R1,"y")

  R2p <- Deriv(R2,"p")
  R2r <- Deriv(R2,"r")
  R2y <- Deriv(R2,"y")

  R3p <- Deriv(R3,"p")
  R3r <- Deriv(R3,"r")
  R3y <- Deriv(R3,"y")

  gradRp <- cbind(R1p(p,r,y),R2p(p,r,y),R3p(p,r,y))
  gradRr <- cbind(R1r(p,r,y),R2r(p,r,y),R3r(p,r,y))
  gradRy <- cbind(R1y(p,r,y),R2y(p,r,y),R3y(p,r,y))

  H1 <- colSums(sapply(1:3,FUN = function(i) {xprod(gradRp[,i],R[,i])}))/2
  H2 <- colSums(sapply(1:3,FUN = function(i) {xprod(gradRr[,i],R[,i])}))/2
  H3 <- colSums(sapply(1:3,FUN = function(i) {xprod(gradRy[,i],R[,i])}))/2

  Hinv <- cbind(H1,H2,H3)
  H <- solve(Hinv)

  return(H)
}

#' Cross Product
#'
#' @param x Vector (3x1)
#' @param y Vector (3x1)
#'
#' @return Cross product of two vectors.
#' @export
#'
#' @examples
xprod <- function(x,y, i = 1:3) {
  # Project inputs into 3D, since the cross product only makes sense in 3D.
  To3D <- function(x) head(c(x, rep(0, 3)), 3)
  x <- To3D(x)
  y <- To3D(y)

  # Indices should be treated cyclically (i.e., index 4 is "really" index 1, and
  # so on).  Index3D() lets us do that using R's convention of 1-based (rather
  # than 0-based) arrays.
  Index3D <- function(i) (i - 1) %% 3 + 1

  # The i'th component of the cross product is:
  # (x[i + 1] * y[i + 2]) - (x[i + 2] * y[i + 1])
  # as long as we treat the indices cyclically.
  return (x[Index3D(i + 1)] * y[Index3D(i + 2)] -
            x[Index3D(i + 2)] * y[Index3D(i + 1)])
}

