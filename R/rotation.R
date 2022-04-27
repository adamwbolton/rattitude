#' Find optimal rotation
#' @param X matrix of vectors to be aligned with Y
#' @param Y matrix of vectors
#' @param method  method for finding optimal rotation The default is Davenport's q-method.
#' @param sd standard deviation
#' @param output format for rotation. Output can be set to: (1) "euler" (pitch, roll, yaw), (2) "rotation" (3x3 rotation matrix) or (3) quaternion (rotation quaternion). Defaults to euler angles.
#' @export
#' @return
#' @return Rotation in the form request by output. Returns covariance estimate is \code{sd} is provided.
#'
#' @examples

find_rotation <- function(X, Y, method = "q", output = "euler", sd = NA) {

  if (method == "q") {
    # attiude profile matrix (APM)
    n <- nrow(X)
    # B <- (t(X)  %*% diag(1/n,nrow = n,ncol = n) %*%  Y)

    B <- matrix(rowMeans(sapply(1:n, FUN = function(i) Y[i,] %*% t(X[i,]))),3,3, byrow = T)

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
    q <- eigenK$vectors[,1]
  }


  if (output == "euler") {
    out <- quaterion_to_euler(q)
  }

  if (output == "rotation") {
    out <- quaterion_to_rotation(q)
  }

  if (output == "quaterion") {
    out <- q
  }

  if (!is.na(sd)) {
    cov <- get_covariance(X,Y, units = output, sd= sd)
  }

  if(is.na(sd)) {
    simpleWarning("Must provide standard deviation for covariance estimate.")
  }

  return(list(rotation = out, cov = cov))
}




#' Get rotation matrix based on provided pitch, roll and yaw angles.
#'
#' @param pitch,roll,yaw euler angles of rotation.
#' @param order order of rotation. only currently takes "xyz"
#' @param units units for angles: degrees or radians.
#'
#' @return Rotation matrix.
#' @export
#'
#' @examples
euler_to_rotation <- function(pitch, roll, yaw, order = "xyz", units = "degrees") {

  if (units == "degrees") {

    cosp <- cos(pitch*pi/180)
    sinp <- sin(pitch*pi/180)

    cosr <- cos(roll*pi/180)
    sinr <- sin(roll*pi/180)

    cosy <- cos(yaw*pi/180)
    siny <- sin(yaw*pi/180)

  }

  if (units == "radians") {

    cosp <- cos(pitch)
    sinp <- sin(pitch)

    cosr <- cos(roll)
    sinr <- sin(roll)

    cosy <- cos(yaw)
    siny <- sin(yaw)

  }

  if(units  %in% c("radians", "degrees") != T) {
    stop("undefined units.")
  }
  Rp <- matrix(c(cosr,0, sinr,
                 0,1,0,
                 -sinr,0,cosr), nrow = 3, ncol = 3,byrow = T)

  Rr <- matrix(c(1,0,0,
                 0, cosp, - sinp,
                 0,sinp,cosp), nrow = 3,ncol = 3,byrow =T)

  Ry <- matrix(c(cosy, -siny,0,
                 siny,cosy,0,
                 0,0,1), nrow = 3,ncol = 3,byrow = T)

  R <- Rp %*% Rr %*% Ry
  return(R)

}



#' Quaterion to euler angles
#'
#' @param q quaternion of rotation.
#'
#' @return sequence of euler angles specifying rotation.
#' @export
#'
#' @examples
quaterion_to_euler <- function(q) {
  q <- Re(q)
  # roll <- atan2(2*(q[1]*q[2] + q[3]*q[4]), 1- 2*(q[2]^2 + q[3]^2))
  yaw <- (atan(2*(q[1]*q[2] + q[3]*q[4])/ (1- 2*(q[2]^2 + q[3]^2)))*(180/pi) + 360) %% 180
  roll <- (asin(2*(q[1]*q[3]-q[4]*q[2]))*(180/pi) + 360) %% 180
  # yaw <- atan2(2*(q[1]*q[4] + q[2]*q[3]), 1- 2*(q[3]^2 + q[4]^2))
  pitch <- (atan(2*(q[1]*q[4] + q[2]*q[3])/(1- 2*(q[3]^2 + q[4]^2)))*(180/pi) + 360) %% 180

  euler <- list(pitch = pitch, roll = roll, yaw = yaw)

  return(euler)
}


#' Rotation matrix to quaternion of rotation.
#'
#' @param R rotation matrix
#'
#' @return quaternion of rotation.
#' @export
#'
#' @examples
rotation_to_quaterion <- function(R) {

    eigenR <- eigen(R)
    idx <- which(eigenR$values == 1)
    qtilde <- eigenR$vectors[,idx]
    theta <- acos((sum(diag(R))-1)/2)
    q <- c(qtilde*sin(theta/2),cos(theta/2))
    q <- Re(q)
    return(q)
}

quaterion_to_rotation <- function(q) {

  # qtilde <- q[1:3]
  # qr <- q[4]
  # Q <- matrix(c(0,-qtilde[3], qtilde[2],
  #               qtilde[3],0,-qtilde[1],
  #               -qtilde[2],qtilde[1],0),
  #             nrow = 3, ncol = 3, byrow = T)

  # R <- (qr^2 - norm(qtilde,"2"))*diag(3) + 2*(qtilde %*% t(qtilde)) + 2*qr*Q

  q1 <- q[1]
  q2 <- q[2]
  q3 <- q[3]
  q4 <- q[4]

  # First row of the rotation matrix
  r11 = 1 - 2*q2^2 - 2*q3^2
  r12 = 2*q1*q2 - 2*q3*q4
  r13 = 2*q1*q3 + 2*q2*q4

  # Second row of the rotation matrix
  r21 = 2*q1*q2 + 2*q3*q4
  r22 = 1 - 2*q1^2 - 2*q3^2
  r23 = 2*q2*q3 - 2*q1*q4

  # Third row of the rotation matrix
  r31 = 2*q1*q3 - 2*q2*q4
  r32 = 2*q2*q3 + 2*q1*q4
  r33 = 1 - 2*q1^2 - 2*q2^2

  # 3x3 rotation matrix
  R <- Re(matrix(c(r11,r12,r13,r21,r22,r23,r31,r32,r33), byrow = T, nrow = 3, ncol = 3))



  return(R)

}

#' Converts rotation matrix to euler angles.
#'
#' @param R rotation matrix
#'
#' @return sequence of euler angles specifying rotation.
#' @export
#'
#' @examples
rotation_to_euler <- function(R) {
  pitch<- (acos(R[3,3])*(180/pi) + 360) %% 180
  roll <- (atan2(R[3,1],R[3,2])*(180/pi) + 360) %% 180
  yaw <- (-atan2(R[1,3],R[2,3])*(180/pi) + 360) %% 180
  return(list(pitch=pitch, roll=roll, yaw=yaw))

}


#' Euler angles to quaternion of rotation
#'
#' @param pitch,roll,yaw euler angles
#'
#' @return quaternion of rotation
#' @export
#'
#' @examples
euler_to_quaterion <- function(pitch, roll, yaw) {

  q <- rotation_to_quaterion(euler_to_rotation(pitch,roll,yaw))

  return(q)
}


#' QUEST method for finding optimal rotation.
#'
#' @param X
#' @param Y
#'
#' @return
#' @export
#'
#' @examples
quest <- function(X,Y) {

  n <- nrow(X)
  B <- (t(X)  %*% Y)/n

  S <- B + t(B)
  sigma <- sum(diag(S))/2
  z  <- c(B[2,3] - B[3,2],B[3,1] - B[1,3], B[1,2] - B[2,1])

  # get adj. of S
  Sadj <- det(S)*solve(S)

  kappa <- sum(diag(Sadj))
  Delta <- det(S)
  # define quantities for characteristic poly.
  a <- as.numeric(sigma^2 + kappa)
  b <- as.numeric(sigma^2 + t(z) %*% z)
  c <- as.numeric(Delta + t(z) %*% S %*% z)
  d <- as.numeric(t(z) %*% S^2 %*% z)

  roots <- polyroot(c((a*b+c*sigma-d), -c , -(a+b),0,1))
  lambda <- max(Re(roots))
  # polyroot(c(1,0,-(a+b),-c,(a*b+c*sigma-d)))


  # opt_fn <- function(x) {
  #   x^4 - (a+b)*x^2 - c*x + (a*b+c*sigma-d)
  # }
  # res <- optim(par = 1,fn = opt_fn, method = "Brent", lower = 0, upper = 1)
  # opt <- res$par

  Gibbs <- solve((lambda + sigma)*diag(3) - S) %*% z
  q <- 1/(sqrt(1 + norm(Y,"2")))*c(Gibbs,1)
  return(q)
}




