#' Find optimal rotation matrix R aligning set of vectors X and Y
#'
#'
#' @param X n x 3 matrix of vectors to be aligned with Y
#' @param Y n x 3 matrix of vectors
#' @param method the method denotes how rotation matrix R is found. The default is Davenport's q-method.
#'
#' @return
#' @export
#'
#' @examples
find_rotation <- function(X, Y, method = "q", output = "euler") {
  # attiude profile matrix (APM)
  n <- nrow(X)
  B <- (t(X)  %*% Y)/n

  # B <- matrix(rowMeans(sapply(1:n, FUN = function(i) Y[i,] %*% t(X[i,]))),3,3, byrow = T)

  S <- B + t(B)
  mu <- sum(diag((B)))
  z <- c(B[2,3] - B[3,2],B[3,1] - B[1,3], B[1,2] - B[2,1])
  # davenport matrix K
  K <- matrix(NA, nrow = 4, ncol = 4)
  K[1:3,1:3] <- S - mu*diag(3)
  K[4,4] <- mu
  K[4,1:3] <- z
  K[1:3,4] <- z

  q <- eigen(K)$vectors[,1]
  if (output == "euler") {
    out <- quaterion_to_euler(q)
  }

  if (output == "rotation") {
    out <- quaterion_to_rotation(q)
  }

  if (output == "quaterion") {
    out <- q
  }

  return(out)
}




#' Get rotation matrix based on provided pitch, roll and yaw angles.
#'
#' @param pitch
#' @param roll
#' @param yaw
#' @param order Order of rotation.
#' @param units Degrees or radians.
#'
#' @return Rotation matrix.
#' @export
#'
#' @examples
euler_to_rotation <- function(pitch, roll, yaw, order = "zyx", units = "degrees") {

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


  Rx <- matrix(c(1,0,0,
                 0, cosp, - sinp,
                 0,sinp,cosp), nrow = 3,ncol = 3)
  Ry <- matrix(c(cosr,0, sinr,0,1,0,-sinr,0,cosr), nrow = 3, ncol = 3)
  Rz <- matrix(c(cosy, -siny,0,siny,cosy,0,0,0,1), nrow = 3,ncol = 3)
  R <- Rx %*% Ry %*% Rz
  return(R)

}



quaterion_to_euler <- function(q) {

  # roll <- atan2(2*(q[1]*q[2] + q[3]*q[4]), 1- 2*(q[2]^2 + q[3]^2))
  yaw <- atan(2*(q[1]*q[2] + q[3]*q[4])/ (1- 2*(q[2]^2 + q[3]^2)))
  pitch <- -asin(2*(q[1]*q[3]-q[4]*q[2]))
  # yaw <- atan2(2*(q[1]*q[4] + q[2]*q[3]), 1- 2*(q[3]^2 + q[4]^2))
  roll <- -atan(2*(q[1]*q[4] + q[2]*q[3])/(1- 2*(q[3]^2 + q[4]^2)))

  euler <- list(pitch = pitch*(180/pi) %% 180, roll = roll*(180/pi) %% 180, yaw = yaw*(180/pi) %% 180)

  return(euler)
}


quaterion_to_rotation <- function(q) {

  qtilde <- q[1:3]
  qr <- q[4]
  Q <- matrix(c(0,-qtilde[3], qtilde[2],
                qtilde[3],0,-qtilde[1],
                -qtilde[2],qtilde[1],0),
              nrow = 3, ncol = 3)

  R <- (qr^2 - norm(qtilde,"2")^2)*diag(3) + 2*(qtilde %*% t(qtilde)) + 2*qr*Q
  return(R)

}

euler_to_quaterion <- function(pitch, roll, yaw) {

  cosp <- cos(pitch*0.5*pi/180)
  sinp <- sin(pitch*0.5*pi/180)

  cosr <- cos(roll*0.5*pi/180)
  sinr <- sin(roll*0.5*pi/180)

  cosy <- cos(yaw*0.5*pi/180)
  siny <- sin(yaw*0.5*pi/180)
  q <- c()

  q[1] = sinr * cosp * cosy - cosr * sinp * siny
  q[2] = cosr * sinp * cosy + sinr * cosp * siny
  q[3] = cosr * cosp * siny - sinr * sinp * cosy
  q[4] = cosr * cosp * cosy + sinr * sinp * siny

  return(q)
}


