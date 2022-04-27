# functions to estimate calibration matrix


#' Calibration model for
#'
#' @param V Matrix (n x 3) of observed response where each row denotes v_i = (v_{ix}, v_{iy}, v_{iz})
#' @param G Matrix (n x 3) of fixed gravity vectors where each row denotes g_i = (g_{ix}, g_{iy}, g_{iz}). Each row should have unit norm.
#'
#' @return List of estimated coefficients from linear model: BetaHat (3x3) matrix of first-order effects, Beta0Hat (1x3) vector of intercepts,
#' Ghat (n x3) matrix of estimated gravity vectors, sd is estimated standard error of responses.
#' @export
#'
#' @examples
calibration <- function(V, G) {

  BetaHat <- matrix(NA, nrow= 3, ncol =3)
  Beta0Hat <- c()
  sumres <- c()
  n <- nrow(V)
  for (i in 1:3) {
    mod <- lm(V[,i] ~ 1 + G)
    coef <- coef(mod)
    Beta0Hat[i] <- coef[1]
    BetaHat[,i] <- coef[2:4]
    sumres[i] <- sum(mod$residuals^2)
  }


  Ghat <- (V - Beta0Hat) %*% solve(BetaHat)
  sd <- sqrt(sum(sumres)/(n*3-2))

  return(list(coef = list(BetaHat = BetaHat, Beta = Beta), Ghat = Ghat, V = V, G=G ,sd = sd))

}

