n <- 10
sd <- 1e-100
pitch <- 90
roll <- 0
yaw <- 0
eulertrue <- cbind(pitch, roll,yaw)
Rtrue <- euler_to_rotation(pitch, roll,yaw, units = "degrees")
Rtrue


qtrue <- rotation_to_quaterion(Rtrue)
qtrue
quaterion_to_rotation(qtrue)
euler_to_quaterion(pitch,roll,yaw)
qtrue

rotation_to_euler(Rtrue)






#
# X <- matrix(rnorm(3 * n, sd = 1),ncol = 3)
# X <- t(apply(X,1,FUN = function(x) {x/norm(x, "2")^2}))
# X
# apply(Y,2,norm,"2")
#
# Y <- t(apply(X,1,function(x) {Rtrue %*% x})) + matrix(rnorm(3*n, mean =0, sd = sd), nrow = n)
# Y
# Y <- t(apply(Y,1,FUN = function(x) {x/norm(x, "2")}))
# Y
# apply(Y,1,norm,"2")
#
#
#
# out <- find_rotation(X,Y, output = "rotation", method = "q", sd = 1e-5)
# R <- out$rotation
# cov <- out$cov
# R
# Rtrue
# Rtrue-R
#
# rotation_to_euler(out$rotation)
# rotation_to_euler(Rtrue)
#
#
#
# cov
#
