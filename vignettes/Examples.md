---
title: "Examples"
output: 
  rmarkdown::github_document
vignette: >
  %\VignetteIndexEntry{Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(rattitude)
```

# Rotations

Rotations can be specified either as rotation matrices, sequences of Euler angles, or rotation quaternions. One can convert between the three with the following functions: `euler_to_rotation`, `euler_to_quaternion`, `rotation_to_quaternion`, `rotation_to_euler`, `quaternion_to_euler`, `quaternion_to_rotation`. Let's start with a rotation of 90 degrees about the x-axis and convert this to a rotation matrix:


```r
# toy data
pitch <- 90
roll <- 0
yaw <- 0
R <- euler_to_rotation(pitch, roll,yaw, units = "degrees")
# rotation matrix
R
#>      [,1]         [,2]          [,3]
#> [1,]    1 0.000000e+00  0.000000e+00
#> [2,]    0 6.123234e-17 -1.000000e+00
#> [3,]    0 1.000000e+00  6.123234e-17
```

converting this back to Euler angles


```r
# convert back to euler angles
rotation_to_euler(R)
#> $pitch
#> [1] 90
#> 
#> $roll
#> [1] 0
#> 
#> $yaw
#> [1] 0
```

and then to a rotation quaternion

```r
# rotation to quaternion
q <- rotation_to_quaterion(R)
q
#> [1] 0.7071068 0.0000000 0.0000000 0.7071068
```

and back to Euler angles

```r
quaterion_to_euler(q)
#> $pitch
#> [1] 90
#> 
#> $roll
#> [1] 0
#> 
#> $yaw
#> [1] 0
```



# Calibration

The calibration for the inertial-based three-axis accelerometers is done through the following first-order model:

$$
v_x = \beta_{0_x} + \beta_{x_x} g_x + \beta_{y_x} g_y + \beta_{z_x} g_z + \epsilon_x
$$


$$
v_y = \beta_{0_y} + \beta_{y_x} g_x + \beta_{y_y} g_y + \beta_{y_y} g_z + \epsilon_y
$$


$$
v_z = \beta_{0_z} + \beta_{z_x} g_x + \beta_{z_y} g_y + \beta_{z_z} g_z + \epsilon_z
$$

where $\epsilon = (\epsilon_x,\epsilon_y, \epsilon_z)^T \sim N(0,\sigma^2 I_3)$. For the x-axis accelerometer, $\beta_{0_x}$ is the bias, $\beta_{1_x}$ is the sensitivity, $\beta_{1_y}$ is the coefficient for the cross-axis sensitivity of $g_y$, $\beta_{1_z}$ is the coefficient for the cross-axis sensitivity of $g_z$. The sensitivity matrix, $\mathbf{\beta}$, and the bias vector, $\mathbf{\beta_0}$, are given by 

$$
\mathbf{\beta} = \left(\begin{array}{rrr}
\beta_{x_x} & \beta_{x_y} & \beta_{x_z} \\
\beta_{y_x} & \beta_{y_y} & \beta_{y_z} \\
\beta_{z_x} & \beta_{z_y} & \beta_{z_z} \\
\end{array}\right), \quad \mathbf{\beta_0} = \left(\begin{array}{r} \beta_{0_x} \\ \beta_{0_y}\\ \beta_{1_z} \end{array}\right)
$$


Calibration is done through the `calibration` function which takes matrices $\mathbf{V}$ and $\mathbf{G}$ are arguments and returns the estimates of the sensitivity and bias terms of the model above. Each row of $\mathbf{G}$ should have unit norm. 

$$
\mathbf{V} = \left(\begin{array}{rrr}
v_{1x} &  v_{1y} & v_{1z} \\ 
\vdots & \vdots & \vdots \\ 
v_{nx} &  v_{ny} & v_{nz} 
\end{array}\right), \quad \mathbf{G} = \left(\begin{array}{rrr}
g_{1x} &  g_{1y} & g_{1z} \\ 
\vdots & \vdots & \vdots \\ 
g_{nx} &  g_{ny} & g_{nz}
\end{array}\right)
$$


```r
# create sensitivity matrix and bias vector
Beta0 <- c(0,0,0)
Beta <- diag(3) + matrix(rnorm(9,sd = 1e-4),nrow = 3)

# number of observations
n <- 10
# create matrix of gravity vectors
G <- matrix(rnorm(3*n), nrow =n)
G <- t(apply(G,1,FUN = function(x) x/norm(x,"2")))

# simulate responses based on calibration model with some noise
V <- G %*% Beta + Beta0 + matrix(rnorm(3*n, sd = 1e-2), nrow =n)
```


To fit the calibration model....

```r
cal <- calibration(V,G)
```

Estimates of the sensitivity and bias terms in the model


```r
cal$coef
#> $BetaHat
#>              [,1]          [,2]         [,3]
#> [1,]  1.005119179  0.0052260559 -0.006503218
#> [2,]  0.006232635  1.0089611587  0.009343589
#> [3,] -0.002010252 -0.0006487523  0.996105434
#> 
#> $Beta
#>               [,1]          [,2]          [,3]
#> [1,]  9.999463e-01 -0.0001055828 -1.132899e-04
#> [2,]  2.250167e-04  1.0000872805  7.231114e-05
#> [3,] -3.801792e-05 -0.0001087613  9.999114e-01
```


Estimates of the gravity vectors given by $\hat{g}_i = (v_i - \mathbf{\beta}_0)^T \mathbf{\beta}^{-1}$

```r
cal$Ghat
#>             [,1]       [,2]        [,3]
#>  [1,] -0.2961968 -0.1805839  0.92836277
#>  [2,] -0.8650530 -0.1261805 -0.48201291
#>  [3,] -0.5180814  0.8509841  0.06618255
#>  [4,] -0.8711858 -0.3571561 -0.33160558
#>  [5,] -0.7916695  0.5301957  0.30598846
#>  [6,] -0.9159560  0.2897141 -0.20857129
#>  [7,]  0.5674436 -0.1598501 -0.80851028
#>  [8,] -0.5697024  0.7606461  0.31864690
#>  [9,]  0.7344055 -0.6956871  0.02399297
#> [10,] -0.5959213 -0.4290545  0.69452389
```

Estimate of $\sigma$ from calibration model


```r
cal$sd
#> [1] 0.007640512
```



# Attitude Estimation


Wahba's problem seeks to minimize the cost-function 

$$
L(R) = \frac{1}{n} \sum_{i=1}^n  \| w_i - \mathbf{R} v_i \|^2, \quad n \ge 2
$$
where $\mathbf{w}_{i}$ is the u-th 3-vector measurement in the reference frame, $\mathbf{v}_i$ is the corresponding i-th 3-vector measurement in the body frame $\mathbf{R}$ is a $3 \times 3$ rotation matrix between the reference frame and body frame. The optimal rotation,  $\mathbf{R}_{opt}$ can be found using Davenport's q-method:

Define the corresponding gain function:

```format
\begin{align}
G(R) & = 1- L(R) \\ & = \sum_{i=1}^n tr(\mathbf{w}_i^T \mathbf{R} \mathbf{v}_i) \\ & = 
tr(\mathbf{W}^T \mathbf{R} \mathbf{V}) \\ & = 
tr(\mathbf{R} \mathbf{B}^T)
\end{align}
```

where $\mathbf{B} = \sum_{i=1}^n \mathbf{W} \mathbf{V}^T$

Parameterizing the rotation matrix in terms of quaternion, $\mathbf{q} = (\mathbf{q}_v, q_r)$  we have 

$$
\mathbf{R}(q) = (q_r - \| \mathbf{q}_v \|^2) \mathbf{I}_3 + 2 \mathbf{q}_v \mathbf{q}_v^T - 2 q_r [\mathbf{q}_v]_x
$$

where $[\mathbf{q}_v]_x$ is the skew-symmetric matrix of vector ${q}_v$. The gain function then becomes

$$
g(\mathbf{q}) = (q_r - \| \mathbf{q}_v \|^2) tr(\mathbf{B}^T) + 2 tr (\mathbf{q}_v \mathbf{q}_v^T \mathbf{B}^T) + 2 q_r tr( [\mathbf{q}_v]_x \mathbf{B}^T)
$$
Putting this expression in terms of its bilinear form, we have 

$$
g(\mathbf{q}) = \mathbf{q}^T \mathbf{K} \mathbf{q}
$$

where 

$$
\mathbf{K}_{4 \times 4} = \left[
\begin{array}{rr}
\sigma & \mathbf{z}^T \\
\mathbf{z} & \mathbf{S} - \sigma \mathbf{I}_3
\end{array}
\right]
$$


where 
$$
\sigma = tr(\mathbf{B}), \quad \mathbf{S} = \mathbf{B} + \mathbf{B}^T, \quad \mathbf{z} = \left(
\begin{array}{r}
B_{23} -B{32} \\ 
B_{31} -B{13} \\ 
B_{12} -B{21} \\ 
\end{array}
\right)
$$


The optimal quaternion which maximimizes the gain function is the eigenvector associated with the largest eigenvalue, $\lambda_{max}$, of matrix $\mathbf{K}$ 

$$
\mathbf{K} \mathbf{q}_{opt} = \lambda_{max} \mathbf{q}_{opt}
$$




```r
# create some vector in the reference frame
n <- 3
sd <- 1e-3
W <- matrix(rnorm(3 * n, sd = ),ncol = 3)
W <- t(apply(W,2,FUN = function(x) {x/norm(x, "2")^2}))
W
#>            [,1]       [,2]        [,3]
#> [1,] 0.12888665  0.3189730 -0.25308253
#> [2,] 0.07445148  0.2711011 -0.05165450
#> [3,] 0.06906870 -0.3335964  0.01405484
```



```r
# create some rotation matrix and
pitch <- 90
roll <- 0
yaw <- 0
R <- euler_to_rotation(pitch, roll,yaw, units = "degrees")
# rotation matrix
R
#>      [,1]         [,2]          [,3]
#> [1,]    1 0.000000e+00  0.000000e+00
#> [2,]    0 6.123234e-17 -1.000000e+00
#> [3,]    0 1.000000e+00  6.123234e-17
```




```r
# create vector in body frame and add noise
sd <- 1e-2

V <- t(apply(W,1,function(x) {R %*% x})) + matrix(rnorm(3*n, mean = 0 ,sd = sd),nrow = n)
V
#>            [,1]        [,2]       [,3]
#> [1,] 0.11984215  0.24416946  0.3152047
#> [2,] 0.08777154  0.04547898  0.2824531
#> [3,] 0.07127344 -0.01509245 -0.3121855
```


The optimal rotation is found by `find_rotation`


```r
out <- find_rotation(W,V, output = "euler", sd = NA)
out$rotation
#> $pitch
#> [1] 89.4698
#> 
#> $roll
#> [1] 0.3348523
#> 
#> $yaw
#> [1] 179.8388
```

If standard deviation estimate is given, an estimate of the covariance matrix for the rotation will be returned. 


```r
out <- find_rotation(W,V, output = "euler", sd = sd)
out$rotation
#> $pitch
#> [1] 89.4698
#> 
#> $roll
#> [1] 0.3348523
#> 
#> $yaw
#> [1] 179.8388
out$cov
#>               pitch         roll           yaw
#> pitch  0.0006338298 0.0008534247 -0.0004656837
#> roll   0.0008534247 0.0038904571  0.0004698057
#> yaw   -0.0004656837 0.0004698057  0.0024645931
```

