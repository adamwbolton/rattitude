
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rattitude

<!-- badges: start -->

<!-- badges: end -->

rattitude has two main functions: (1) inertial-based three-axis
accelerometer calibration and (2) attitude determination.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("adamwbolton/rattitude")
```

## Example

``` r
library(rattitude)
```

# Rotations

Rotations can be specified either as rotation matrices, sequences of
Euler angles, or rotation quaternions. One can convert between the three
with the following functions: `euler_to_rotation`,
`euler_to_quaternion`, `rotation_to_quaternion`, `rotation_to_euler`,
`quaternion_to_euler`, `quaternion_to_rotation`. Let’s start with a
rotation of 90 degrees about the x-axis and convert this to a rotation
matrix:

``` r
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

``` r
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

``` r
# rotation to quaternion
q <- rotation_to_quaterion(R)
q
#> [1] 0.7071068 0.0000000 0.0000000 0.7071068
```

and back to Euler angles

``` r
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

The calibration for the inertial-based three-axis accelerometers is done
through the following first-order model:

  
![&#10;v\_x = \\beta\_{0\_x} + \\beta\_{x\_x} g\_x + \\beta\_{y\_x} g\_y
+ \\beta\_{z\_x} g\_z +
\\epsilon\_x&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Av_x%20%3D%20%5Cbeta_%7B0_x%7D%20%2B%20%5Cbeta_%7Bx_x%7D%20g_x%20%2B%20%5Cbeta_%7By_x%7D%20g_y%20%2B%20%5Cbeta_%7Bz_x%7D%20g_z%20%2B%20%5Cepsilon_x%0A
"
v_x = \\beta_{0_x} + \\beta_{x_x} g_x + \\beta_{y_x} g_y + \\beta_{z_x} g_z + \\epsilon_x
")  

  
![&#10;v\_y = \\beta\_{0\_y} + \\beta\_{y\_x} g\_x + \\beta\_{y\_y} g\_y
+ \\beta\_{y\_y} g\_z +
\\epsilon\_y&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Av_y%20%3D%20%5Cbeta_%7B0_y%7D%20%2B%20%5Cbeta_%7By_x%7D%20g_x%20%2B%20%5Cbeta_%7By_y%7D%20g_y%20%2B%20%5Cbeta_%7By_y%7D%20g_z%20%2B%20%5Cepsilon_y%0A
"
v_y = \\beta_{0_y} + \\beta_{y_x} g_x + \\beta_{y_y} g_y + \\beta_{y_y} g_z + \\epsilon_y
")  

  
![&#10;v\_z = \\beta\_{0\_z} + \\beta\_{z\_x} g\_x + \\beta\_{z\_y} g\_y
+ \\beta\_{z\_z} g\_z +
\\epsilon\_z&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Av_z%20%3D%20%5Cbeta_%7B0_z%7D%20%2B%20%5Cbeta_%7Bz_x%7D%20g_x%20%2B%20%5Cbeta_%7Bz_y%7D%20g_y%20%2B%20%5Cbeta_%7Bz_z%7D%20g_z%20%2B%20%5Cepsilon_z%0A
"
v_z = \\beta_{0_z} + \\beta_{z_x} g_x + \\beta_{z_y} g_y + \\beta_{z_z} g_z + \\epsilon_z
")  

where ![\\epsilon = (\\epsilon\_x,\\epsilon\_y, \\epsilon\_z)^T \\sim
N(0,\\sigma^2
I\_3)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%20%3D%20%28%5Cepsilon_x%2C%5Cepsilon_y%2C%20%5Cepsilon_z%29%5ET%20%5Csim%20N%280%2C%5Csigma%5E2%20I_3%29
"\\epsilon = (\\epsilon_x,\\epsilon_y, \\epsilon_z)^T \\sim N(0,\\sigma^2 I_3)").
For the x-axis accelerometer,
![\\beta\_{0\_x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7B0_x%7D
"\\beta_{0_x}") is the bias,
![\\beta\_{1\_x}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7B1_x%7D
"\\beta_{1_x}") is the sensitivity,
![\\beta\_{1\_y}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7B1_y%7D
"\\beta_{1_y}") is the coefficient for the cross-axis sensitivity of
![g\_y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g_y
"g_y"),
![\\beta\_{1\_z}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_%7B1_z%7D
"\\beta_{1_z}") is the coefficient for the cross-axis sensitivity of
![g\_z](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g_z
"g_z"). The sensitivity matrix,
![\\mathbf{\\beta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B%5Cbeta%7D
"\\mathbf{\\beta}"), and the bias vector,
![\\mathbf{\\beta\_0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7B%5Cbeta_0%7D
"\\mathbf{\\beta_0}"), are given by

  
![&#10;\\mathbf{\\beta} = \\left(\\begin{array}{rrr}&#10;\\beta\_{x\_x}
& \\beta\_{x\_y} & \\beta\_{x\_z} \\\\&#10;\\beta\_{y\_x} &
\\beta\_{y\_y} & \\beta\_{y\_z} \\\\&#10;\\beta\_{z\_x} & \\beta\_{z\_y}
& \\beta\_{z\_z} \\\\&#10;\\end{array}\\right), \\quad
\\mathbf{\\beta\_0} = \\left(\\begin{array}{r} \\beta\_{0\_x} \\\\
\\beta\_{0\_y}\\\\ \\beta\_{1\_z}
\\end{array}\\right)&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7B%5Cbeta%7D%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Brrr%7D%0A%5Cbeta_%7Bx_x%7D%20%26%20%5Cbeta_%7Bx_y%7D%20%26%20%5Cbeta_%7Bx_z%7D%20%5C%5C%0A%5Cbeta_%7By_x%7D%20%26%20%5Cbeta_%7By_y%7D%20%26%20%5Cbeta_%7By_z%7D%20%5C%5C%0A%5Cbeta_%7Bz_x%7D%20%26%20%5Cbeta_%7Bz_y%7D%20%26%20%5Cbeta_%7Bz_z%7D%20%5C%5C%0A%5Cend%7Barray%7D%5Cright%29%2C%20%5Cquad%20%5Cmathbf%7B%5Cbeta_0%7D%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Br%7D%20%5Cbeta_%7B0_x%7D%20%5C%5C%20%5Cbeta_%7B0_y%7D%5C%5C%20%5Cbeta_%7B1_z%7D%20%5Cend%7Barray%7D%5Cright%29%0A
"
\\mathbf{\\beta} = \\left(\\begin{array}{rrr}
\\beta_{x_x} & \\beta_{x_y} & \\beta_{x_z} \\\\
\\beta_{y_x} & \\beta_{y_y} & \\beta_{y_z} \\\\
\\beta_{z_x} & \\beta_{z_y} & \\beta_{z_z} \\\\
\\end{array}\\right), \\quad \\mathbf{\\beta_0} = \\left(\\begin{array}{r} \\beta_{0_x} \\\\ \\beta_{0_y}\\\\ \\beta_{1_z} \\end{array}\\right)
")  

Calibration is done through the `calibration` function which takes
matrices
![\\mathbf{V}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BV%7D
"\\mathbf{V}") and
![\\mathbf{G}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BG%7D
"\\mathbf{G}") are arguments and returns the estimates of the
sensitivity and bias terms of the model above. Each row of
![\\mathbf{G}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BG%7D
"\\mathbf{G}") should have unit norm.

  
![&#10;\\mathbf{V} = \\left(\\begin{array}{rrr}&#10;v\_{1x} & v\_{1y} &
v\_{1z} \\\\ &#10;\\vdots & \\vdots & \\vdots \\\\ &#10;v\_{nx} &
v\_{ny} & v\_{nz} &#10;\\end{array}\\right), \\quad \\mathbf{G} =
\\left(\\begin{array}{rrr}&#10;g\_{1x} & g\_{1y} & g\_{1z} \\\\
&#10;\\vdots & \\vdots & \\vdots \\\\ &#10;g\_{nx} & g\_{ny} &
g\_{nz}&#10;\\end{array}\\right)&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BV%7D%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Brrr%7D%0Av_%7B1x%7D%20%26%20%20v_%7B1y%7D%20%26%20v_%7B1z%7D%20%5C%5C%20%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%20%0Av_%7Bnx%7D%20%26%20%20v_%7Bny%7D%20%26%20v_%7Bnz%7D%20%0A%5Cend%7Barray%7D%5Cright%29%2C%20%5Cquad%20%5Cmathbf%7BG%7D%20%3D%20%5Cleft%28%5Cbegin%7Barray%7D%7Brrr%7D%0Ag_%7B1x%7D%20%26%20%20g_%7B1y%7D%20%26%20g_%7B1z%7D%20%5C%5C%20%0A%5Cvdots%20%26%20%5Cvdots%20%26%20%5Cvdots%20%5C%5C%20%0Ag_%7Bnx%7D%20%26%20%20g_%7Bny%7D%20%26%20g_%7Bnz%7D%0A%5Cend%7Barray%7D%5Cright%29%0A
"
\\mathbf{V} = \\left(\\begin{array}{rrr}
v_{1x} &  v_{1y} & v_{1z} \\\\ 
\\vdots & \\vdots & \\vdots \\\\ 
v_{nx} &  v_{ny} & v_{nz} 
\\end{array}\\right), \\quad \\mathbf{G} = \\left(\\begin{array}{rrr}
g_{1x} &  g_{1y} & g_{1z} \\\\ 
\\vdots & \\vdots & \\vdots \\\\ 
g_{nx} &  g_{ny} & g_{nz}
\\end{array}\\right)
")  

``` r
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

To fit the calibration model….

``` r
cal <- calibration(V,G)
```

Estimates of the sensitivity and bias terms in the model

``` r
cal$coef
#> $BetaHat
#>             [,1]         [,2]        [,3]
#> [1,] 0.999601350 -0.016203819 0.002389500
#> [2,] 0.002422397  0.991612219 0.004893546
#> [3,] 0.005512963  0.007760368 1.002375200
#> 
#> $Beta
#>               [,1]          [,2]         [,3]
#> [1,]  1.000076e+00  8.059277e-05 1.135893e-04
#> [2,]  5.456315e-05  9.997552e-01 4.076539e-06
#> [3,] -1.038208e-04 -6.119115e-06 9.999143e-01
```

Estimates of the gravity vectors given by ![\\hat{g}\_i = (v\_i -
\\mathbf{\\beta}\_0)^T
\\mathbf{\\beta}^{-1}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bg%7D_i%20%3D%20%28v_i%20-%20%5Cmathbf%7B%5Cbeta%7D_0%29%5ET%20%5Cmathbf%7B%5Cbeta%7D%5E%7B-1%7D
"\\hat{g}_i = (v_i - \\mathbf{\\beta}_0)^T \\mathbf{\\beta}^{-1}")

``` r
cal$Ghat
#>               [,1]        [,2]        [,3]
#>  [1,] -0.708042358  0.59581656  0.38242550
#>  [2,] -0.115753373  0.96577065  0.18156004
#>  [3,]  0.049395167 -1.02113780 -0.01229979
#>  [4,]  0.122635690  0.05509749 -0.98369083
#>  [5,] -0.003190392  0.70095190  0.71947303
#>  [6,] -0.859939419  0.30872718  0.44181834
#>  [7,] -0.256087371 -0.94040855  0.18621783
#>  [8,] -0.927600534 -0.15038225 -0.28208025
#>  [9,]  0.130101505 -0.64665419  0.74144121
#> [10,]  0.971819066  0.27117651  0.07182864
```

Estimate of
![\\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma
"\\sigma") from calibration model

``` r
cal$sd
#> [1] 0.01080902
```

# Attitude Estimation

Wahba’s problem seeks to minimize the cost-function

  
![&#10;L(R) = \\frac{1}{n} \\sum\_{i=1}^n \\| w\_i - \\mathbf{R} v\_i
\\|^2, \\quad n
\\ge 2&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AL%28R%29%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%20%5Csum_%7Bi%3D1%7D%5En%20%20%5C%7C%20w_i%20-%20%5Cmathbf%7BR%7D%20v_i%20%5C%7C%5E2%2C%20%5Cquad%20n%20%5Cge%202%0A
"
L(R) = \\frac{1}{n} \\sum_{i=1}^n  \\| w_i - \\mathbf{R} v_i \\|^2, \\quad n \\ge 2
")  
where
![\\mathbf{w}\_{i}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bw%7D_%7Bi%7D
"\\mathbf{w}_{i}") is the u-th 3-vector measurement in the reference
frame,
![\\mathbf{v}\_i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bv%7D_i
"\\mathbf{v}_i") is the corresponding i-th 3-vector measurement in the
body frame
![\\mathbf{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D
"\\mathbf{R}") is a ![3
\\times 3](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3%20%5Ctimes%203
"3 \\times 3") rotation matrix between the reference frame and body
frame. The optimal rotation,
![\\mathbf{R}\_{opt}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BR%7D_%7Bopt%7D
"\\mathbf{R}_{opt}") can be found using Davenport’s q-method:

Define the corresponding gain function:

``` format
\begin{align}
G(R) & = 1- L(R) \\ & = \sum_{i=1}^n tr(\mathbf{w}_i^T \mathbf{R} \mathbf{v}_i) \\ & = 
tr(\mathbf{W}^T \mathbf{R} \mathbf{V}) \\ & = 
tr(\mathbf{R} \mathbf{B}^T)
\end{align}
```

where ![\\mathbf{B} = \\sum\_{i=1}^n \\mathbf{W}
\\mathbf{V}^T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BB%7D%20%3D%20%5Csum_%7Bi%3D1%7D%5En%20%5Cmathbf%7BW%7D%20%5Cmathbf%7BV%7D%5ET
"\\mathbf{B} = \\sum_{i=1}^n \\mathbf{W} \\mathbf{V}^T")

Parameterizing the rotation matrix in terms of quaternion, ![\\mathbf{q}
= (\\mathbf{q}\_v,
q\_r)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7Bq%7D%20%3D%20%28%5Cmathbf%7Bq%7D_v%2C%20q_r%29
"\\mathbf{q} = (\\mathbf{q}_v, q_r)") we have

  
![&#10;\\mathbf{R}(q) = (q\_r - \\| \\mathbf{q}\_v \\|^2) \\mathbf{I}\_3
+ 2 \\mathbf{q}\_v \\mathbf{q}\_v^T - 2 q\_r
\[\\mathbf{q}\_v\]\_x&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BR%7D%28q%29%20%3D%20%28q_r%20-%20%5C%7C%20%5Cmathbf%7Bq%7D_v%20%5C%7C%5E2%29%20%5Cmathbf%7BI%7D_3%20%2B%202%20%5Cmathbf%7Bq%7D_v%20%5Cmathbf%7Bq%7D_v%5ET%20-%202%20q_r%20%5B%5Cmathbf%7Bq%7D_v%5D_x%0A
"
\\mathbf{R}(q) = (q_r - \\| \\mathbf{q}_v \\|^2) \\mathbf{I}_3 + 2 \\mathbf{q}_v \\mathbf{q}_v^T - 2 q_r [\\mathbf{q}_v]_x
")  

where
![\[\\mathbf{q}\_v\]\_x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5B%5Cmathbf%7Bq%7D_v%5D_x
"[\\mathbf{q}_v]_x") is the skew-symmetric matrix of vector
![{q}\_v](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7Bq%7D_v
"{q}_v"). The gain function then becomes

  
![&#10;g(\\mathbf{q}) = (q\_r - \\| \\mathbf{q}\_v \\|^2)
tr(\\mathbf{B}^T) + 2 tr (\\mathbf{q}\_v \\mathbf{q}\_v^T \\mathbf{B}^T)
+ 2 q\_r tr( \[\\mathbf{q}\_v\]\_x
\\mathbf{B}^T)&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ag%28%5Cmathbf%7Bq%7D%29%20%3D%20%28q_r%20-%20%5C%7C%20%5Cmathbf%7Bq%7D_v%20%5C%7C%5E2%29%20tr%28%5Cmathbf%7BB%7D%5ET%29%20%2B%202%20tr%20%28%5Cmathbf%7Bq%7D_v%20%5Cmathbf%7Bq%7D_v%5ET%20%5Cmathbf%7BB%7D%5ET%29%20%2B%202%20q_r%20tr%28%20%5B%5Cmathbf%7Bq%7D_v%5D_x%20%5Cmathbf%7BB%7D%5ET%29%0A
"
g(\\mathbf{q}) = (q_r - \\| \\mathbf{q}_v \\|^2) tr(\\mathbf{B}^T) + 2 tr (\\mathbf{q}_v \\mathbf{q}_v^T \\mathbf{B}^T) + 2 q_r tr( [\\mathbf{q}_v]_x \\mathbf{B}^T)
")  
Putting this expression in terms of its bilinear form, we have

  
![&#10;g(\\mathbf{q}) = \\mathbf{q}^T \\mathbf{K}
\\mathbf{q}&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Ag%28%5Cmathbf%7Bq%7D%29%20%3D%20%5Cmathbf%7Bq%7D%5ET%20%5Cmathbf%7BK%7D%20%5Cmathbf%7Bq%7D%0A
"
g(\\mathbf{q}) = \\mathbf{q}^T \\mathbf{K} \\mathbf{q}
")  

where

  
![&#10;\\mathbf{K}\_{4 \\times 4} =
\\left\[&#10;\\begin{array}{rr}&#10;\\sigma & \\mathbf{z}^T
\\\\&#10;\\mathbf{z} & \\mathbf{S} - \\sigma
\\mathbf{I}\_3&#10;\\end{array}&#10;\\right\]&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BK%7D_%7B4%20%5Ctimes%204%7D%20%3D%20%5Cleft%5B%0A%5Cbegin%7Barray%7D%7Brr%7D%0A%5Csigma%20%26%20%5Cmathbf%7Bz%7D%5ET%20%5C%5C%0A%5Cmathbf%7Bz%7D%20%26%20%5Cmathbf%7BS%7D%20-%20%5Csigma%20%5Cmathbf%7BI%7D_3%0A%5Cend%7Barray%7D%0A%5Cright%5D%0A
"
\\mathbf{K}_{4 \\times 4} = \\left[
\\begin{array}{rr}
\\sigma & \\mathbf{z}^T \\\\
\\mathbf{z} & \\mathbf{S} - \\sigma \\mathbf{I}_3
\\end{array}
\\right]
")  

where   
![&#10;\\sigma = tr(\\mathbf{B}), \\quad \\mathbf{S} = \\mathbf{B} +
\\mathbf{B}^T, \\quad \\mathbf{z} =
\\left(&#10;\\begin{array}{r}&#10;B\_{23} -B{32} \\\\ &#10;B\_{31}
-B{13} \\\\ &#10;B\_{12} -B{21} \\\\
&#10;\\end{array}&#10;\\right)&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Csigma%20%3D%20tr%28%5Cmathbf%7BB%7D%29%2C%20%5Cquad%20%5Cmathbf%7BS%7D%20%3D%20%5Cmathbf%7BB%7D%20%2B%20%5Cmathbf%7BB%7D%5ET%2C%20%5Cquad%20%5Cmathbf%7Bz%7D%20%3D%20%5Cleft%28%0A%5Cbegin%7Barray%7D%7Br%7D%0AB_%7B23%7D%20-B%7B32%7D%20%5C%5C%20%0AB_%7B31%7D%20-B%7B13%7D%20%5C%5C%20%0AB_%7B12%7D%20-B%7B21%7D%20%5C%5C%20%0A%5Cend%7Barray%7D%0A%5Cright%29%0A
"
\\sigma = tr(\\mathbf{B}), \\quad \\mathbf{S} = \\mathbf{B} + \\mathbf{B}^T, \\quad \\mathbf{z} = \\left(
\\begin{array}{r}
B_{23} -B{32} \\\\ 
B_{31} -B{13} \\\\ 
B_{12} -B{21} \\\\ 
\\end{array}
\\right)
")  

The optimal quaternion which maximimizes the gain function is the
eigenvector associated with the largest eigenvalue,
![\\lambda\_{max}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_%7Bmax%7D
"\\lambda_{max}"), of matrix
![\\mathbf{K}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BK%7D
"\\mathbf{K}")

  
![&#10;\\mathbf{K} \\mathbf{q}\_{opt} = \\lambda\_{max}
\\mathbf{q}\_{opt}&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BK%7D%20%5Cmathbf%7Bq%7D_%7Bopt%7D%20%3D%20%5Clambda_%7Bmax%7D%20%5Cmathbf%7Bq%7D_%7Bopt%7D%0A
"
\\mathbf{K} \\mathbf{q}_{opt} = \\lambda_{max} \\mathbf{q}_{opt}
")  

``` r
# create some vector in the reference frame
n <- 3
sd <- 1e-3
W <- matrix(rnorm(3 * n, sd = ),ncol = 3)
W <- t(apply(W,2,FUN = function(x) {x/norm(x, "2")^2}))
W
#>            [,1]        [,2]        [,3]
#> [1,] -0.3876940 -0.08612243 -0.42998760
#> [2,] -0.2685413  0.47500568  0.03848487
#> [3,] -0.2655570  0.69677187 -0.32195384
```

``` r
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

``` r
# create vector in body frame and add noise
sd <- 1e-2

V <- t(apply(W,1,function(x) {R %*% x})) + matrix(rnorm(3*n, mean = 0 ,sd = sd),nrow = n)
V
#>            [,1]        [,2]        [,3]
#> [1,] -0.3791908  0.43733015 -0.07604239
#> [2,] -0.2631449 -0.03920862  0.48399961
#> [3,] -0.2604135  0.31035092  0.68793813
```

The optimal rotation is found by `find_rotation`

``` r
out <- find_rotation(W,V, output = "euler", sd = NA)
out$rotation
#> $pitch
#> [1] 89.26877
#> 
#> $roll
#> [1] 179.7568
#> 
#> $yaw
#> [1] 179.2509
```

If standard deviation estimate is given, an estimate of the covariance
matrix for the rotation will be returned.

``` r
out <- find_rotation(W,V, output = "euler", sd = sd)
out$rotation
#> $pitch
#> [1] 89.26877
#> 
#> $roll
#> [1] 179.7568
#> 
#> $yaw
#> [1] 179.2509
out$cov
#>            pitch       roll        yaw
#> pitch 0.01583568 0.01204993 0.02902866
#> roll  0.01204993 0.01242081 0.02620081
#> yaw   0.02902866 0.02620081 0.05957948
```
