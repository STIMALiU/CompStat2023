#########################################################################
### Computational Statistics                                          ###
### Program for two-dimensional steepest ascent example from lectures ###
### Frank Miller, 2023-03-15                                          ###
#########################################################################
# Some constants:
s1  <- 0.6
s2  <- 0.5
mu1 <- 1.5
mu2 <- 1.2

# The function g to be maximised, partial derivatives dg1, dg2, and gradient
g <- function(x)
{
  (exp(-(x[1]^2+x[2]^2)/(2*s1))/s1 + exp(-((x[1]-mu1)^2+(x[2]-mu2)^2)/(2*s2))/s2) / (4*pi)
}
dg1 <- function(x)
{
  ((-x[1]/s1^2)*exp(-(x[1]^2+x[2]^2)/(2*s1)) + (-(x[1]-mu1)/s2^2)*exp(-((x[1]-mu1)^2+(x[2]-mu2)^2)/(2*s2))) / (4*pi)
}
dg2 <- function(x)
{
  ((-x[2]/s1^2)*exp(-(x[1]^2+x[2]^2)/(2*s1)) + (-(x[2]-mu2)/s2^2)*exp(-((x[1]-mu1)^2+(x[2]-mu2)^2)/(2*s2))) / (4*pi)
}
gradient <- function(x)
{
  c(dg1(x), dg2(x))
}
# Second order partial derivatives and Hessian (not needed for steepest ascent)
dg11 <- function(x)
{
  x1 <- x[1]
  x2 <- x[2]
  ((x1^2/s1^3-1/s1^2)*exp(-(x1^2+x2^2)/(2*s1)) + ((x1-mu1)^2/s2^3-1/s2^2)*exp(-((x1-mu1)^2+(x2-mu2)^2)/(2*s2))) / (4*pi)
}
dg12 <- function(x)
{
  x1 <- x[1]
  x2 <- x[2]
  (x1*x2/s1^3*exp(-(x1^2+x2^2)/(2*s1)) + (x1-mu1)*(x2-mu2)/s2^3*exp(-((x1-mu1)^2+(x2-mu2)^2)/(2*s2))) / (4*pi)
}
dg22 <- function(x)
{
  x1 <- x[1]
  x2 <- x[2]
  ((x2^2/s1^3-1/s1^2)*exp(-(x1^2+x2^2)/(2*s1)) + ((x2-mu2)^2/s2^3-1/s2^2)*exp(-((x1-mu1)^2+(x2-mu2)^2)/(2*s2))) / (4*pi)
}
hessian <- function(x)
{
  matrix(c(dg11(x), dg12(x), dg12(x), dg22(x)), nrow=2, ncol=2)
}

# Produce a contour plot; define first a grid where function is evaluated
x1grid <- seq(-2, 2.5, by=0.05)
x2grid <- seq(-2, 3, by=0.05)
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx  <- dx1*dx2
gx  <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
for (j in 1:dx2)
{
  gx[i,j] <- g(c(x1grid[i], x2grid[j]))
}
mgx <- matrix(gx, nrow=dx1, ncol=dx2)
contour(x1grid, x2grid, mgx, nlevels=34)  # Note: For other functions g, you might need to choose another nlevels-value to get a good impression of the function 

#Steepest ascent function:
steepestasc <- function(x0, eps=1e-8, alpha0=1)
{
  xt   <- x0
  conv <- 999
  points(xt[1], xt[2], col=2, pch=4, lwd=3)
  while(conv>eps)
  {
    alpha <- alpha0
    xt1   <- xt
    xt    <- xt1 + alpha*gradient(xt1)
    while (g(xt)<g(xt1))
    {
      alpha <- alpha/2
      xt    <- xt1 + alpha*gradient(xt1)
    }
    points(xt[1], xt[2], col=2, pch=4, lwd=1)
    conv <- sum((xt-xt1)*(xt-xt1))
  }
  points(xt[1], xt[2], col=4, pch=4, lwd=3)
  xt
}
steepestasc(c(1, 0))

