################################
### Computational Statistics ###
### Code from math-lecture 1 ###
################################

### Analytic optimization    ###
g <- function(x1, x2){
  -3*x1^2 - 4*x2^2 + x1*x2^3
}
# symbolic partial derivatives
D(expression(-3*x1^2 - 4*x2^2 + x1*x2^3), "x1")
D(expression(-3*x1^2 - 4*x2^2 + x1*x2^3), "x2")

# prepare contour/3d-plot
x1grid <- -20:50/20
x2grid <- -20:50/20
dx1 <- length(x1grid)
dx2 <- length(x2grid)
dx <- dx1*dx2

gx <- matrix(rep(NA, dx), nrow=dx1)
for (i in 1:dx1)
for (j in 1:dx2)
{
  x <- x1grid[i]
  y <- x2grid[j]
  gx[i, j] <- g(x, y)
}
mgx <- matrix(gx, nrow=dx1, ncol=dx2)

# 3d plot
mgxcut <- mgx/(mgx > -15)   # set values below -15 to -Inf to avoid them in the plot
persp(x1grid, x2grid, mgxcut, xlab="x", ylab="y", zlab="z", theta=45, phi=20, zlim=c(-15, 0))

dev.new(width=6, height=6, unit="cm")
par(oma=c(1, 1, 0, 1), mar=c(2, 3, 0.5, 0))
contour(x1grid, x2grid, mgx, nlevels=20)

H <- function(x){
  matrix(c(-6, 3*x[2]^2, 3*x[2]^2, -8+6*x[1]*x[2]), ncol=2)
}

cand1  <- c(0, 0)
ecand1 <- eigen(H(cand1))
ecand1
cand2  <- c(4/3, 2)
ecand2 <- eigen(H(cand2))
ecand2

# draw eigenvectors in contour plot
line1 <- rbind(cand2 - ecand2$vectors[, 1], cand2 + ecand2$vectors[, 1])
line2 <- rbind(cand2 - ecand2$vectors[, 2], cand2 + ecand2$vectors[, 2])
lines(line1, lwd=5, col=4)
lines(line2, lwd=5, col=2)


### Basic matrix algebra     ###
X <- matrix(c(1, 2,
              1, 5,
              1, 8), ncol=2, byrow=TRUE)

V    <- t(X) %*% X
Vinv <- solve(V)
Vinv
Vinv %*% V

A <- matrix(c(1, 0, 0, 0), ncol=2)
solve(A)

# Example how the LS estimator can be coded
y    <- c(2, 2, 4)
beta <- solve(t(X) %*% X) %*% t(X) %*% y
beta <- solve(crossprod(X)) %*% t(X) %*% y
beta <- solve(crossprod(X), t(X) %*% y)

# investigate running time for 10^5 LS-calculations
stime <- proc.time()
for (i in 1:1e5) beta <- solve(t(X) %*% X) %*% t(X) %*% y
proc.time()-stime

stime <- proc.time()
for (i in 1:1e5) beta <- solve(crossprod(X)) %*% t(X) %*% y
proc.time()-stime

stime <- proc.time()
for (i in 1:1e5) beta <- solve(crossprod(X), t(X) %*% y)
proc.time()-stime
# around 40% saving in running time


### About determinants       ###
# density of multivariate normal distribution
f <- function(x, mu, Sigma){
  p <- length(mu)
  1/(2*pi^p*det(Sigma)) * exp(-1/2 * t(x-mu) %*% solve(Sigma) %*% (x-mu))
}
f(c(0, 0), c(0, 0), diag(2))   # Sigma = identity matrix I = diag(2)

# det(solve(A))=1/det(A), the latter is faster to compute
set.seed(2021)
size <- 5
rep <- 100000
A <- matrix(rnorm(size^2), ncol=size)
det(solve(A))
#[1] -0.5090438
1/det(A)
#[1] -0.5090438

stime <- proc.time()
for (i in 1:rep){
  s <- det(solve(A))
}
proc.time() - stime
#user system elapsed
#5.00 0.06 5.06

stime <- proc.time()
for (i in 1:rep){
  s <- 1/det(A)
}
proc.time() - stime
#user system elapsed
#1.52 0.03 1.55

