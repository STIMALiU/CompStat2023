################################
### Computational Statistics ###
### One dimensional Newton   ###
### Example from lecture     ###
################################

# Example function to be maximized
g <- function(x){
  co <- sqrt(2)/4
  1/sqrt(2*3.141592) * exp(-(x-co)^2/2)
}
# First derivative
dg <- function(x){
  co <- sqrt(2)/4
  -(x-co)/sqrt(2*3.141592) * exp(-(x-co)^2/2)
}
# Second derivative
d2g <- function(x){
  co <- sqrt(2)/4
  (-1+(x-co)^2)/sqrt(2*3.141592) * exp(-(x-co)^2/2)
}

# Newton function
# x0 in function-call is starting value
# with graphical output
newton <- function(x0, eps=0.0001){
  xt  <- x0
  xt1 <- x0 + 2  # xt1 is x(t-1); starting value here just to ensure that while-condition is met 
  while(abs(xt-xt1)>eps)
  {
    xt1 <- xt
    xt  <- xt1 - dg(xt1)/d2g(xt1)
    print(xt)   # iteration output
  }
  xt
}

xs   <- -0.2  # starting value - test other values and check what happens
xmax <- newton(xs)

# graphical output
xv   <- seq(-4, 4, by=0.01)
plot(xv, g(xv), type="l", xlab="x", ylab="y")
points(xs, g(xs), col=2, pch=4, lwd=3)      # starting value in plot
points(xmax, g(xmax), col=4, pch=4, lwd=3)  # result in plot
 