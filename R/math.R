#' Mean of positive values
#' 
#' Returns the arithmetic mean of values greater than 0.
#' @param x (Numeric Vector) Input values
#' @export mean.pos
mean.pos <- function(x) {
  y <- which(x > 0)
  mean(x[y])
}

#' Arithmetic mean of log-transformed values
#' 
#' Un-logs, takes the arithmetic mean, and re-logs
#' @param x (Numeric vector) Input values
#' @param base (Numeric) Log base to use (default: 2)
#' @export mean.of.logs
mean.of.logs <- function(x, base=2) {
  log(mean((base ^ x) - 1) + 1, base = base)
}

#' Sum of log-transformed values
#' 
#' Un-logs, takes the sum, and re-logs
#' @param x (Numeric vector) Input values
#' @param base (Numeric) Log base to use (default: 2)
#' @export sum.of.logs
sum.of.logs <- function(x, base=2) {
  log(sum((base ^ x) - 1) + 1, base = base)
}

#' Arithmetic mean of positive log-transformed values
#' 
#' Un-logs, takes the arithmetic mean of those values greater than 0, and re-logs
#' @param x (Numeric vector) Input values
#' @param base (Numeric) Log base to use (default: 2)
#' @export mean.of.logs.pos
mean.of.logs.pos <- function(x, base=2) {
  y <- x[x>0]
  if (length(y) > 0) {
    return(log(mean((base ^ y) - 1) + 1, base = base))
  } else {
    return(0)
  }
}

#' Return the signed value with greater magnitude
#' 
#' Determines which value has greater magnitude (i.e. by absolute value) and returns the original signed value.
#' Compares vectors x and y pairwise.
#' @param x (Numeric vector) Input values
#' @param y (Numeric vector) Input values
#' @export pmax.abs
pmax.abs <- function(x, y) {
  z <- y
  x.bigger <- (abs(x) > abs(y))
  z[x.bigger] <- x [x.bigger]
  return(z)
}

#' Returns proportion of expressing cells
#' 
#' Determines proportion of input values > 0.
#' @param x (Numeric vector) Input values
#' @export prop.exp
prop.exp <- function(x) {
  return(length(which(x>0)) / length (x))
}

#' Returns proportion of non-expressing cells
#' 
#' Determines proportion of input values <= 0.
#' @param x (Numeric vector) Input values
#' @export prop.nonexp
prop.nonexp <- function(x) {
  return(length(which(x<=0)) / length (x))
}

#' Logistic function
#' 
#' @param x (Numeric vector) Input values
#' @param x0 (Numeric) Inflection point
#' @param k (Numeric) Slope
#' @param c (Numeric) Max value (default: 1)
logistic <- function(x, x0, k, c=1) {
  c / (1 + exp(-1*k*(x-x0)))
}

#' Inverse logistic function
#' 
#' @param x (Numeric vector) Input values
#' @param x0 (Numeric) Inflection point
#' @param k (Numeric) Slope
#' @param c (Numeric) Max value (default: 1)
inv.logistic <- function(x, x0, k, c=1) {
  log(c/x-1)/(-1*k)+x0
}

#' Determining preference between a pair of values
#' 
#' @param x (Numeric vector) 
#' @param y (Numeric vector) 
#' @param signed (Logical)
preference <- function(x, y, signed=F) {
  z <- as.data.frame(cbind(x, y))
  if (signed) {
    z$p <- apply(z, 1, function(q) (q[1]-q[2])/(q[1]+q[2]))
  } else {
    z$p <- apply(z, 1, function(q) abs(q[1]-q[2])/(q[1]+q[2]))
  }
  # If both x and y are 0, preference will be NA, but want to return 0.
  z[is.na(z$p),"p"] <- 0
  return(z$p)
}

#' Arithmetic Mean of As Numeric
#' 
num.mean <- function(x) mean(as.numeric(x))

sigmoid <- function(x, a, b, c, d) {
  c / (1 + exp(-a*(x-b))) + d
}

d.dx.sigmoid <- function(x, a, b, c, d=NULL) {
  a*c*exp(a*(x+b)) / ((exp(a*x) + exp(a*b))^2)
}

inv.sigmoid <- function(y, a, b, c, d) {
  -log(c/(y-d) - 1)/a + b
}