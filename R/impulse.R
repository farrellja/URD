#' Log-likelihood ratio test
#'
#' @param D (Numeric)
#' @param M0 (Numeric)
#' @param M (Numeric)
#' @param sd.bg (Numeric) Standard deviation of background genes for calculating null
log.likelihood.ratio <- function(D, M0, M, sd.bg) {
  y0 <- (D - M0) / sd.bg # p(D|M0)
  p0 <- pnorm(y0)
  y <- (D - M) / sd.bg # p(D|M)
  p <- pnorm(y)
  
  e <- sum(log2(p) - log2(p0))
  return(e)
}

#' Evaluate the impulse model (double sigmoid) function, given parameters
#' 
#' Chechik et al. 2008 Nat Biotech, with Michal's help
#'
#' @param x (Numeric)
#' @param b1 (Numeric) slope of onset (should be negative for impulse)
#' @param b2 (Numeric) slope of offset (should be positive for impulse)
#' @param h0 (Numeric) initial level
#' @param h1 (Numeric) peak level
#' @param h2 (Numeric) after offset level ("steady-state")
#' @param t1 (Numeric) time of onset
#' @param t2 (Numeric) time of offset
#' 
#' @return Double-sigmoid impulse function with parameters b1, b2, h0, h1, h2, t1, t2, as evaluated for the vector x.
impulse.double <- function(x, b1, b2, h0, h1, h2, t1, t2) {
  f1 <- h0 + (h1 - h0)/(1 + exp(b1*(x - t1)))
  f2 <- h2 + (h1 - h2)/(1 + exp(b2*(x - t2)))
  sg1 <- f1*f2/h1
  return(sg1)
}

#' Evaluate the impulse model (single sigmoid) function, given parameters
#' 
#' Chechik et al. 2008 Nat Biotech, with Michal's help
#'
#' @param x (Numeric)
#' @param b1 (Numeric) slope of onset (should be negative for onset)
#' @param h0 (Numeric) initial level
#' @param h1 (Numeric) after onset level ("steady-state")
#' @param t1 (Numeric) time of onset
#' 
#' @return Double-sigmoid impulse function with parameters b1, b2, h0, h1, h2, t1, t2, as evaluated for the vector x.
impulse.single <- function(x, b1, h0, h1, t1) {
  h0 + (h1 - h0)/(1 + exp(b1*(x - t1)))
}

#' Starting Conditions for impulse model (double sigmoid)
#' 
#' Function to determine k sets of initial starting conditions for impulse fitting
#' 
impulse.start.double <- function(x, y, k, limit.shape=c("none","concave","convex")) {
  df <- data.frame(x = x, y = y)
  
  # Estimate something 'intelligent' for h0, h1, h2
  tenth <- ceiling(length(y)/10)
  start.tenth <- mean(y[1:tenth])
  end.tenth <- mean(y[length(y)-tenth:length(y)])
  if (limit.shape[1]=="concave") {
    poss.h <- sort(runif(3*(k-1), min(y), max(y)), decreasing = F)
    h1 <- c(0, poss.h[1:(k-1)])
    poss.h <- sample(sort(poss.h)[1:(2*(k-1))], replace=F)
    h0 <- c(start.tenth, poss.h[1:(k-1)])
    h2 <- c(end.tenth, poss.h[k:(2*k-2)])
  } else {
    poss.h <- sort(runif(3*(k-1), min(y), max(y)), decreasing = T)
    h1 <- c(1, poss.h[1:(k-1)])
    poss.h <- sample(sort(poss.h)[1:(2*(k-1))], replace=F)
    h0 <- c(start.tenth, poss.h[1:(k-1)])
    h2 <- c(end.tenth, poss.h[k:(2*k-2)])
  }
  
  # Estimate something intelligent for t1, t2
  if (limit.shape[1]=="concave") inflect <- x[which.min(y)] else inflect <- x[which.max(y)]
  middling <- x[which(abs(y-mean(y)) < 0.1)]
  t20 <- mean(middling[middling > inflect])
  t10 <- mean(middling[middling < inflect])
  possiblet <- sort(runif(2*k-2, min=min(x), max=max(x)))
  t1 <- c(t10, possiblet[1:k-1])
  t2 <- c(t20, possiblet[k:(2*k-2)])
  
  # Estimate b based on the other random stuff?
  starts <- data.frame(
    h0=h0,
    h1=h1,
    h2=h2,
    t1=t1,
    t2=t2,
    stringsAsFactors=F
  )
  starts$b1 <- (starts$h0-starts$h1)/(apply(starts[,c("t1","t2")], 1, mean)-starts$t1)
  starts$b2 <- (starts$h1-starts$h2)/(starts$t2 - apply(starts[,c("t1","t2")], 1, mean))
  
  # Reorder to the correct variable order.
  starts <- starts[,c("b1","b2","h0","h1","h2","t1","t2")]
  return(starts)
}

#' Starting Conditions for impulse model (single sigmoid)
#' 
#' Function to determine k sets of initial starting conditions for impulse fitting
#' 
impulse.start.single <- function(x, y, k, limit.slope) {
  ## Estimate two 'intelligent' pairs for h0, h1 from the data
  # Mean expression across first and last tenth of pseudotime
  mean.first.tenth.pt <- mean(y[which(x <= (range(x)[1] + 0.1*diff(range(x))))])
  mean.last.tenth.pt <- mean(y[which(x >= (range(x)[2] - 0.1*diff(range(x))))])
  # Mean expression across lowest 5% of data and highest 5% of data
  mean.lowest.data <- mean(sort(y)[1:round(length(y)/20)])
  mean.highest.data <- mean(sort(y, decreasing = T)[1:round(length(y)/20)])
  
  ## Estimate inflection points for t1, based on where expression has middle value
  # Find all stretches of x that have scaled expression between 0.33 & 0.67,
  # and try those as inflection points.
  middle.rle <- rle(y > 0.33 & y < 0.67)
  middle.rle.df <- data.frame(
    length=middle.rle$lengths,
    start=cumsum(c(0,head(middle.rle$lengths,-1)))+1,
    end=cumsum(middle.rle$lengths),
    value=middle.rle$values,
    stringsAsFactors=F
  )
  t10 <- apply(middle.rle.df[middle.rle.df$value==T,], 1, function(i) mean(x[i[2]:i[3]]))

  ## Estimate slopes
  ## Estimate a b1 that matches slope within inflection stretches
  b10 <- -1*apply(middle.rle.df[middle.rle.df$value==T,], 1, function(i) (y[i[3]]-y[i[2]]) / (x[i[3]] - x[i[2]]))

  # Make a data.frame of potential starting conditions
  l <- length(t10)
  starts <- data.frame(
    h0=c(rep(mean.first.tenth.pt, 2*l), rep(mean.lowest.data, 2*l), runif(k, min(y), max(y))),
    h1=c(rep(mean.last.tenth.pt, 2*l), rep(mean.highest.data, 2*l), runif(k, min(y), max(y))),
    t1=c(rep(t10, 4), runif(k, min(x), max(x))),
    b1=c(b10, rep(NA, l), b10, rep(NA, l), rep(NA, k)),
    stringsAsFactors=F
  )
  l <- length(which(is.na(starts$b1)))
  starts[is.na(starts$b1), "b1"] <- seq(1,100,length.out=l)*with(starts[is.na(starts$b1),], (h1-h0) / (min(x)-max(x)))
  
  # If limit slope, then reverse all the ones that need reversing
  if (limit.slope=="on") {
    starts.2 <- as.data.frame(t(apply(starts, 1, function(x) {
      if (x[1] > x[2]) y <- c(x[2], x[1]) else y <- c(x[1], x[2]) # Make sure h0 < h1
      if (x[4] > 0) z <- -1 * x[4] else z <- x[4] # Make sure slope negative (because onset slopes are neg)
      return(c(y, x[3], z))
    })))
  } else if (limit.slope=="off") {
    starts.2 <- as.data.frame(t(apply(starts, 1, function(x) {
      if (x[1] < x[2]) y <- c(x[2], x[1]) else y <- c(x[1], x[2]) # Make sure h0 > h1
      if (x[4] < 0) z <- -1 * x[4] else z <- x[4] # Make sure slope positive (because onset slopes are neg)
      return(c(y, x[3], z))
    })))
  } else {
    starts.2 <- starts
  }
  names(starts.2) <- names(starts)
  
  # Reorder to match the right variable order.
  starts.2 <- starts.2[,c("b1","h0","h1","t1")]
  return(starts)
}

#' Fit Impulse Model (Double sigmoid)
#' 
#' @importFrom minpack.lm nlsLM nls.lm.control
#' 
#' @param x (Numeric)
#' @param y (Numeric)
#' @param k (Numeric) Number of starting conditions to try
#' @param limit.shape ("none", "convex", "concave") 
impulse.fit.double <- function(x, y, k=20, limit.shape=c("none","convex","concave")) {
  # Determine parameter limits
  if (limit.shape[1] == "convex") {
    b1.low=-Inf
    b1.high=0
    b2.low=0
    b2.high=Inf
  } else if (limit.shape[1]=="concave") {
    b1.low=0
    b1.high=Inf
    b2.low=-Inf
    b2.high=0
  } else {
    b1.low=-Inf
    b1.high=Inf
    b2.low=-Inf
    b2.high=Inf
  }
  df <- data.frame(x = x, y = y)
  # Get starting conditions
  starts <- impulse.start.double(x, y, k, limit.shape = limit.shape)
  # Try the fits!
  con <- nls.lm.control(maxiter=200)
  fits <- lapply(1:k, function(f) try(nlsLM(y ~ impulse.double(x, b1, b2, h0, h1, h2, t1, t2), data=df, start=starts[f,], lower=c(b1.low, b2.low, 0, 0, 0, -Inf, -Inf), upper=c(b1.high, b2.high, 1, 1, 1, Inf, Inf), control=con), silent=T))
  # Figure out the error of the fits and choose the best one.
  best.fit <- which.min(unlist(lapply(fits, function(f) {
    if (class(f) == "try-error") return(NA)
    sum(resid(f)^2)
  })))
  # Return values
  if (length(best.fit)==0) return(rep(NA, 8))
  to.ret <- c(coef(fits[[best.fit]]), sum(resid(fits[[best.fit]])^2))
  names(to.ret)[8] <- "sum.resid2"
  return(to.ret)
}

#' Fit Impulse Model (Single sigmoid)
#' 
#' @importFrom minpack.lm nlsLM nls.lm.control
#' 
#' @param x (Numeric)
#' @param y (Numeric)
#' @param k (Numeric) Number of starting conditions to try
#' @param limit.slope ("none", "on", "off") 
impulse.fit.single <- function(x, y, k=20, limit.slope=c("none","on","off")) {
  if (length(limit.slope) > 1) limit.slope <- limit.slope[1]
  if (!(limit.slope %in% c("none","on","off"))) stop("limit.slope must be 'none', 'on', or 'off'.")
  df <- data.frame(x = x, y = y)
  # Get starting conditions
  starts <- impulse.start.single(x, y, k, limit.slope=limit.slope)
  # Limit slopes?
  if (limit.slope == "on") {
    slope.low=-Inf
    slope.high=0
  } else if (limit.slope == "off") {
    slope.low=0
    slope.high=Inf
  } else {
    slope.low=-Inf
    slope.high=Inf
  }
  
  # Try the fits!
  con <- nls.lm.control(maxiter=200)
  fits <- lapply(1:nrow(starts), function(f) try(nlsLM(y ~ impulse.single(x, b1, h0, h1, t1), data=df, start=starts[f,], control=con, lower=c(slope.low, 0, 0, -Inf), upper=c(slope.high, 1, 1, Inf)), silent=T))
  # Figure out the error of the fits and choose the best one.
  best.fit <- which.min(unlist(lapply(fits, function(f) {
    if (class(f) == "try-error") return(NA)
    sum(resid(f)^2)
  })))
  # Return values
  if (length(best.fit)==0) return(rep(NA, 5))
  to.ret <- c(coef(fits[[best.fit]]), sum(resid(fits[[best.fit]])^2))
  names(to.ret)[5] <- "sum.resid2"
  return(to.ret)
}

#' Log likelihood ratio test for impulse fits
#' 
#' Determines whether additional degrees of freedom in double sigmoid impulse fit are
#' warranted by improved fit to the data
#' 
#' @param x (Numeric) x-values
#' @param y (Numeric)
#' @param I1 (Numeric)
#' @param I2 (Numeric)
#' @param sd.bg (Numeric) Standard deviation of background data for estimating noise model
#' @param df (Numeric) Difference in degrees of freedom (i.e. number of parameters)
impulse.llrtest<- function(x, y, I1, I2, sd.bg, df=3) {
  M0 <- impulse.single(x, b1=as.numeric(I1['b1']), h0=as.numeric(I1['h0']), h1=as.numeric(I1['h1']), t1=as.numeric(I1['t1']))
  M <- impulse.double(x, b1=as.numeric(I2['b1']), b2=as.numeric(I2['b2']), h0=as.numeric(I2['h0']), h1=as.numeric(I2['h1']), h2=as.numeric(I2['h2']), t1=as.numeric(I2['t1']), t2=as.numeric(I2['t2']))
  
  e <- log.likelihood.ratio(y, M0, M, sd.bg)
  p <- 1 - pchisq(2*log(2)*e, df)
  return(p)
}

#' Log likelihood ratio test for impulse fits
#' 
#' Determines whether additional degrees of freedom in double sigmoid impulse fit are
#' warranted by improved fit to the data
#' 
#' @param y (Numeric) Values of fit data
#' @param M0 (Numeric) Values from fit function M0 (null)
#' @param M (Numeric) Values from fit function M (test)
#' @param sd.bg (Numeric) Standard deviation of background data for estimating noise model
#' @param df (Numeric) Difference in degrees of freedom (i.e. number of parameters) between M and M0.
llrtest.dof <- function(y, M0, M, sd.bg, df) {
  e <- log.likelihood.ratio(y, M0, M, sd.bg)
  p <- 1 - pchisq(2*log(2)*e, df)
  return(p)
}

#' Fit gene expression data with an impulse model
#' 
#' @param x (Numeric)
#' @param y (Numeric)
#' @param limit.single.slope ("none", "on", "off")
#' @param sd.bg (Numeric) Standard deviation of background data
#' @param a (Numeric) 
#' @param k (Numeric) Number of sets of starting conditions to try
#' @param onset.thresh (Numeric)
#' 
#' @export
impulseFit <- function(x, y, limit.single.slope=c("none","on","off"), sd.bg, a=0.05, k=50, onset.thresh=0.1) {
  ## Fit linear regression, single sigmoid, and double sigmoids (concave and convex)
  linear <- lm(y~x, data=list(x=x, y=y))
  linear <- c(coef(linear), sum(resid(linear)^2))
  names(linear) <- c("b", "a", "sum.resid2")
  single <- impulse.fit.single(x=x, y=y, k=k, limit.slope = limit.single.slope)
  double.convex <- impulse.fit.double(x=x, y=y, k=k, limit.shape = "convex")
  double.concave <- impulse.fit.double(x=x, y=y, k=k, limit.shape = "concave")
  
  # Which fits are options?
  option.linear <- T # Linear regression will always work
  if (any(is.na(single))) option.single <- F else option.single <- T
  
  # Choose between the double fits - if only one worked, it moves on, otherwise choose best error.
  if (any(is.na(double.convex)) & any(is.na(double.concave))) {
    option.double <- F
    double <- NULL
  } else if (any(is.na(double.convex))) {
    double <- c(type=3, double.concave)
    option.double <- T
  } else if (any(is.na(double.concave))) {
    double <- c(type=2, double.convex)
    option.double <- T
  } else if (double.concave['sum.resid2'] < double.convex['sum.resid2']) {
    double <- c(type=3, double.concave)
    option.double <- T
  } else {
    double <- c(type=2, double.convex)
    option.double <- T
  }
  
  ## Choose the best fit
  
  # If single sigmoid, exists, check whether it fits better than linear model
  if (option.single) {
    M0 <- linear["b"] + linear["a"]*x
    M <- impulse.single(x, b1=as.numeric(single['b1']), h0=as.numeric(single['h0']), h1=as.numeric(single['h1']), t1=as.numeric(single['t1']))
    llr <- llrtest.dof(y, M0=M0, M=M, sd.bg = sd.bg, df=2)
    if (llr > a) option.single <- F
  }
  # Check whether double sigmoid fits better than the currently best model (single sigmoid or linear)
  if (option.double) {
    if (option.single) {
      # If double sigmoid exists and single sigmoid exists (and is better than linear), check whether double fits better than single
      M0 <- impulse.single(x, b1=as.numeric(single['b1']), h0=as.numeric(single['h0']), h1=as.numeric(single['h1']), t1=as.numeric(single['t1']))
      M <- impulse.double(x, b1=as.numeric(double['b1']), b2=as.numeric(double['b2']), h0=as.numeric(double['h0']), h1=as.numeric(double['h1']), h2=as.numeric(double['h2']), t1=as.numeric(double['t1']), t2=as.numeric(double['t2']))
      llr <- llrtest.dof(y, M0=M0, M=M, sd.bg = sd.bg, df=3)
      if (llr > a) option.double <- F
    } else {
      # Single sigmoid not good, but maybe double is still better than linear
      M0 <- linear["b"] + linear["a"]*x
      M <- impulse.double(x, b1=as.numeric(double['b1']), b2=as.numeric(double['b2']), h0=as.numeric(double['h0']), h1=as.numeric(double['h1']), h2=as.numeric(double['h2']), t1=as.numeric(double['t1']), t2=as.numeric(double['t2']))
      llr <- llrtest.dof(y, M0=M0, M=M, sd.bg = sd.bg, df=5)
      if (llr > a) option.double <- F
    }
  }
 
  if (option.double) {
    best.fit <- "double"
  } else if (option.single) {
    best.fit <- "single"
  } else if (option.linear) {
    best.fit <- "linear"
  }
  
  ## Figure out time of onset and offset
  if (best.fit=="linear") {
    if (linear["a"] > 0) {
      # Positive slope
      time.off <- Inf
      time.on <- -1*linear["b"]/linear["a"]
    } else {
      # Negative slope
      time.on <- -Inf
      time.off <- -1*linear["b"]/linear["a"]
    }
  } else if (best.fit=="single") {
    y.pred <- impulse.single(x, b1 = single['b1'], h0 = single['h0'], h1 = single['h1'], t1=single['t1'])
    onset <- abs((single['h1'] - single['h0'])) * onset.thresh + min(single['h0'], single['h1'])
    # Single, onset
    if (single['b1'] < 0) {
      time.on <- x[min(which(y.pred > onset))]
      time.off <- Inf
    } else {
      # Single, offset
      time.on <- -Inf
      time.off <- x[min(which(y.pred < onset))]
    }
  } else if (best.fit=="double") {
    y.pred <- impulse.double(x, b1 = double['b1'], b2=double['b2'], h0 = double['h0'], h1 = double['h1'], h2=double['h2'], t1=double['t1'], t2=double['t2'])
    # Is it convex or concave? (You still can't trust the previous function to fit them correctly.)
    outer <- mean(c(head(y.pred, 5), tail(y.pred, 5)))
    inner <- mean(y.pred[ceiling((length(y.pred)/2) + c(-5,5))])
    if (inner > outer) {
      # Double, convex
      init <- mean(y.pred[1:5])
      top <- sort(y.pred, decreasing=T)[5]
      onset <- (top-init)*onset.thresh+init
      time.on <- x[min(which(y.pred > onset))]
      time.off <- x[max(which(y.pred > onset))]
    } else {
      # Double, concave
      end <- mean(tail(y.pred, 5))
      min <- sort(y.pred, decreasing=F)[5]
      onset <- (end-min)*onset.thresh+min
      inflection <- mean(x[order(y.pred)[1:5]])
      time.on <- x[min(intersect(which(y.pred > onset), which(x > inflection)))]
      time.off <- x[max(intersect(which(y.pred > onset), which(x > inflection)))]
    }
  }  
  
  ## Format and return the results
  if (best.fit == "linear") {
    return(c(type=0, time.on=time.on, time.off=time.off, linear))
  } else if (best.fit == "single") {
    return(c(type=1, time.on=time.on, time.off=time.off, single))
  } else if (best.fit == "double") {
    return(c(time.on=time.on, time.off=time.off, double))
  } else {
    return(NA)
  }
}