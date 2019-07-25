

#' Derivative of single sigmoid
#' 
#' @param x (Numeric)
#' @param b1 (Numeric) Rate parameter of sigmoid
#' @param h0 (Numeric) Initial level
#' @param h1 (Numeric) Final level
#' @param t1 (Numeric) Midpoint of sigmoid function
#' 
#' @return Returns numeric vector of first derivative of single sigmoid with given parameters at each x.
#' 
#' @keywords internal
impulse.deriv.single <- function(x, b1, h0, h1, t1){
  f1 <- impulse.single(x, b1, h0, h1, t1)
  d_f1 <- b1*(f1-h0)*(f1-h1)/(h1-h0)
  return(d_f1)
}

#' Second derivative of single sigmoid
#' 
#' @param x (Numeric)
#' @param b1 (Numeric) Rate parameter of sigmoid
#' @param h0 (Numeric) Initial level
#' @param h1 (Numeric) Final level
#' @param t1 (Numeric) Midpoint of sigmoid function
#' 
#' @return Returns numeric vector of second derivative of single sigmoid with given parameters at each x.
#' 
#' @keywords internal
impulse.dderiv.single <- function(x, b1, h0, h1, t1){
  f1 <- impulse.single(x, b1, h0, h1, t1)
  dd_f1 <- b1^2*(f1-h0)*(f1-h1)*(2*f1-h1-h0)/(h0-h1)^2
  return(dd_f1)
}

#' Single derivative of double sigmoid ("impulse")
#' 
#' @param x (Numeric)
#' @param b1 (Numeric) Rate parameter of first sigmoid (should be negative for impulse)
#' @param b2 (Numeric) Rate parameter of second sigmoid (should be negative for impulse)
#' @param h0 (Numeric) Initial level
#' @param h1 (Numeric) Peak level
#' @param h2 (Numeric) Final level
#' @param t1 (Numeric) Midpoint of first sigmoid function
#' @param t2 (Numeric) Midpoint of second sigmoid function
#' 
#' @return Returns numeric vector of first derivative of double sigmoid ("impulse") with given parameters at each x.
#' 
#' @keywords internal
impulse.deriv.double <- function(x, b1, b2, h0, h1, h2, t1, t2){
  f1 <- impulse.single(x, b1, h0, h1, t1)
  f2 <- impulse.single(x, b2, h2, h1, t2)
  dd_sg1 <- (b2*f1*(f2-h1)*(h1-h0)*(f2-h2)+b1*f2*(f1-h0)*(h1-h2)*(f1-h1))/(h1*(h1-h0)*(h1-h2))
  return(dd_sg1)
}

#' Second derivative of double sigmoid ("impulse")
#' 
#' @param x (Numeric)
#' @param b1 (Numeric) Rate parameter of first sigmoid (should be negative for impulse)
#' @param b2 (Numeric) Rate parameter of second sigmoid (should be negative for impulse)
#' @param h0 (Numeric) Initial level
#' @param h1 (Numeric) Peak level
#' @param h2 (Numeric) Final level
#' @param t1 (Numeric) Midpoint of first sigmoid function
#' @param t2 (Numeric) Midpoint of second sigmoid function
#' 
#' @return Returns numeric vector of second derivative of double sigmoid ("impulse") with given parameters at each x.
#' 
#' @keywords internal
impulse.dderiv.double <- function(x, b1, b2, h0, h1, h2, t1, t2){
  f1 <- impulse.single(x, b1, h0, h1, t1)
  f2 <- impulse.single(x, b2, h2, h1, t2)
  dd_sg1 <- (b2^2*f1*(f2-h1)*(h0-h1)^2*(f2-h2)*(2*f2-h1-h2) + 2*b1*b2*(f1-h0)*(f1-h1)*(h0-h1)*(h1-f2)*(f2-h2)*(h1-h2) + b1^2*f2*(f1-h0)*(f1-h1)*(2*f1-h0-h1)*(h1-h2)^2)/((h0 - h1)^2*h1*(h1 - h2)^2)
  return(dd_sg1)
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
#' 
#' @keywords internal
impulse.single <- function(x, b1, h0, h1, t1) {
  h0 + (h1 - h0)/(1 + exp(b1*(x - t1)))
}

#' Evaluate the impulse model (double sigmoid) function, given parameters
#'
#' @param x (Numeric)
#' @param b1 (Numeric) Rate parameter of first sigmoid (should be negative for impulse)
#' @param b2 (Numeric) Rate parameter of second sigmoid (should be negative for impulse)
#' @param h0 (Numeric) Initial level
#' @param h1 (Numeric) Peak level
#' @param h2 (Numeric) Final level
#' @param t1 (Numeric) Midpoint of first sigmoid function
#' @param t2 (Numeric) Midpoint of second sigmoid function
#' 
#' @return Double-sigmoid impulse function with parameters b1, b2, h0, h1, h2, t1, t2, as evaluated for the vector x.
#' 
#' @keywords internal
impulse.double <- function(x, b1, b2, h0, h1, h2, t1, t2) {
  f1 <- h0 + (h1 - h0)/(1 + exp(b1*(x - t1)))
  f2 <- h2 + (h1 - h2)/(1 + exp(b2*(x - t2)))
  sg1 <- f1*f2/h1
  return(sg1)
}

#' Starting Conditions for Impulse Model (double sigmoid)
#' 
#' This function determines \code{k} sets of initial starting conditions for impulse fitting by choosing combinations of points from the data randomly. h0 and h2 (initial and final values) are taken from points chosen at random from data at the beginning and end of the data series (along the x-axis). h1 is chosen from points in the middle of the data (along the x-axis) that are either greater than both h0 and h2 or less than both h0 and h2. t1 and t2 are then estimated as the midpoints (along x-axis) between points randomly chosen to determine h0, h1, and h2. b1 and b2 are estimated as linear slopes between points chosen to determine h0, h1, and h2. If \code{k} is larger than the number of possible combinations that fit these conditions, additional sets of starting conditions are chosen totally randomly that fit within the bounds of the existing data with a small margin of padding.
#' 
#' @param x (Numeric vector) Expression data to fit (x)
#' @param y (Numeric vector) Expression data to fit (y)
#' @param k (Numeric) Number of starting conditions to determine
#' @param h0.frac (Numeric, 0 to 1) Portion of data (along x-axis, from beginning) to use to estimate h0
#' @param h2.frac (Numeric, 0 to 1) Portion of data (along x-axis, from end) to use to estimate h2
#' @param interpolate (Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points. Default (\code{NULL}) does not do any interpolation.
#' 
#' @return Returns a data frame of potential starting conditions for fitting functions for double sigmoid ("impulse") model
#' 
#' @keywords internal
impulse.start.double <- function(x, y, k, h0.frac=0.2, h2.frac=0.2, interpolate=NULL) {
  if(!is.null(interpolate)){
    n_xy <- approx(x, y, n=interpolate)
    x <- n_xy$x
    y <- n_xy$y
  }
  x_range <- range(x)
  ind.first.tenth <- which(x <= (x_range[1] + h0.frac*diff(x_range)))
  ind.last.tenth <- which(x >= (x_range[2] - h2.frac*diff(x_range)))
  mean.first.tenth <- mean(y[ind.first.tenth])
  mean.last.tenth <- mean(y[ind.last.tenth])
  mid.range <- sort(c(mean.first.tenth, mean.last.tenth), decreasing=F)
  mid.ind <- setdiff(c(1:length(x)),c(ind.first.tenth,ind.last.tenth))
  mid.points <- mid.ind[union(which(y[mid.ind]>(mid.range[2])), which(y[mid.ind]<(mid.range[1])))]
  ind.tbl <- c()
  num_pair <- 0
  ind_pair <- rbind(rep(ind.first.tenth,length(ind.last.tenth)),rep(ind.last.tenth,each=length(ind.first.tenth)))
  if(length(mid.points)>0){
    ind_trp <- rbind(ind_pair[,rep(1:dim(ind_pair)[2],length(mid.points))],rep(mid.points,each=dim(ind_pair)[2]))
    y_trp <- rbind(y[ind_trp[1,]],y[ind_trp[2,]],y[ind_trp[3,]])
    order_trp <- apply(y_trp, 2, order)
    rm_ind <- which(order_trp[2,]==3)
    if(length(rm_ind)>0){
      ind_trp <- ind_trp[,-rm_ind]
    }
    if(length(ind_trp[1,])>0){
      h0s <- y[ind_trp[1,]]
      h1s <- y[ind_trp[3,]]
      h2s <- y[ind_trp[2,]]
      t1s <- 0.5*(x[ind_trp[1,]]+x[ind_trp[3,]])
      t2s <- 0.5*(x[ind_trp[2,]]+x[ind_trp[3,]])
      b1s <- -(h0s-h1s)/(x[ind_trp[1,]]-x[ind_trp[3,]])
      b2s <- -(h2s-h1s)/(x[ind_trp[2,]]-x[ind_trp[3,]])
      starts <- data.frame(
        h0=h0s,
        h1=h1s,
        h2=h2s,
        t1=t1s,
        t2=t2s,
        b1=b1s,
        b2=b2s,
        stringsAsFactors=F
      )
      starts <- starts[,c("b1","b2","h0","h1","h2","t1","t2")]
      num_pair <- dim(starts)[1] 
      if(num_pair>=k){
        starts <- starts[sample(1:dim(starts)[1],k),]
        return(starts)
      }
    }
  }
  #print(paste0("generating ", k-num_pair, " random initial conditions"))
  #### if couldn't guess enough starting points, generate some random ones
  y_range <- range(y)
  y_range.use <- c(max(0,y_range[1]-0.2*diff(y_range)),y_range[2]+0.2*diff(y_range))
  x_range.use <- c(max(0,x_range[1]-0.2*diff(x_range)),x_range[2]+0.2*diff(x_range))
  num_points <- k-num_pair
  y_trp <- rbind(runif(num_points, y_range.use[1], y_range.use[2]),runif(num_points, y_range.use[1], y_range.use[2]),runif(num_points, y_range.use[1], y_range.use[2]))
  x_trp <- rbind(runif(num_points, x_range.use[1], x_range.use[2]),runif(num_points, x_range.use[1], x_range.use[2]),runif(num_points, x_range.use[1], x_range.use[2]))
  x_trp <- apply(x_trp,2,sort)
  h0s <- y_trp[1,]
  h1s <- y_trp[2,]
  h2s <- y_trp[3,]
  t1s <- 0.5*(x_trp[1,]+x_trp[2,])
  t2s <- 0.5*(x_trp[2,]+x_trp[3,])
  b1s <- -(h0s-h1s)/(x_trp[1,]-x_trp[2,])
  b2s <- -(h2s-h1s)/(x_trp[3,]-x_trp[2,])

  starts.2 <- data.frame(
    h0=h0s,
    h1=h1s,
    h2=h2s,
    t1=t1s,
    t2=t2s,
    b1=b1s,
    b2=b2s,
    stringsAsFactors=F
  )
  starts.2 <- starts.2[,c("b1","b2","h0","h1","h2","t1","t2")]
  if(num_pair==0){
    starts.2$b1 <- -abs(starts.2$b1)
    starts.2$b2 <- -abs(starts.2$b2)
    return(starts.2)
  }else{
    starts.2$b1 <- -abs(starts.2$b1)
    starts.2$b2 <- -abs(starts.2$b2)
    starts$b1 <- -abs(starts$b1)
    starts$b2 <- -abs(starts$b2)
    return(rbind(starts,starts.2))
  }
}

#' Starting conditions for single sigmoid model
#' 
#' This function determines \code{k} sets of initial starting conditions for single sigmoid fitting by choosing combinations of points from the data randomly. h0 and h1 (initial and final values) are taken from points chosen at random from data with the highest and lowest y-axis values. t1 and b1 are then estimated from two randomly chosen points in the middle of the data (along the y-axis). If \code{k} is larger than the number of possible combinations that fit these conditions, additional sets of starting conditions are chosen totally randomly that fit within the bounds of the existing data with a small margin of padding.
#' 
#' @param x (Numeric vector) Expression data to fit (x)
#' @param y (Numeric vector) Expression data to fit (y)
#' @param k (Numeric) Number of starting conditions to determine
#' @param h.frac (Numeric, 0 to 1) Portion of data (highest along y-axis) to use to estimate h0/h1
#' @param l.frac (Numeric, 0 to 1) Portion of data (lowest along y-axis) to use to estimate h0/h1
#' @param interpolate (Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points. Default (\code{NULL}) does not do any interpolation.
#' 
#' @return Returns a data frame of potential starting conditions for fitting functions for double sigmoid ("impulse") model
#' 
#' @keywords internal
impulse.start.single <- function(x, y, k, h.frac=0.2, l.frac=0.2, interpolate=NULL) {
  if(!is.null(interpolate)){
    n_xy <- approx(x, y, n=k)
    x <- n_xy$x
    y <- n_xy$y
  }
  num_st <- 0
  y_range <- c(min(y),max(y))
  starts <- data.frame()

  mean.lowest.data <- y_range[1]+diff(y_range)*l.frac
  mean.highest.data <- y_range[2]-diff(y_range)*h.frac
  high.ind <- which(y>=mean.highest.data)
  low.ind <- which(y<=mean.lowest.data)
  ind.pairs <- rbind(rep(high.ind,length(low.ind)),rep(low.ind,each=length(high.ind)))
  ind.pairs <- apply(ind.pairs,2,sort)
  too.close <- which((ind.pairs[2,]-ind.pairs[1,])<3)
  if(length(too.close)>0){
    ind.pairs <- ind.pairs[,-too.close]
  }
  if(dim(ind.pairs)[2]>0){
    ind.tbl <- c()
    for(i in 1:dim(ind.pairs)[2]){
      inds.mid <- ind.pairs[1,i]:ind.pairs[2,i]
      inds.l <- inds.mid[which(x[inds.mid]<0.5*(x[ind.pairs[1,i]]+x[ind.pairs[2,i]]))]
      inds.r <- inds.mid[which(x[inds.mid]>0.5*(x[ind.pairs[1,i]]+x[ind.pairs[2,i]]))]
      mid_pair <- rbind(rep(inds.l, length(inds.r)), rep(inds.r, each=length(inds.l)))
      ind.tbl <- cbind(ind.tbl, rbind(rep(ind.pairs[1,i],dim(mid_pair)[2]),rep(ind.pairs[2,i],dim(mid_pair)[2]), mid_pair))
    }
    if(!is.null(ind.tbl)){
      num_st <- dim(ind.tbl)[2]
      starts <- data.frame(
        h0=y[ind.tbl[1,]],
        h1=y[ind.tbl[2,]],
        t1=(x[ind.tbl[3,]]+x[ind.tbl[4,]])/2,
        b1=-(y[ind.tbl[4,]]-y[ind.tbl[3,]])/(x[ind.tbl[4,]]-x[ind.tbl[3,]]),
        slope.low=rep(-Inf,dim(ind.tbl)[2]),
        slope.high=rep(0,dim(ind.tbl)[2]),
        stringsAsFactors=F
      )
    }
  }

  #### if couldn't guess enough starting points, generate some random ones
  if(num_st<k){
    y_range.use <- c(max(0,y_range[1]-0.2*diff(y_range)),y_range[2]+0.2*diff(y_range))
    x_range <- range(x)
    x_range.use <- c(max(0,x_range[1]-0.2*diff(x_range)),x_range[2]+0.2*diff(x_range))
    num_points <- k-num_st
    x_mid <- 0.5*(x_range[1]+x_range[2])
    ymean_l <- mean(y[which(x<x_mid)])
    ymean_r <- mean(y[which(x>x_mid)])

    y_tbl <- matrix(runif(num_points*4,y_range.use[1], y_range.use[2]),nrow=4)
    x_tbl <- matrix(runif(num_points*4,x_range.use[1], x_range.use[2]),nrow=4)
    y_tbl <- apply(y_tbl, 2, sort)
    x_tbl <- apply(x_tbl, 2, sort)
    if(ymean_l>ymean_r){
      y_tbl <- y_tbl[c(4,1,3,2),,drop=0]
    }else{
      y_tbl <- y_tbl[c(1,4,2,3),,drop=0]
    }
    x_tbl <- x_tbl[c(1,4,2,3),,drop=0]
    starts.2 <- data.frame(
      h0=y_tbl[1,],
      h1=y_tbl[2,],
      t1=(x_tbl[3,]+x_tbl[4,])/2,
      b1=-(y_tbl[4,]-y_tbl[3,])/(x_tbl[4,]-x_tbl[3,]),
      slope.low=rep(-Inf,num_points),
      slope.high=rep(0,num_points),
      stringsAsFactors=F
      )
    starts <- rbind(starts,starts.2)
    starts$b1 <- -abs(starts$b1)
    return(starts[,c("b1","h0","h1","t1","slope.low","slope.high")])
  }else{
    starts$b1 <- -abs(starts$b1)
    return(starts[,c("b1","h0","h1","t1","slope.low","slope.high")])
  }
}


#' Fit Impulse Model (Double sigmoid)
#' 
#' Attempts to fit a double sigmoid ("impulse" model) to given expression data using expectation maximization.
#' 
#' @importFrom minpack.lm nlsLM nls.lm.control
#' 
#' @param x (Numeric) Expression data
#' @param y (Numeric) Expression data
#' @param k (Numeric) Number of starting conditions to try
#' @param maxiter (Numeric) Maximum number of iterations to try fitting
#' @param interpolate (Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points for choosing potential starting conditions. Default (\code{NULL}) does not do any interpolation. Interpolation is NOT used during actual function fitting.
#' 
#' @return (Named vector, length 8) Returns the coefficients and sum of squared residuals of the best fitting double sigmoid out of \code{k} conditions. If all attempts to fit fail to converge, returns vector of same length, with all values \code{NA}.
#' 
#' @keywords internal
impulse.fit.double <- function(x, y, k=20, maxiter=200, interpolate=NULL) {
  y_range <- c(min(y),max(y))
  low_bd <- y_range[1]-diff(y_range)
  high_bd <- y_range[2]+diff(y_range)
  df <- data.frame(x = x, y = y)
  # Get starting conditions
  starts <- impulse.start.double(x, y, k, interpolate=interpolate)
  
  # Try the fits!
  con <- nls.lm.control(maxiter=maxiter)
  fits <- lapply(1:k, function(f) try(nlsLM(y ~ impulse.double(x, b1, b2, h0, h1, h2, t1, t2), data=df, start=starts[f,], lower=c(-Inf, -Inf, low_bd, low_bd, low_bd, min(x)-1, min(x)-1), upper=c(Inf, Inf, high_bd, high_bd, high_bd, max(x)+1, max(x)+1), control=con), silent=T))
  #print("fit.double.ed")

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

#' Fit Single Sigmoid Model
#' 
#' Attempts to fit a single sigmoid to given expression data using expectation maximization.
#' 
#' @importFrom minpack.lm nlsLM nls.lm.control
#' 
#' @param x (Numeric) Expression data
#' @param y (Numeric) Expression data
#' @param k (Numeric) Number of starting conditions to try
#' @param maxiter (Numeric) Maximum number of iterations to try fitting
#' @param interpolate (Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points for choosing potential starting conditions. Default (\code{NULL}) does not do any interpolation. Interpolation is NOT used during actual function fitting.
#' 
#' @return (Named vector, length 5) Returns the coefficients and sum of squared residuals of the best fitting double sigmoid out of \code{k} conditions. If all attempts to fit fail to converge, returns vector of same length, with all values \code{NA}.
#' 
#' @keywords internal
impulse.fit.single <- function(x, y, k=20, maxiter=200, interpolate=NULL) {
  y_range <- c(min(y),max(y))
  low_bd <- y_range[1]-diff(y_range)*0.5
  high_bd <- y_range[2]+diff(y_range)*0.5
  df <- data.frame(x = x, y = y)
  # Get starting conditions
  #print("start.single.st")
  starts <- impulse.start.single(x, y, k, interpolate=interpolate)
  #print(starts)
  
  # Try the fits!
  con <- nls.lm.control(maxiter=maxiter)
  #print("fit.single.st")
  fits <- lapply(1:nrow(starts), function(f) try(nlsLM(y ~ impulse.single(x, b1, h0, h1, t1), data=df, start=starts[f,c("b1","h0","h1","t1")], control=con, lower=c(starts[f,"slope.low"], low_bd, low_bd, min(x)-1), upper=c(starts[f,"slope.high"], high_bd, high_bd, max(x)+1)), silent=T))
  # Figure out the error of the fits and choose the best one.
  #print("fit.single.ed")
  best.fit <- which.min(unlist(lapply(fits, function(f) {
    if (class(f) == "try-error") return(NA)
    sum(resid(f)^2)
  })))
  # Return values
  if (length(best.fit)==0) return(rep(NA, 5))
  to.ret <- c(coef(fits[[best.fit]]), sum(resid(fits[[best.fit]])^2))
  names(to.ret)[5] <- "sum.resid2"
  #print(to.ret)
  return(to.ret)
}



#' Fit gene expression data with an impulse model
#' 
#' This fits either a linear, single sigmoid, or double sigmoid ("impulse") model to given x and y coordinates (e.g. expression data). This function is primarily used by \code{\link{geneCascadeProcess}}.
#' 
#' Thanks to Yiqun Wang for considerable improvements to the impulse fitting functions.
#' 
#' @references Chechik G, Oh E, Rando O, Weissman J, Regev A, Koller D. "Activity motifs reveal principles of timing in transcriptional control of the yeast metabolic network." Nat Biotechnol. 2008 Nov;26(11):1251-9. doi: 10.1038/nbt.1499
#' 
#' @param x (Numeric) Expression data to fit (e.g. pseudotime)
#' @param y (Numeric) Expression data to fit (e.g. expression level)
#' @param k (Numeric) Number of sets of starting conditions to try
#' @param interpolate (Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points for choosing potential starting conditions. Default (\code{NULL}) does not do any interpolation. Interpolation is NOT used during actual function fitting.
#' @param pulse.only (Logical) If \code{TRUE}, filters out double sigmoid functions that are monotonous and prefers either single sigmoid or linear fit for them instead.
#' @param min.slope (Numeric) For \code{pulse.only} filtering, first derivative of double sigmoid function must include positive and negative values with absolute value > \code{min.slope}.
#' 
#' @return (Named List, length 4) Parameters & information about best fitting model:
#' \itemize{ 
#'   \item{\strong{\code{type}}: 0/1/2 = linear/single sigmoid/double sigmoid}
#'   \item{\strong{\code{time.on}}: Gene onset pseudotimes}
#'   \item{\strong{\code{time.off}}: Gene offset pseudotimes}
#'   \item{\strong{\code{model}}: Parameters of the best model}
#' }
#' 
#' @export

impulseFit <- function(x, y, k=50, interpolate=NULL, pulse.only=T, min.slope=0.1) {
  ## Fit linear regression, single sigmoid, and double sigmoids
  linear <- lm(y~x, data=list(x=x, y=y))
  linear <- c(coef(linear), sum(resid(linear)^2))
  names(linear) <- c("b", "a", "sum.resid2")
  single <- impulse.fit.single(x=x, y=y, k=k, interpolate=interpolate)
  double <- impulse.fit.double(x=x, y=y, k=k, interpolate=interpolate)

  num.x <- max(c(200,length(x)))
  ## Choose the best fit according to residual
  resid2 <- c(linear["sum.resid2"],single["sum.resid2"],double["sum.resid2"])
  names(resid2) <- c("linear","single","double")
  best.fit <- names(sort(resid2))[1]
  ## modified by YW @2018-09-07 to reject double sigmoid fit if curve is monotonous
  if (pulse.only && best.fit =='double'){
    if(length(sort(resid2))>1){
      x_ <- seq(min(x),max(x),length.out = num.x)
      dy <- impulse.deriv.double(x_, b1=double['b1'], b2=double['b2'], h0=double['h0'], h1=double['h1'], h2=double['h2'], t1=double['t1'], t2=double['t2'])
      dy_ <- dy
      if(any(abs(dy) < min.slope)){
        dy_[which(abs(dy) < min.slope)] <- 0
      }
      if(all(dy_>=0) || all(dy_<=0)){
        best.fit <- names(sort(resid2))[2]
      }else{
        dy_sign_switch <- which(diff(dy>=0)!=0)
        y_ <- impulse.double(x_, b1=double['b1'], b2=double['b2'], h0=double['h0'], h1=double['h1'], h2=double['h2'], t1=double['t1'], t2=double['t2'])
        y_range <- max(y_)-min(y_)
        if(length(dy_sign_switch)>1){
          y_diff <- diff(y_[dy_sign_switch])
          if(any(abs(y_diff)<0.25*y_range)){
            best.fit <- names(sort(resid2))[2]
          }
        }else if(length(dy_sign_switch)==1){
          if(any(abs(diff(c(y_[1],y_[dy_sign_switch],y_[length(y_)])))<0.1*y_range)){
            best.fit <- names(sort(resid2))[2]
          }
        }
      }
    }
  }
  ## Figure out time of onset and offset
  x_ <- seq(min(x),max(x),length.out = num.x)
  if (best.fit%in%c('linear','single','double')) {
    i <- get(best.fit)
    if (best.fit=="linear") {
      if (i["a"] > 0) {
        # Positive slope
        time.off <- Inf
        time.on <- as.numeric(-1*i["b"]/i["a"])
      } else {
        # Negative slope
        time.on <- -Inf
        time.off <- as.numeric(-1*i["b"]/i["a"])
      }
    }else if (best.fit == "single") {
      y_ <- impulse.single(x_, b1=i['b1'], h0=i['h0'], h1=i['h1'], t1=i['t1'])
      dy_max <- impulse.deriv.single(i['t1'], b1=i['b1'], h0=i['h0'], h1=i['h1'], t1=i['t1'])
      dy <- impulse.deriv.single(x_, b1=i['b1'], h0=i['h0'], h1=i['h1'], t1=i['t1'])
      y_hlf_max <- c(((2 + sqrt(2))*i['h0'] - (-2 + sqrt(2))*i['h1'])/4, (-(-2 + sqrt(2))*i['h0'] + (2 + sqrt(2))*i['h1'])/4)
      if(dy_max >0){
        time.on <- x_[which.min(abs(y_-min(y_hlf_max)))]
        time.off <- Inf
      }else{
        time.off <- x_[which.min(abs(y_-max(y_hlf_max)))]
        time.on <- -Inf
      }
    } else {
      time.on <- c()
      time.off <- c()
      y_ <- impulse.double(x_, b1=i['b1'], b2=i['b2'], h0=i['h0'], h1=i['h1'], h2=i['h2'], t1=i['t1'], t2=i['t2'])
      y_range <- max(y_)-min(y_)
      ddy <- impulse.dderiv.double(x_, b1=i['b1'], b2=i['b2'], h0=i['h0'], h1=i['h1'], h2=i['h2'], t1=i['t1'], t2=i['t2'])
      ddy_sign_switch <- which(diff(ddy>=0)!=0)
      if(length(ddy_sign_switch)==0){
        print("no local maximum/minimum found for the first derivative")
      }
      dy <- impulse.deriv.double(x_, b1=i['b1'], b2=i['b2'], h0=i['h0'], h1=i['h1'], h2=i['h2'], t1=i['t1'], t2=i['t2'])
      dy_sign_switch <- which(diff(dy>=0)!=0)
      if(length(dy_sign_switch)==0){
        print("no peak found for double sigmoid fit")
        ## this should not happen if pulse.only=T
        if(all(dy>=0) || all(dy<=0)){
          max_ind <- which.max(abs(dy))
          hlf_max <- which.min(abs(dy[1:max_ind]-0.5*dy[max_ind]))
          if(dy[max_ind]>0){
            time.on <- x_[hlf_max]
            time.off <- Inf
          }else{
            time.off <- x_[hlf_max]
            time.on <- -Inf
          }
        }
      }else{
        if(dy_sign_switch[1]!=1){
          seg_ind <- 1:dy_sign_switch[1]
          dy_seg <- dy[seg_ind]
          max_ind <- which.max(abs(dy_seg))
          if(abs(y_[1]-y_[dy_sign_switch[1]])>0.25*y_range){
            hlf_max <- which.min(abs(dy_seg[1:max_ind]-0.5*dy_seg[max_ind]))
            hlf_max <- seg_ind[hlf_max]
            if(dy_seg[max_ind] >0){
              time.on <- x_[hlf_max]
            }else{
              time.off <- x_[hlf_max]
            }
          }
        }
        if(length(dy_sign_switch)>1){
          for(j in 2:length(dy_sign_switch)){
            seg_ind <- (dy_sign_switch[j-1]+1):dy_sign_switch[j]
            if(length(seg_ind)>2){
              dy_seg <- dy[seg_ind]
              max_ind <- which.max(abs(dy_seg))
              #h0 <- y_[dy_sign_switch[j-1]]
              #h1 <- y_[dy_sign_switch[j]]
              #y_hlf_max <- c(((2 + sqrt(2))*h0 - (-2 + sqrt(2))*h1)/4, (-(-2 + sqrt(2))*h0 + (2 + sqrt(2))*h1)/4)
              hlf_max <- which.min(abs(dy_seg[1:max_ind]-0.5*dy_seg[max_ind]))
              hlf_max <- seg_ind[hlf_max]
              if(dy_seg[max_ind] >0){
                #time.on <- c(time.on, x_[seg_ind[which.min(abs(y_[seg_ind]-min(y_hlf_max)))]])
                time.on <- c(time.on, x_[hlf_max])
              }else{
                #time.off <- c(time.off, x_[seg_ind[which.min(abs(y_[seg_ind]-max(y_hlf_max)))]])
                time.off <- c(time.off, x_[hlf_max])
              }
            }
          }
        }
        #print(dy_sign_switch)
        if(dy_sign_switch[length(dy_sign_switch)]<(length(x_)-2)) {
          seg_ind <- (dy_sign_switch[length(dy_sign_switch)]+1):length(x_)
          dy_seg <- dy[seg_ind]
          max_ind <- which.max(abs(dy_seg))
          if((max(y_[seg_ind])-min(y_[seg_ind]))>0.25*y_range){
            # if(abs(dy_seg[length(dy_seg)])<0.1){
            #   h0 <- y_[seg_ind[1]]
            #   h1 <- y_[length(y_)]
            #   y_hlf_max <- c(((2 + sqrt(2))*h0 - (-2 + sqrt(2))*h1)/4, (-(-2 + sqrt(2))*h0 + (2 + sqrt(2))*h1)/4)
            #  if(dy_seg[max_ind] >0){
            #     time.on <- c(time.on, x_[seg_ind[which.min(abs(y_[seg_ind]-min(y_hlf_max)))]])
            #   }else{
            #     time.off <- c(time.off, x_[seg_ind[which.min(abs(y_[seg_ind]-max(y_hlf_max)))]])
            #   }
            # }else{
            hlf_max <- which.min(abs(dy_seg[1:max_ind]-0.5*dy_seg[max_ind]))
            hlf_max <- seg_ind[hlf_max]
            if(dy_seg[max_ind] >0){
              time.on <- c(time.on, x_[hlf_max])
            }else{
              time.off <- c(time.off, x_[hlf_max])
            }
            # }
          }
        }
      }
      if(is.null(time.on)){
        time.on <- -Inf
      }
      if(is.null(time.off)){
        time.off <- Inf
      }
    }
  }
  
  ## Format and return the results
  
  if (best.fit == "linear") {
    #return(c(type=0, time.on=time.on, time.off=time.off, linear))
    return(list(type=0, time.on=time.on, time.off=time.off, model=linear))
  } else if (best.fit == "single") {
    #return(c(type=1, time.on=time.on, time.off=time.off, single))
    return(list(type=1, time.on=time.on, time.off=time.off, model=single))
  } else if (best.fit == "double") {
    return(list(type=2, time.on=time.on, time.off=time.off, model=double))
  } else {
    return(NA)
  }
}