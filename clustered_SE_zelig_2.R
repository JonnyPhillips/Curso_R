library(rms)
library(tidyverse)
library(zeligverse)
library(nycflights13)

model1 <- flights %>% zelig(dep_delay~dep_time + origin,data=.,model="ls")
summary(model1)

model2 <- flights %>% zelig(dep_delay~dep_time + origin,data=.,model="ls")
summary(model1)

model3 <- flights %>% ols(dep_delay~dep_time+origin, data=.,x=T,y=T)
summary(model3)

fit <- model1
cluster <- flights$origin
method <- "huber"

model3$residuals
model1$get_residuals()
names(model1)
model1$
residuals(model3, type="score")

robcov_zelig <- function(fit,cluster){
  method <- match.arg(method)
  var <- vcov(fit)
  vname <- dimnames(var)[[1]]
  if (fit$name=="ls")
    var <- fit$get_df_residual * var/sum(fit$get_residuals()^2)
  X <- as.matrix(fit$get_residuals())
  X <- as.matrix(residuals(model3, type = if (method == "huber") 
    "score"
    else "hscore"))
  n <- nrow(X)
  if (missing(cluster)) {
    clusterInfo <- NULL
    cluster <- 1:n
  }
  else {
    if (any(is.na(cluster))) 
      stop("cluster contains NAs")
    clusterInfo <- list(name = deparse(substitute(cluster)))
  }
  if (length(cluster) != n) 
    stop("length of cluster does not match number of observations used in fit")
  cluster <- as.factor(cluster)
  p <- ncol(var)
  j <- is.na(X %*% rep(1, ncol(X)))
  if (any(j)) {
    X <- X[!j, , drop = FALSE]
    cluster <- cluster[!j, drop = TRUE]
    n <- length(cluster)
  }
  j <- order(cluster)
  X <- X[j, , drop = FALSE]
  clus.size <- table(cluster)
  if (length(clusterInfo)) 
    clusterInfo$n <- length(clus.size)
  clus.start <- c(1, 1 + cumsum(clus.size))
  nc <- length(levels(cluster))
  clus.start <- clus.start[-(nc + 1)]
  storage.mode(clus.start) <- "integer"
  W <- matrix(.Fortran(F_robcovf, n, p, nc, clus.start, clus.size, 
                       X, double(p), double(p * p), w = double(p * p))$w, nrow = p)
  adjvar <- var %*% W %*% var
  fit$orig.var <- var
  fit$var <- adjvar
  fit$clusterInfo <- clusterInfo
  fit
}


###
flights2 <- flights %>% filter(!is.na(dep_delay) & !is.na(dep_time) & !is.na(origin))

model1 <- flights2 %>% zelig(dep_delay~dep_time + origin,data=.,model="ls")
summary(model1)

model_lm <- flights2 %>% 
  lm(dep_delay~dep_time + origin,data=.)

cl   <- function(dat,fm, cluster){
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
  coeftest(fm, vcovCL) }

summary(model_lm)
cl(flights2,model_lm,flights2$origin)

model <- model1
cluster <- flights2$origin

estfun_zelig <- function (x, ...) {
  xmat <- model.matrix(x$formula,data=x$data)
  xmat <- naresid(NULL, xmat) #Changed, probably not working
  if (any(alias <- is.na(coef(x)))) 
    xmat <- xmat[, !alias, drop = FALSE]
  wts <- weights(x)
  if (is.null(wts)) 
    wts <- 1
  res <- x$get_residuals()[[1]]
  rval <- as.vector(res) * wts * xmat
  attr(rval, "assign") <- NULL
  attr(rval, "contrasts") <- NULL
  if (is.zoo(res)) 
    rval <- zoo(rval, index(res), attr(res, "frequency"))
  if (is.ts(res)) 
    rval <- ts(rval, start = start(res), frequency = frequency(res))
  return(rval)
}
x <- model_lm
meat=crossprod(uj)/N

sandwich_zelig <- function(model,meat){
    estd <- estfun_zelig(model)
    sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/sum(dim(estd))) #Double check, but about right
    bread <- model$get_vcov()[[1]]*sum(dim(estd))*sigma^2 #dim not exact - 4 rows duplicated for some reason...
    n <- nrow(estd)
    return(1/n * (bread %*% meat %*% bread..))
  }
names(model1)


temp <- summary(model1)

sx$cov.unscaled * as.vector(sum(sx$df[1:2])) *sx$sigma^2 #Yes! Is identical to with zelig - so just need to find sigme (residual SE)
methods(bread)
getAnywhere(bread.lm)

function (x, ...) 
{
  if (!is.null(x$na.action)) 
    class(x$na.action) <- "omit"
  sx <- summary.lm(x)
  return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))
}
names(sx)
cl_zelig   <- function(model, cluster){
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- dim(vcov(model1)[[1]])[1] #Probably wrong!
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun_zelig(model),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich(model, meat=crossprod(uj)/N)
  coeftest(model, vcovCL) }

cl_zelig(model1,flights2$origin)


