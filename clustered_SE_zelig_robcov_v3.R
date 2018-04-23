library(rms)
library(tidyverse)
library(zeligverse)
library(nycflights13)

###
flights2 <- flights %>% filter(!is.na(dep_delay) & !is.na(dep_time) & !is.na(origin))

model1 <- flights2 %>% zelig(dep_delay~dep_time + origin,data=.,model="ls")
summary(model1)

model_ols <- flights2 %>% ols(dep_delay~dep_time + origin,data=.,x=T,y=T)

model1_logit <- flights2 %>% mutate(Late=ifelse(dep_delay>0,1,0)) %>% 
  zelig(Late~dep_time + origin,data=.,model="logit")
summary(model1_logit)

model_lrm <- flights2 %>% mutate(Late=ifelse(dep_delay>0,1,0)) %>% 
  lrm(Late~dep_time + origin,data=.,x=T,y=T)

model_glm <- flights2 %>% mutate(Late=ifelse(dep_delay>0,1,0)) %>% 
  glm(Late~dep_time + origin,data=.,x=T,y=T,family="binomial")

model_ols %>% robcov(flights2$origin)

fit <- model_ols
cluster <- flights2$origin
method <- "huber"

model <- model1_logit

residuals(model_ols, "working")

residuals(model_glm)
model1_logit$get_residuals()
residuals(model1_logit)
residuals.glm(model_glm,"working")

names(model1)
model1$data[[]]

y <- object$y
r <- object$residuals
mu <- object$fitted.values
wts <- object$prior.weights
switch(type, deviance = , pearson = , response = if (is.null(y)) {
  mu.eta <- object$family$mu.eta
  eta <- object$linear.predictors
  y <- mu + r * mu.eta(eta)
})
res <- switch(type, deviance = if (object$df.residual > 0) {
  d.res <- sqrt(pmax((object$family$dev.resids)(y, mu, 
                                                wts), 0))
  ifelse(y > mu, d.res, -d.res)
} else rep.int(0, length(mu)), pearson = (y - mu) * sqrt(wts)/sqrt(object$family$variance(mu)), 
working = r, response = y - mu, partial = r)
if (!is.null(object$na.action)) 
  res <- naresid(object$na.action, res)
if (type == "partial") 
  res <- res + predict(object, type = "terms")
res


estfun_zelig <- function (model, ...) {
  if (model$fn=="stats::lm") {
    xmat <- model.matrix(model$formula,data=model$data)
    xmat <- naresid(NULL, xmat) #Changed, probably not working
    if (any(alias <- is.na(coef(model)))) 
      xmat <- xmat[, !alias, drop = FALSE]
    wts <- weights(model)
    if (is.null(wts)) 
      wts <- 1
    res <- model$get_residuals()[[1]]
    rval <- as.vector(res) * wts * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
  } 
  else if (model$fn=="stats::glm") {
    xmat <- model.matrix(model$formula,data=model$data)
    if(any(alias <- is.na(coef(model)))) xmat <- xmat[, !alias, drop = FALSE]
    wres <- as.vector(residuals(model, "working")) * weights(model, "working")
    dispersion <- if(substr(x$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
    else sum(wres^2, na.rm = TRUE)/sum(weights(x, "working"), na.rm = TRUE)
    rval <- wres * xmat / dispersion
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    res <- residuals(x, type = "pearson")
    if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
    if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
    return(rval)
  }
  
  
  
}

model <- model1_logit
cluster <- flights2$origin

robcov_zelig <- function(model,cluster) {
var <- vcov(model)[[1]]
vname <- dimnames(var)[[1]]
if (model$fn=="stats::lm")
  var <- model$get_df_residual()[[1]] * var/sum(model$get_residuals()[[1]]^2)
X <- estfun_zelig(model)
n <- nrow(X)
if (any(is.na(cluster))) 
  stop("cluster contains NAs")
clusterInfo <- list(name = deparse(substitute(cluster)))
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
W <- matrix(.Fortran(rms:::F_robcovf, n, p, nc, clus.start, clus.size, 
                     X, double(p), double(p * p), w = double(p * p))$w, nrow = p)
adjvar <- var %*% W %*% var
fit$orig.var <- var
fit$var <- adjvar
fit$clusterInfo <- clusterInfo
fit
}

fit <- model_lrm

robcov <- function (fit, cluster, method = c("huber", "efron")) {
  method <- match.arg(method)
  var <- vcov(fit, intercepts = "all")
  vname <- dimnames(var)[[1]]
  if (inherits(fit, "ols")) 
    var <- fit$df.residual * var/sum(fit$residuals^2)
  else if (method == "efron") 
    stop("method=\"efron\" only works for ols fits")
  X <- as.matrix(residuals(fit, type = if (method == "huber") 
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

#Problem is estfun works differently for each model...Need to make methods for each type of zelig function...


model1 %>% robcov_zelig(flights2$origin)
model_ols %>% robcov(flights2$origin)



model_lrm %>% robcov(flights2$origin)

model1_logit %>% robcov_zelig(flights2$origin) #Wrong
