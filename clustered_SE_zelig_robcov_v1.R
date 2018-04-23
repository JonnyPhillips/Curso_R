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


model_ols %>% robcov(flights2$origin)

fit <- model_ols
cluster <- flights2$origin
method <- "huber"

estfun_zelig <- function (model, ...) {
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

model1 %>% robcov_zelig(flights2$origin)
model_ols %>% robcov(flights2$origin)



model_lrm %>% robcov(flights2$origin)

model1_logit %>% robcov_zelig(flights2$origin) #Wrong
