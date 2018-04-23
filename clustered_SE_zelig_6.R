library(rms)
library(tidyverse)
library(zeligverse)
library(nycflights13)

###
flights2 <- flights %>% filter(!is.na(dep_delay) & !is.na(dep_time) & !is.na(origin))

model1 <- flights2 %>% zelig(dep_delay~dep_time + origin,data=.,model="ls")
summary(model1)

model_lm <- flights2 %>% lm(dep_delay~dep_time + origin,data=.)

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


zelig_cluster_se <- function(model,cluster){
  
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
    if (is.zoo(res)) 
      rval <- zoo(rval, index(res), attr(res, "frequency"))
    if (is.ts(res)) 
      rval <- ts(rval, start = start(res), frequency = frequency(res))
    return(rval)
  }
  
  sandwich_zelig <- function(model,meat){
    estd <- estfun_zelig(model)
    K <- dim(vcov(model1)[[1]])[1]
    sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/dim(estd)[1]-K)#Double check, but about right
    bread <- (model$get_vcov()[[1]]/sigma^2)*dim(estd)[1]#dim not exact - 4 rows duplicated for some reason...
    n <- nrow(estd)
    return(1/n * (bread %*% meat %*% bread))
  }
  
  cl_zelig   <- function(model, cluster){
    require(sandwich, quietly = TRUE)
    require(lmtest, quietly = TRUE)
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- dim(vcov(model1)[[1]])[1] #Probably wrong!
    dfc <- (M/(M-1))*((N-1)/(N-K))
    uj  <- apply(estfun_zelig(model),2, function(x) tapply(x, cluster, sum));
    vcovCL <- dfc*sandwich_zelig(model, meat=crossprod(uj)/N)
    coeftest(model, vcovCL) 
  }

  cl_zelig(model,cluster)
  
}

flights2 %>% zelig(dep_delay~dep_time + origin,data=.,model="ls") %>% 
  zelig_cluster_se(flights2$origin)

flights2 %>% mutate(Late=ifelse(dep_delay>0,1,0)) %>% 
  ols(dep_delay~dep_time + origin,data=.,x=T,y=T) %>% 
  robcov(flights2$origin)

flights2 %>% mutate(Late=ifelse(dep_delay>0,1,0)) %>% 
  zelig(Late~dep_time + origin,data=.,model="logit") %>% 
  zelig_se(flights2$origin)


#Works but seems to produce smaller SEs so check carefully... Better to integrate my adaptations into robcov function


flights2 %>% mutate(Late=ifelse(dep_delay>0,1,0)) %>% 
  lrm(Late~dep_time + origin,data=.,x=T,y=T) %>% 
  robcov(flights2$origin)


zelig_se
cl(flights2,model_lm,flights2$origin)
#Pretty close
