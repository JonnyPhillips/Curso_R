library(rms)
library(tidyverse)
library(zeligverse)
library(nycflights13)

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

sandwich_zelig <- function(model,meat){
    estd <- estfun_zelig(model)
    sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/sum(dim(estd))) #Double check, but about right
    bread <- model$get_vcov()[[1]]*sum(dim(estd))*sigma^2 #dim not exact - 4 rows duplicated for some reason...
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

cl_zelig(model1,flights2$origin)
cl(flights2,model_lm,flights2$origin)

dat <- flights2
fm <- model_lm
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M_z <- length(unique(cluster))
  N_z <- length(cluster)
  K_z <- fm$rank
  dfc_z <- (M_z/(M_z-1))*((N_z-1)/(N_z-K_z))
  uj_z  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
  vcovCL_z <- dfc_z*sandwich(fm, meat=crossprod(uj_z)/N_z)
  coeftest(fm, vcovCL_z)

model <- model1
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- dim(vcov(model1)[[1]])[1] #Probably wrong!
  dfc <- (M/(M-1))*((N-1)/(N-K))
  uj  <- apply(estfun_zelig(model),2, function(x) tapply(x, cluster, sum));
  vcovCL <- dfc*sandwich_zelig(model, meat=crossprod(uj)/N) #This line still wrong
  coeftest(model, vcovCL) 

#So need to compare sandwich(fm, meat=crossprod(uj_z)/N_z) with sandwich_zelig(model, meat=crossprod(uj)/N)

  
#Sandwich (verified):
  x <- model_lm
  meat_z <- crossprod(uj_z)/N_z
  bread_z <- bread(x)
  n_z <- NROW(estfun(x))
  1/n_z * (bread_z %*% meat_z %*% bread_z)
  
#Sandwich_zelig
  model <- model1
  meat <- crossprod(uj)/N
  estd <- estfun_zelig(model)
  sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/sum(dim(estd))) #Double check, but about right
  bread <- model$get_vcov()[[1]]*sum(dim(estd))*sigma^2 #dim not exact - 4 rows duplicated for some reason...
  n <- nrow(estd)
  1/n * (bread %*% meat %*% bread)
  
#So difference is just in bread
bread_z #should be

bread <- model$get_vcov()[[1]]*sigma^2

bread/bread_z #Fixed difference, but unclear where comes from...

#bread.lm
if (!is.null(x$na.action)) 
  class(x$na.action) <- "omit"
sx <- summary.lm(x)
return(sx$cov.unscaled * as.vector(sum(sx$df[1:2])))

sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/sum(dim(estd)))
(model$get_vcov()[[1]]/sigma^2)/summary(x)$cov.unscaled

summary(x)$df[1:2]

summary(x)$cov.unscaled* as.vector()

x <- model_lm
model <- model1
meat <- crossprod(uj)/N

sandwich <- function (x, bread. = bread, meat. = meat, ...) 
{
  if (is.list(x) && !is.null(x$na.action)) 
    class(x$na.action) <- "omit"
  if (is.function(bread.)) 
    bread. <- bread(x)
  if (is.function(meat.)) 
    meat. <- meat.(x, ...)
  n <- NROW(estfun(x))
  return(1/n * (bread. %*% meat. %*% bread.))
}

sx <- summary.lm(model_lm)
sx$cov.unscaled * as.vector(sum(sx$df[1:2]))
sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/dim(model$data)[[1]]) #+/-1 to df???
model$get_vcov()[[1]]/sigma^2 #Not exactly the same...hmmm...

sandwich_zelig <- function(model,meat){
  estd <- estfun_zelig(model)
  sigma <- sqrt(sum(model1$get_residuals()[[1]]^2)/dim(model$data)[[1]])
  bread_use <- model$get_vcov()[[1]]*(dim(estd)[[1]])/sigma^2 #dim not exact - 4 rows duplicated for some reason...
  n <- nrow(estd)
  return(1/n * (bread_use %*% meat %*% bread_use))
}
sandwich_zelig(model1,meat)
sandwich(model_lm,meat)
#Difference in bread - small but makes huge dif...