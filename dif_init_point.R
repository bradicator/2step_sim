library(MASS)

# Study the effect of different starting point
# 1. always converges
# 2. always converge to the same solution: check betamat diff and b[1]
# diff
# setup
M <- 100 
p <- 6
set.seed(123)
count <- round(runif(100, 30, 300))
betamat <- array(c(1:6, 6:1, c(1,3,5,5,3,1)), c(6,3))
Sigma1 <- diag(p)
Sigma2 <- 2*Sigma1
Sigma3 <- 0.9^abs(outer(1:p, 1:p, "-")) 
Sigma4 <- 2 * Sigma3

Blist <- vector("list", M)
for (i in 1:M){
  Blist[[i]] <- array(runif(3*p, -1, 1), c(p, 3))
}

xlist <- vector("list", M)
for (i in 1:M){
  xlist[[i]] <- array(runif(3*count[i], -3, 3), c(count[i], 3))
}

ylist <- vector("list", M)
for (i in 1:M){
  if (i < M/4 + 1) temp <- Sigma1
  if (i > M/4 & i <  M/2 + 1) temp <- Sigma2
  if (i > M/2 & i < 3*M/4 + 1) temp <- Sigma3
  if (i > 3*M/4) temp <- Sigma4
  err <- t(mvrnorm(count[i], rep(0, p), temp))
  ylist[[i]] <- (betamat + Blist[[i]]) %*% t(xlist[[i]]) + err
}

# add column
for (i in 1:M){
  xlist[[i]] <- cbind(rep(1, count[i]), xlist[[i]])
}

# 100 reps see if converge to the same point.
reptime <- 5
beta_his <- rep(0, reptime)
b_his <- rep(0, reptime)
result <- new_shrinkage(xlist, ylist, 5000, 2000, c(2,3,4))
betamat_0 <- result$betamat
bmat_0 <- result$Blist[[1]]
for (i in 1:reptime){
  result <- new_shrinkage(xlist, ylist, 5000, 2000, c(2,3,4))
  beta_his[i] <- log(norm(result$betamat - betamat_0, "F"))
  b_his[i] <- log(norm(result$Blist[[1]] - bmat_0, "F"))
}
