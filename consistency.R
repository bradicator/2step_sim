library(MASS)

# Study if b and beta are consistent
# 1. always converges
# 2. always converge to the same solution: check betamat diff and b[1]
# diff
# setup
consistency <- function(M, lb, ub){
  p <- 6
  set.seed(123)
  count <- round(runif(M, lb, ub))
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
  
  # run procedure
  result <- new_shrinkage(xlist, ylist, 5000, 2000, c(2,3,4))
  betamat0 <- array(c(rep(0, 6), 1:6, 6:1, c(1,3,5,5,3,1)), c(6,4))
  beta_dis <- norm(result$betamat - betamat0, "F")^2
  
  b_dis <- 0 
  for (i in 1:M){
    b_dis <- b_dis + norm(result$Blist[[i]][, 1, drop = F] + result$betamat[,2, drop = F] - 
                            Blist[[i]][, 1, drop = F] - betamat[, 1, drop = F], "F")^2
  }
  b_dis <- b_dis / M
  
  return(list(beta_dis = beta_dis, b_dis = log(b_dis)))
}

# change no of drivers
beta_dis_rec <- b_dis_rec <- rep(0, 7)
Mvec <- c(100, 300, 500, 1000, 2000, 4000, 8000)
for (i in 1:7){
  temp <- consistency(Mvec[i], 30, 300)
  beta_dis_rec[i] <- temp$beta_dis
  b_dis_rec[i] <- temp$b_dis
}


# change obs per driver
beta_dis_rec <- b_dis_rec <- rep(0, 5)
Nvec <- c(0, 120, 320, 820, 1820)
for (i in 1:5){
  temp <- consistency(40, 30 + Nvec[i], 330 + Nvec[i])
  beta_dis_rec[i] <- temp$beta_dis
  b_dis_rec[i] <- temp$b_dis
}
