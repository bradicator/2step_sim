library(MASS)

# Large Scale Simulation
# setup
M <- 100 
p <- 6
set.seed(123)
count <- round(runif(100, 30, 300))
betamat <- 2 * array(c(1:6, 6:1, c(1,3,5,5,3,1)), c(6,3))
Sigma1 <- diag(p)
Sigma2 <- 2*Sigma1
Sigma3 <- 0.9^abs(outer(1:p, 1:p, "-")) 
Sigma4 <- 2 * Sigma3

Blist <- vector("list", M)
for (i in 1:M){
  Blist[[i]] <- array(runif(3*p, -3, 3), c(p, 3))
}

xlist <- vector("list", M)
for (i in 1:M){
  xlist[[i]] <- array(runif(3*count[i], -3, 3), c(count[i], 3))
}

ylist <- vector("list", M)
for (i in 1:M){
  if (i < 26) temp <- Sigma1
  if (i > 25 & i < 51) temp <- Sigma2
  if (i > 50 & i < 76) temp <- Sigma3
  if (i > 75) temp <- Sigma4
  err <- t(mvrnorm(count[i], rep(0, p), temp))
  ylist[[i]] <- (betamat + Blist[[i]]) %*% t(xlist[[i]]) + err
}

# add column
for (i in 1:M){
  xlist[[i]] <- cbind(rep(1, count[i]), xlist[[i]])
}

# run our procedure
result <- new_shrinkage(xlist, ylist, 5000, 2000, c(2,3,4))
result2 <- new_shrinkage(xlist, ylist, 5000, 5, c(2,3,4))


#
beta_dis <- rep(0, 81)
for (i in 1:81){
  beta_dis[i] <- log(norm(result$betahis[[i]] - result$betahis[[81]], "F"))
}
plot(1:81, beta_dis)



# MSE
cons <- -betamat+result$betamat[,2:4]
mse <- function(list1, list2, cons){
  n <- length(list1)
  sse <- 0
  for (i in 1:n){
    sse <- sse + norm(list1[[i]]-list2[[i]]+cons, "F")^2
  }
  sse/(3*n)
}

mse(result$Blist, Blist, cons)

ind <- which(count >= 200)
mse(result$Blist[ind], Blist[ind], cons)

ind <- which(count < 40)
mse(result$Blist[ind], Blist[ind], cons)

ind <- which(count >= 40 & count < 100)
mse(result$Blist[ind], Blist[ind], cons)

ind <- which(count >= 100 & count < 200)
mse(result$Blist[ind], Blist[ind], cons)
# mse for covariance structure.
sum(sapply(sig$Sigma_output[1:25], function(x) norm(x-Sigma1, "F")^2))/25
sum(sapply(sig$Sigma_output[26:50], function(x) norm(x-Sigma2, "F")^2))/25
sum(sapply(sig$Sigma_output[51:75], function(x) norm(x-Sigma3, "F")^2))/25
sum(sapply(sig$Sigma_output[76:100], function(x) norm(x-Sigma4, "F")^2))/25



ga <- 5000/ result$count
ind <- which(result$count > 120)
ga[ind] <- 0
sig <- optimize_sigma(result, ga, F)



# covariance recovery? mds
# mds of all sigma matrices. import mds fun from pilotstudy.r
sigout <- c(sig$Sigma_output, list(Sigma1, Sigma2, Sigma3, Sigma4, result$Sigma))
fit <- mds(sigout)
x <- fit$points[,1]
y <- fit$points[,2]
col <- as.character(c(rep(1:4, each = 25),rep("Benchmark", 5)))

aa <- ggplot(aes(x=x, y=y, colour = col, size = size, shape = sp), data = data.frame(x=x,y=y,col=col,size=c(count,rep(300,5)),sp = c(rep("a",100), rep("b",5))), 
       size = size) +geom_point() +guides(shape=F) + ggtitle(expression(paste("MDS of ", hat(Sigma[i]))))
aa<-aa + scale_colour_discrete(name="Real Covariance",labels = c(expression(paste(Sigma[i]," = I")),expression(paste(Sigma[i]," = 2I")),expression(paste(Sigma[i]," = R")),
                                    expression(paste(Sigma[i]," = 2R")),"Benchmark"))
aa <- aa + theme(legend.text.align=0)





# fit one single driver at a time 
oned_record <- rep(0, 100)
for (i in 1:100){
  special_dri <- i
  res2 <- new_shrinkage(xlist[special_dri], ylist[special_dri], 0, 3000, c(2,3,4))
  x1 <- res2$betamat[,2:4] + res2$Blist[[1]]
  xt <- betamat + Blist[[special_dri]]
  oned_record[i] <- norm(xt-x1, "F")^2
}

# MSE
cons <- -betamat+result$betamat[,2:4]
mse <- function(list1, list2, cons){
  n <- length(list1)
  sse <- 0
  for (i in 1:n){
    sse <- sse + norm(list1[[i]]-list2[[i]]+cons, "F")^2
  }
  sse/(3*n)
}

mse(result$Blist, Blist, cons)
mean(oned_record)/3

ind <- which(count >= 200)
mse(result$Blist[ind], Blist[ind], cons)
mean(oned_record[ind])/3

ind <- which(count < 40)
mse(result$Blist[ind], Blist[ind], cons)
mean(oned_record[ind])/3

ind <- which(count >= 40 & count < 100)
mse(result$Blist[ind], Blist[ind], cons)
mean(oned_record[ind])/3

ind <- which(count >= 100 & count < 200)
mse(result$Blist[ind], Blist[ind], cons)
mean(oned_record[ind])/3




## two rounds of fitting
r2res <- cov_to_mean(xlist, ylist, 5000, 3000, result, sig,c(2,3,4))

cons <- -betamat+r2res$betamat[,2:4]

mse(r2res$Blist, Blist, cons)

ind <- which(count >= 200)
mse(r2res$Blist[ind], Blist[ind], cons)

ind <- which(count < 40)
mse(r2res$Blist[ind], Blist[ind], cons)

ind <- which(count >= 40 & count < 100)
mse(r2res$Blist[ind], Blist[ind], cons)

ind <- which(count >= 100 & count < 200)
mse(r2res$Blist[ind], Blist[ind], cons)