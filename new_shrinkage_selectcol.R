new_shrinkage <- function(xlist, ylist, lambda_input, iter, Bselect, if_count = T){
  
  # initialization
  nvar <- ncol(xlist[[1]])
  p <- nrow(ylist[[1]])
  n <- length(xlist)
  count <- sapply(ylist, ncol)
  
  betamat <- array(1, c(p, nvar))
  Blist <- vector("list", n)
  
  if (if_count){
    lambda <- lambda_input / count
  } else{
    lambda <- rep(lambda_input, n)
  }
  
  for (i in 1:n) Blist[[i]] <- betamat[, Bselect]
  
  Sigma = diag(p)
  
  for (i in 1:iter){
    print(i)
    # update Blist
    for (j in 1:length(Bselect)){
      for (k in 1:n){
        #print(k)
        #print(betamat)
        #print(t(xlist[[k]]))
        temp <- ylist[[k]] - betamat %*% t(xlist[[k]]) - 
                Blist[[k]][, -j, drop = F] %*% t(xlist[[k]][, Bselect, drop = F][, -j, drop = F])
        temp2 <- temp %*% xlist[[k]][, Bselect, drop = F][, j, drop = F]
        temp_A <- lambda[k]*Sigma + sum(xlist[[k]][, Bselect, drop = F][, j]^2)*diag(p)
        #print(Blist[[k]][, j, drop = F])
        #print(solve(temp_A, temp2))
        Blist[[k]][, j] <- solve(temp_A, temp2)       
      }
    }
    
    # update betamat
    for (j in 1:nvar){
      temp3 <- array(0, c(p, 1))
      temp4 <- 0 
      for (k in 1:n){
        temp <- ylist[[k]] - betamat[, -j, drop = F] %*% 
                t(xlist[[k]][, -j, drop = F]) - Blist[[k]] %*% 
                t(xlist[[k]][, Bselect, drop = F])
        temp2 <- temp %*% xlist[[k]][, j, drop = F]
        temp3 <- temp3 + temp2
        temp4 <- temp4 + sum(xlist[[k]][,j]^2)
      }
      betamat[, j] <- temp3 / temp4
    }
    
    # update Sigma
    Slist <- vector("list", n)
    tempres <- array(0, c(p, p))
    count <- 0
    for (k in 1:n){
      temp <- ylist[[k]] - betamat %*% t(xlist[[k]]) - 
              Blist[[k]] %*% t(xlist[[k]][, Bselect, drop = F])
      Slist[[k]] <- temp %*% t(temp)
      tempres <- tempres + Slist[[k]]
      count <- count + ncol(ylist[[k]])
    }
    Sigma <- tempres / count
    
    # calculate loglikehood
    loglik <- count * log(det(Sigma))
    
    for (item in Blist) loglik <- loglik + lambda[k] * norm(item, type="F")^2
    print(loglik)
    if (i > 1){
      if (loglik_his - loglik < 0.01) break
    }
    loglik_his <- loglik
    # print(Sigma)
  }
  
  
  return(list(Blist = Blist, betamat = betamat, loglik = loglik,
              Sigma = Sigma, Slist = Slist, count = sapply(ylist, ncol)))
}