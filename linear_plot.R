


P<-50

lambda <- rep(Inf, P)
mu <- rep(Inf, P)

for( i in 1:P ){
 
    mu[i] <- (P - i) / P # linear
  
  lambda[i] <- 1 - mu[i]
  
}

par( cex =.9, col="blue")

plot(1:P, lambda, type="o", xlab="rank of pop", ylab = "lambda", main = "liner migration curve")
plot(1:P, mu, type="o", xlab="rank of pop", ylab = "mu", main = "liner migration curve")



norm_rank <- rep(Inf, P)
for( i in 1:P ){
  
  norm_rank[i] <- (i - 1) / (P-1) # linear
}

#par( cex =.9, col="blue")
#plot(1:50, tuned_val_a, type="o", xlab="rank of pop", ylab = "lambda", main = "evolved migration curve")
#par( cex =.9, col="blue")
#plot(1:50, tuned_val_b, type="o", xlab="rank of pop", ylab = "mu", main = "evolved migration curve")