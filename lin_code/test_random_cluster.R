n <- 10^5
ks <- round(10^seq(1, 4, 0.2))

nrun <- 50
propA <- matrix(nrow=nrun, ncol=length(ks))
for (j in 1:length(ks)){
  k <- ks[j]
  for (i in 1:nrun){
    d1 <- sample(n, k)
    d2 <- sample(n, k)
    propA[i,j] <- (2*k-length(unique(c(d1, d2))))/k
  }
}
boxplot(propA[,6:13])


nrun <- 100
statA <- matrix(nrow=nrun, ncol=length(ks))
for (j in 1:length(ks)){
  k <- ks[j]
  for (i in 1:nrun){
    d1 <- sample(n, k)
    d2 <- sample(n, k)
    p <- k/n*k/n
    statA[i,j] <- ((2*k-length(unique(c(d1, d2)))) - n*p)/sqrt(n*p*(1-p))
  }
}
boxplot(statA)

nrun <- 100
statA <- matrix(nrow=nrun, ncol=length(ks))
for (j in 1:length(ks)){
  k <- ks[j]
  for (i in 1:nrun){
    d1 <- sample(n, k)
    d2 <- sample(n, k)
    p <- k/n*k/n
    statA[i,j] <- ((2*k-length(unique(c(d1, d2)))) - n*p)/sqrt(n*p*(1-p))
  }
}

rd <- sapply(numregions[which(upper.tri(numregions))], function(k, n){return((2*k-length(unique(c(sample(n, k), sample(n, k)))))/k)}, n)
rd <- 1-rd

hist(distm[which(upper.tri(distm))], breaks=200)
hist(rd, breaks=200, add=T, col="red")

rdm <- 0*distm
rdm[which(upper.tri(rdm))] <- rd
rdm <- rdm + t(rdm)

## some artifitial clusters
power <- 30
heatmap.2(rdm^power, Colv="Rowv",trace='none')
