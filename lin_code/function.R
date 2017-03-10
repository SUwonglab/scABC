getClusterCount <- function(cluster, samples, cells){
  clusterCount <- matrix(0, nrow=length(cells), ncol=length(cells))
  for (i in 1:length(unique(cluster))){
    tmp <- samples[which(cluster==i)]
    for (j in 1:length(cells)){
      cell <- cells[j]
      clusterCount[i, j] <- sum(tmp==cell)
    }
  }
  ## rearrange it, sometimes get an error
  tmp <- apply(clusterCount, 2, which.max)
  if (max(table(tmp))>1){
    seqs <- c()
    for (i in 1:nrow(clusterCount)){
      if (i==1){
        seqs <- c(seqs, which.max(clusterCount[,i])[1])
      } 
      if (i>1 & i<nrow(clusterCount)){
        seqs <- c(seqs, c(c(1:nrow(clusterCount))[-seqs])[which.max(clusterCount[-seqs,i])[1]])  
      }
      if (i==nrow(clusterCount)){
        seqs <- c(seqs, c(1:nrow(clusterCount))[-seqs]) 
      }
    }
  } else {
    seqs <- tmp
  }
  clusterCount <- clusterCount[seqs,]
  row.names(clusterCount) <- paste("cluster", 1:6)
  colnames(clusterCount) <- cells
  return(clusterCount)
}