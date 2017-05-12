getGapStatLin <- function(ForeGround, BackGroundMedian, nClusters=1:10, nPerm=20, quiet=FALSE, topPeak=5000){
  ## sort nClusters from small to large
  nClusters <- sort(nClusters, decreasing=FALSE)
  if (nrow(ForeGround)>=topPeak){
    ForeGround <- ForeGround[order(rowSums(ForeGround), decreasing=T)[1:topPeak],]
  }
  ##
  lambda <- 0.1
  c <- min(8, quantile(BackGroundMedian, 0.5))
  W <- 1/(1+exp(-(BackGroundMedian-c)/(lambda*c)))
  distS <- 1-cor(ForeGround, method="spearman")
  
  ## function for calculating the objective
  calObj <- function(cluster, distMatrix, weight){
    ## may encounter problem when there is singleton
    diag(distMatrix) <- 0
    medoids <- unique(cluster)  
    Obj <- 0
    for (medoid in medoids){
      disttmp <- distMatrix[medoid,which(cluster==medoid)] 
      weighttmp <- weight[which(cluster==medoid)]
      Obj <- Obj + sum(disttmp*weighttmp)
    }
    return(Obj)
  }
  
  ## turn off warning for wcKMedoids
  options(warn=-1)
  ##
  if (!quiet){
    cat("\nCalculate gap statistic for data\n")
  }
  ObjData_nClusters <- c()
  count <- 0 
  for (nCluster in nClusters){
    if (nCluster==1){
      ObjData <- min(rowSums(t(t(distS)*W)))
      ObjData_nClusters <- c(ObjData_nClusters, ObjData) 
    } else {
      resultW <- wcKMedoids(distS, k=nCluster, weights=W)
      clusterW <- resultW$clustering
      ObjData <- calObj(cluster=clusterW, distMatrix=distS, weight=W) 
      ObjData_nClusters <- c(ObjData_nClusters, ObjData) 
    }
    count <- count + 1
    if (!quiet){
      cat("\r", round((count/length(nClusters))*100), "%","completed")
    }
  }
  
  if (!quiet){
    cat("\nCalculate gap statistic for permutation\n")
  }
  
  ObjPerm_nClusters <- matrix(nrow=nPerm, ncol=length(nClusters))
  for (perm in 1:nPerm){
    ForeGroundPerm <- t(apply(ForeGround, 1, sample))
    distP <- 1-cor(ForeGroundPerm, method="spearman")
    ObjPA <- c()
    count <- 0
    for (nCluster in nClusters){
      count <- count + 1
      if (!quiet){
        cat("\r", round((perm-1+count/length(nClusters))/nPerm*100), "%", "completed")
      }
      if (nCluster==1){
        ObjP <- min(rowSums(t(t(distP)*W)))
        ObjPA <- c(ObjPA, ObjP)  
      } else {
        resultP <- wcKMedoids(distP, k=nCluster, weights=W)
        clusterP <- resultP$clustering
        ObjP <- calObj(cluster=clusterP, distMatrix=distP, weight=W) 
        ObjPA <- c(ObjPA, ObjP)  
      }
    }
    ObjPerm_nClusters[perm,] <- ObjPA
  }
  ## turn on warning
  options(warn=0)
  
  ## get the optimal number of clusters
  ymean <- apply(log(ObjPerm_nClusters), 2, mean) - log(ObjData_nClusters)
  ysd <- apply(log(ObjPerm_nClusters), 2, sd)
  
  flag <- ymean[-length(nClusters)] >= (ymean - ysd)[-1]
  if (sum(flag)){
    nClusterOptimal <- min(which(flag)) 
  } else {
    warning("Increase the maximum cluster number")
    nClusterOptimal <- max(nClusters)
  }
  
  return(list(ObjData=ObjData_nClusters, ObjPerm=ObjPerm_nClusters, nClusterOptimal=nClusterOptimal))
}