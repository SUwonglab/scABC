
select_top <-function(x, n_top){
  thresh = sort(x, decreasing = TRUE)[n_top]
  x[x < thresh] = 0;
  return(x);
}

#' compute landmarks
#' @name compute_landmarks
#' @param ForeGround matrix or data frame of Foreground values
#' @param BackGround matrix or data frame of BackGround values
#' @param nCluster number of clusters (default = 2)
#' @param lambda weighting parameter (default = 0.1)
#' @param nTop number of top clusters
#' @import WeightedCluster 
#' @export compute_landmarks
compute_landmarks <- function(ForeGround, BackGround, nCluster = 2, lambda = 0.1, nTop = 2000){
  # check types
  stopifnot(is.matrix(ForeGround) || is.data.frame(ForeGround))
  stopifnot(is.matrix(BackGround) || is.data.frame(BackGround))
  stopifnot(is.numeric(nCluster))
  # ensure both are matrices
  ForeGround = as.matrix(ForeGround);
  BackGround = as.matrix(BackGround);
  # compute weighted k-medioids
  FGdist = 1 - cor(ForeGround, method = "spearman");
  BGmedian = apply(BackGround, 2, median);
  c = min(8, median(BGmedian));
  weights = 1/(1 + exp(-(BGmedian - c)/(c*lambda)));
  kMeds = wcKMedoids(FGdist, k = nCluster, weights = weights);
  kMedsCluster = kMeds$clustering;
  kMedsCluster = as.numeric(factor(kMedsCluster));

  landmarks = c()
  for(i in 1:nCluster){
    tmp = which(kMedsCluster == i);
    if(length(tmp) == 1){
      landmarks = cbind(landmarks, ForeGround[ ,tmp]);
    }
    else{
      landmarks= cbind(landmarks, rowSums(ForeGround[ ,tmp]));
    }
  }

  topLandmarks = apply(landmarks, 2, select_top, nTop);

  return(topLandmarks);
}

#' assign cells/samples to the closest package
#' 
#' @docType package
#' @name assign2landmarks
#' @param ForeGround matrix or data frame of Foreground values
#' @param topLandmarks output from compute_landmarks
#' @export assign2landmarks
assign2landmarks <- function(ForeGround, topLandmarks){
 scor = cor(ForeGround, topLandmarks, method = "spearman")
 return(apply(scor, 1, which.max))
}

#' get cluster specific pvals for peaks
#' 
#' @docType package
#' @name getClusterSpecificPvalue
#' @param data matrix of peaks by cells counts
#' @param cluster matrix of cluster indicator membership 
#' @param background_median median background values for each cell
#' @param landmark optional. If landmark is provided, we only do hypothesis testing on the union of landmark peaks
getClusterSpecificPvalue <- function(ForeGround, cluster_assignments, background_medians, 
                                     landmark=NULL, maxiter=1000, thresMLE=10^-3, 
                                     thresMAP=10^-5, quiet=TRUE){
  ## the main function for peak selection
  ## data is nPeaks(p) by Cells(n)
  ## cluster_assignments is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## background_medians corresponds to h_i, here it is a vector of length n
  ## landmark is optional. If landmark is provided, we only do hypothesis testing on the union of landmark peaks
  if (!is.null(landmark)){
    ## take union of the landmark peaks
    unionlandmark <- which(rowSums(landmark)!=0)
    p <- nrow(ForeGround)
    ForeGround <- ForeGround[unionlandmark,]
  } else {
  }
  
  ForeGround <- t(ForeGround) # make data n by p
  
  ## get beta_MLE
  if (!quiet){
    cat("\nEstimating beta MLE\n")  
  }
  betaMLE_ini <- matrix(0, nrow=length(unique(cluster_assignments)), ncol=ncol(ForeGround))
  betaMLE <- getbetaMLE(data=ForeGround, cluster=cluster_assignments, background_medians=background_medians, beta_ini=betaMLE_ini, maxiter=maxiter, thres=thresMLE, quiet=quiet)
  
  ## get the empirical prior sigma
  sigmas <- getSigmaPrior(betaMLE)
  
  ## get beta_MAP and the p-value
  if (!quiet){
    cat("\nEstimating beta MAP\n")
  }
  betaMAP_ini <- rbind(0, matrix(0, nrow=length(unique(cluster_assignments)), ncol=ncol(ForeGround)))
  result <- getbetaMAP(data=ForeGround, cluster=cluster_assignments, background_medians=background_medians, sigmas=sigmas, beta_ini=betaMAP_ini, maxiter=maxiter, thres=thresMAP, quiet=quiet)
  
  if (!is.null(landmark)){
    betaMAP <- matrix( nrow=length(unique(cluster_assignments))+1, ncol=p )
    pvalue <- matrix( nrow=p, ncol=length(unique(cluster_assignments)) )
    betaMAP[, unionlandmark] <- result$beta
    pvalue[unionlandmark,] <- result$pvalue
  } else {
    betaMAP <- result$beta
    pvalue <- result$pvalue
  }
  colnames(pvalue) <- paste("Cluster", 1:length(unique(cluster_assignments)) )
  return(list(betaMAP=betaMAP, sigmas=sigmas, pvalue=pvalue))
}


getbetaMLE <- function(data, cluster, background_medians, 
                       beta_ini = 1, maxiter = 1000, 
                       thres = 1e-5, quiet = TRUE){
  ## data is n by p
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## background_medians corresponds to h_i, here it is a vector of length n
  p <- ncol(data)
  x <- getdesign(cluster)
  beta <- beta_ini
  betaDiffs <- c() 
  percConverged <- c()
  converged_flag <- rep(0, p)
  
  converged <- 0
  iter <- 1
  while (!converged & iter<maxiter){
    mu <- background_medians*exp(x%*%beta) 
    u <- crossprod(x, data-mu)
    updateBeta <- function(r){
      if (converged_flag[r]==0 & !is.na(beta[1,r]) ){
        Jr <- crossprod(x, mu[,r]*x)
        betaDiff <- try(solve(Jr, u[,r]), silent=T)
        if (class(betaDiff)=="try-error"){
          return( rep(NA, length(beta[,r])) )    
        } else {
          return(beta[,r] + betaDiff)  
        }
      } else {
        return(beta[,r])
      }
    }
    
    betaPre <- beta
    beta <- sapply(1:p, updateBeta)
    betaDiffs <- rbind(betaDiffs, colSums(abs(beta - betaPre)) )
    converged_flag <- (betaDiffs[iter,] <= thres) + 0
    converged_flag[which(is.na(converged_flag))] <- 0
    singular_flag <- is.na(beta[1,]) + 0
    percConverged <- c(percConverged, sum(converged_flag[which(singular_flag==0)])/sum(singular_flag==0) )
    
    if (!quiet){
      cat("\r", round(percConverged[length(percConverged)]*100), "%", "converged")  
    }
    if (percConverged[length(percConverged)]==1){
      converged <- 1
    }
    iter <- iter + 1
  }
  return(beta)
}

getdesign <- function(cluster){
  nCluster <- length(unique(cluster))
  n <- length(cluster)
  x <- matrix(0, nrow=n, ncol=nCluster)
  for (i in 1:nCluster){
    x[which(cluster==i), i] <- 1
  }
  return(x)
}


getSigmaPrior <- function(beta, q=0.05){
  ## get the empirical prior estimate using quantile matching
  ## q is the quantile to match
  
  ## center each column in beta
  betaC <- t(t(beta) - colMeans(beta))
  ## 
  sigmas <- rep(0, nrow(betaC))
  for (i in 1:nrow(betaC)){
    tmp <- betaC[i,]
    tmp <- tmp[which(!is.na(tmp))]
    sigmas[i] <- quantile(abs(tmp), 1-q)/qnorm(1-q/2)
  }
  return(sigmas)
}

getbetaMAP <- function(data, cluster, background_medians, sigmas, beta_ini, maxiter, thres, quiet){
  ## data is n by p
  ## beta_ini includes the intercept
  ## cluster is the cluster membership matrix, length n. Take entries 1 to nCluster
  ## background_medians corresponds to h_i, here it is a vector of length n
  
  p <- ncol(data)
  x <- getdesign(cluster)
  ## add the intercept in x
  x <- cbind(1, x)
  ## get lambda
  lambda <- c(0, 1/sigmas^2)
  
  beta <- beta_ini
  betaDiffs <- c()
  percConverged <- c()
  converged_flag <- rep(0, p)
  
  converged <- 0
  iter <- 1
  while (!converged & iter<maxiter){
    mu <- background_medians*exp(x%*%beta) 
    z <- log(mu/background_medians) + (data-mu)/mu
    
    updateBetaRidge <- function(r){
      if (converged_flag[r]==0){
        JrRidge <- crossprod(x, mu[,r]*x) + diag(lambda)
        betar <- solve(JrRidge, crossprod(x, matrix(mu[,r]*z[,r], ncol=1)) )
        return( as.vector(betar) )
      } else {
        return(beta[,r])
      }
    }
    
    betaPre <- beta
    beta <- sapply(1:p, updateBetaRidge)
    betaDiffs <- rbind(betaDiffs, colSums(abs(beta - betaPre)) )
    converged_flag <- (betaDiffs[iter,] <= thres) + 0
    percConverged <- c(percConverged, sum(converged_flag)/p )
    
    if (!quiet){
      cat("\r", round(percConverged[length(percConverged)]*100), "%", "converged")  
    }
    
    if (percConverged[length(percConverged)]==1){
      converged <- 1
    }
    iter <- iter + 1
  }
  mu <- background_medians*exp(x%*%beta) 
  
  ## get all the contrast matrices and cluster label for each contrast
  nCluster <- length(unique(cluster))
  tmp <- getContrast(nCluster)
  constrasts <- tmp$constrasts
  
  ## pad the intercept with 0 in the constrast
  constrasts <- cbind(0, constrasts)
  constrastsCluster <- tmp$constrastsCluster
  
  calRidgePvalueA <- function(r){
    ## calculates the p-value for a small hypothesis
    if ( !is.na(beta[1,r]) ){
      Jr <- crossprod(x, mu[,r]*x)
      JrRidgeInv <- solve(Jr + diag(lambda))
      covRidge <- JrRidgeInv%*%Jr%*%JrRidgeInv
      ##
      betaC <- constrasts%*%matrix(beta[,r], ncol=1)  
      CcovC <- constrasts%*%covRidge%*%t(constrasts)
      SEbetaC <- sqrt(diag(CcovC))
      ##
      pvs <- pnorm(betaC/SEbetaC, lower.tail = FALSE)
      return(pvs)
    } else {
      return(rep(NA, length(constrastsCluster)))
    }
  }
  if (!quiet){
    cat("\nCalculating the p-values\n")
  }
  pvsA <- sapply(1:p, calRidgePvalueA) 
  ## get the p-value for the large null hypothesis by taking maximum
  if (nCluster>2){
    pvs <- c()
    for (i in 1:nCluster){
      pvs <- cbind(pvs, apply(pvsA[which(constrastsCluster==i),], 2, max))
    }  
  } else {
    pvs <- t(pvsA)
  }
  return(list(beta=beta, pvalue=pvs))
}

getContrast <- function(nCluster){
  ## get all the contrasts that is needed to calculate the p-value
  constrasts <- matrix(0, nrow=(nCluster-1)*nCluster, ncol=nCluster)
  constrastsCluster <- rep(1:nCluster, each=nCluster-1)
  for (i in 1:nCluster){
    constrasttmp <- matrix(0, nrow=nCluster-1, ncol=nCluster)
    constrasttmp[, i] <- 1
    count <- 1
    for (j in 1:(nCluster-1)){
      constrasttmp[count, (c(1:nCluster)[-i])[j]] <- -1
      count <- count + 1
    }
    constrasts[(i*(nCluster-1)-(nCluster-1)+1):(i*(nCluster-1)), ] <- constrasttmp
  }
  return(list(constrasts=constrasts, constrastsCluster=constrastsCluster))
}



getGapStat <- function(ForeGround, BackGroundMedian, nClusters=1:10, 
                       nPerm=10, nTop = 10000, quiet=FALSE){
  ## sort nClusters from small to large
  nClusters <- sort(nClusters, decreasing=FALSE)
  
  ##
  lambda <- 0.1
  c <- quantile(BackGroundMedian, 0.5)
  W <- 1/(1+exp(-(BackGroundMedian-c)/(lambda*c)))
  distS <- 1-cor(ForeGround[head(order(apply(ForeGround, 1, sum), decreasing = TRUE), nTop), ], 
                 method="spearman")
  
  ## function for calculating the objective
  calObj <- function(cluster, distMatrix, weight, nTop){
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

plotGapStat <- function(GapStat, nClusters, main, fold=1, cex=0.3, lty=2){
  ObjPerm <- GapStat$ObjPerm
  ObjData <- GapStat$ObjData
  ymean <- apply(log(ObjPerm), 2, mean) - log(ObjData)
  ysd <- apply(log(ObjPerm), 2, sd)
  
  plot(nClusters, ymean,
       ylim=range(c(ymean-ysd, ymean+ysd)),
       pch=19, xlab="Number of clusters", ylab="Gap", cex=cex,
       main=main)
  lines(nClusters, ymean, col="red", lty=lty)
  # hack: we draw arrows but with very special "arrowheads"
  arrows(nClusters, ymean-ysd, nClusters, ymean+ysd, length=0.05, angle=90, code=3)
  return(1)
}
