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
#' @param lambda weighting parameter (default = 1)
#' @param nTop number of top clusters
#' @import WeightedCluster 
#' @export
compute_landmarks <- function(ForeGround, BackGround, nCluster = 2, lambda = 1, nTop = 2000){
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
  q = median(BGmedian);
  lambda = 1;
  weights = 1/(1 + exp(-lambda*(BGmedian - q)));
  kMeds = wcKMediods(FGdist, k = nCluster, weights = weights);
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
#' @export
assign2landmarks <- function(ForeGround, topLandmarks){
 scor = cor(ForeGround, topLandmarks, method = "spearman")
 return(apply(scor, 1, which.max))
}


