#' scABC automated pipeline
#' @param bamfiles location of bam files
#' @param peakfile location of peak file
#' @param plot boolean variable to make plots
#' @param nClusters number of clusters, if not set manually then the number of clusters is determined by gap stat analysis
#' @import WeightedCluster 
#' @import Rsamtools
#' @import GenomicRanges
#' @output list of cluster assigments, landmarks, peak p-values
#' 
scABC <- function(bamfiles, peakfile, plot = FALSE, nClusters = 0, pValThresh = 2, 
                  nreadsThresh = 2, ncellsThresh = 10, medianBGthresh = 2, 
                  lambda = 1, nTop = 2000){
  peaks = select_peaks(peakfile,thresh = pValThresh)
  ForeGround = get_counts_matrix(bamfiles, peaks)
  ForeGroundFiltered = filter_peaks(ForeGround$ForeGround, ForeGround$peaks,
                                    nreads_thresh = nreadsThresh, ncells_thresh = ncellsThresh)
  peaks = ForeGroundFiltered$peaks
  BackGround = get_background(bamfiles, peaks)
  ForeGroundBackGroundFiltered = filter_background(ForeGround = ForeGroundFiltered$ForeGround, 
                                                   BackGround = BackGround$BackGround, 
                                                   median_bg_thresh = medianBGthresh)
  ForeGroundMatrix = ForeGroundBackGroundFiltered$ForeGround
  BackGroundMatrix = ForeGroundBackGroundFiltered$Background
  if(nClusters > 0){
    InSilicoLandMarks = compute_landmarks(ForeGround = ForeGroundMatrix, 
                                          BackGround = BackGroundMatrix, 
                                          nCluster = nCluster, lambda = lambda, nTop = nTop)
  }
  else{ #nClusters unset, so comput nClusters by GapStat
    
  }
}