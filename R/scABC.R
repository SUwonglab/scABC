#' scABC automated pipeline
#' @param bamfiles location of bam files
#' @param peakfile location of peak file
#' @param plot boolean variable to make plots
#' @param nClusters number of clusters, if not set manually then the number of clusters is determined by gap stat analysis
#' @import WeightedCluster 
#' @import Rsamtools
#' @import GenomicRanges
#' @output list of filtered data and peaks, cluster assigments, landmarks, and peak p-values
#' 
scABC <- function(bamfiles, peakfile, PLOT = FALSE, QUIET = TRUE,
                  nClusters = c(), pValThresh = 2, 
                  nreadsThresh = 2, ncellsThresh = 10, medianBGthresh = 2, 
                  lambda = 1, nTop = 2000, nPerm = 20, 
                  maxiter=1000, thresMLE=10^-3, thresMAP=10^-5){
  peaks = select_peaks(peakfile,thresh = pValThresh)
  if(!QUIET){cat("\nreading in foreground\n")}
    
  ForeGround = get_counts_matrix(bamfiles, peaks)
  ForeGroundFiltered = filter_peaks(ForeGround$ForeGround, ForeGround$peaks,
                                    nreads_thresh = nreadsThresh, ncells_thresh = ncellsThresh)
  peaks = ForeGroundFiltered$peaks
  if(!QUIET){cat("\nreading in background\n")}
  BackGround = get_background(bamfiles, peaks)
  ForeGroundBackGroundFiltered = filter_background(ForeGround = ForeGroundFiltered$ForeGround, 
                                                   BackGround = BackGround$BackGround, 
                                                   median_bg_thresh = medianBGthresh)
  ForeGroundMatrix = ForeGroundBackGroundFiltered$ForeGround
  BackGroundMatrix = ForeGroundBackGroundFiltered$Background
  BackGroundMedian = apply(BackGroundMatrix, 2, median)
  if(!QUIET){cat("\ndone reading in data, beginning clustering\n")}
  if(length(nClusters) == 1){
    LandMarks = compute_landmarks(ForeGround = ForeGroundMatrix, 
                                  BackGround = BackGroundMatrix, 
                                  nCluster = nClusters, lambda = lambda, nTop = nTop)
  }
  else if(length(nClusters) > 1){
    GapStat = getGapStat(ForeGroundMatrix, BackGroundMedian, 
                         nClusters=nClusters, nPerm=nPerm, quiet=QUIET)
    if(PLOT){
      plotGapStat(GapStat, nClusters=nClusters, main = "Gap Stat")
    }
    nClusters = GapStat$nClusterOptimal
    LandMarks = compute_landmarks(ForeGround = ForeGroundMatrix, 
                                  BackGround = BackGroundMatrix,
                                  nCluster = GapStat$nClusterOptimal, 
                                  lambda = lambda, nTop = nTop)
  }
  else {
    GapStat = getGapStat(ForeGroundMatrix, BackGroundMedian, 
                         nClusters = 1:10, nPerm=nPerm, quiet = QUIET)
    if(PLOT){
      plotGapStat(GapStat, nClusters=nClusters, main = "Gap Stat")
    }
    nClusters = GapStat$nClusterOptimal
    LandMarks = compute_landmarks(ForeGround = ForeGroundMatrix, 
                                  BackGround = BackGroundMatrix,
                                  nCluster = nClusters, 
                                  lambda = lambda, nTop = nTop)
  }
  LandMarkAssignments = assign2landmarks(ForeGround = ForeGroundMatrix, LandMarks)
  PeakPvals = getClusterSpecificPvalue(data = ForeGround, cluster = LandMarkAssignments, 
                                       offset = BackGroundMedian, landmark=LandMarks, 
                                       maxiter=maxiter, thresMLE=thresMLE, thresMAP=10^-5, 
                                       quiet=QUIET)
  which_peaks = which(!is.na(PeakPvals$pvalue[,1]))
  
  return(list(ForeGroundMatrix = ForeGroundMatrix[which_peaks, ], peaks = peaks[which_peaks, ],
              nClusters = nClusters, cluster_assignments = LandMarkAssignments,
               LandMarks = LandMarks[which_peaks, ], PeakPVals = PeakPvals$pvalue[which_peaks, ]))
}