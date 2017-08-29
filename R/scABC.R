#' scABC
#' 
#' Core function of the scABC package.  Does processing, clustering, and peak p-value calculation in one step.
#' 
#' @param bamfiles location of bam files
#' @param peakfile location of peak file
#' @param byReadGroup boolean variable indicating whether or not individual experiments are separated by read group
#' @param PLOT boolean variable to make plots, default: TRUE
#' @param QUIET boolean variable for verbosity, default: TRUE
#' @param nClusters number of clusters, if not set manually then the number of clusters is determined by gap stat analysis
#' @param pValThresh the -log10(p) threshhold on the MACS2 computed peak p value, default: 2 (p < 0.01)
#' @param nreadsThresh threshold on the minimum number of reads in ncellsThresh to require when filtering cells, default: 2
#' @param ncellsThresh threshold on the minimum number of cells with nreadsThresh to require when filtering cell, default: 10
#' @param readsFGthresh threshold on the number of reads falling in peaks to filter cells by, default: min(500, nrow(ForeGround)/50)
#' @param lambda weighting parameter for clustering, default: 1
#' @param nTop the number of top peaks to use when forming the landmarks for each cluster, default: 2000
#' @param nPerm the number of permutation to use when selecting peaks by the gap statistic, default: 20
#' @param nTopGapStat the number of top peaks to use when computing the gap statistic, default: 5000
#' @param maxiter maximum number of iterations when computing beta, default: 1000
#' @param thresMLE numerical threshhold for convergence of beta MLE, default: 10^-3
#' @param thresMAP numerical threshhold for convergence of MAP, default: 10^-5
#' @return a list of the processed ForeGround, peaks, number of clusters, cluster assignments, landmarks, and peak p values
#' @import WeightedCluster 
#' @import Rsamtools
#' @import GenomicRanges
#' @export
#' @examples 
#' 
#' # bamfiles is a vector of bamfiles and peakfile is the location of the peaks
#' bamfiles = sapply(paste0("SRX860", c(379:474, 187:282), "Chr12.bam"), function(x) system.file("extdata", x, package = "scABC"))
#' peakfile = system.file("extdata", "chr12Peaks.bed", package = "scABC")
#' scABCcluster = scABC(bamfiles, peakfile, nClusters = 1:5) 
scABC <- function(bamfiles, peakfile, byReadGroup = FALSE,
                  PLOT = TRUE, QUIET = TRUE,
                  nClusters = c(), pValThresh = 2, 
                  nreadsThresh = 2, ncellsThresh = 10, readsFGthresh = NULL, 
                  lambda = 1, nTop = 2000, nPerm = 20, nTopGapStat = 5000, 
                  maxiter=1000, thresMLE=10^-3, thresMAP=10^-5){
  peaks = select_peaks(peakfile,thresh = pValThresh)
  if(!QUIET){cat("\nreading in foreground\n")}
    
  ForeGround = getCountsMatrix(bamfiles, peaks)
  ForeGroundFiltered = filterPeaks(ForeGround$ForeGroundMatrix, ForeGround$peaks,
                                    nreads_thresh = nreadsThresh, ncells_thresh = ncellsThresh)
  peaks = ForeGroundFiltered$peaks
  if(!QUIET){cat("\nreading in background\n")}
  BackGround = get_background(bamfiles, peaks)
  ForeGroundBackGroundFiltered = filterSamples(ForeGround = ForeGroundFiltered$ForeGroundMatrix, 
                                                BackGround = BackGround$BackGroundMatrix, 
                                                readsFGthresh = readsFGthresh)
  ForeGroundMatrix = ForeGroundBackGroundFiltered$ForeGroundMatrix
  BackGroundMatrix = ForeGroundBackGroundFiltered$BackGroundMatrix
  BackGroundMedian = apply(BackGroundMatrix, 2, median)
  if(!QUIET){cat("\ndone reading in data, beginning clustering\n")}
  if(length(nClusters) == 1){
    LandMarks = computeLandmarks(ForeGround = ForeGroundMatrix, 
                                 BackGround = BackGroundMatrix, 
                                 nCluster = nClusters, lambda = lambda, nTop = nTop)
  }
  else if(length(nClusters) > 1){
    GapStat = getGapStat(ForeGroundMatrix, BackGroundMedian, 
                         nClusters=nClusters, nPerm=nPerm, 
                         nTop = nTopGapStat, quiet=QUIET)
    if(PLOT){
      plotGapStat(GapStat, nClusters=nClusters, main = "Gap Stat")
    }
    nClusters = GapStat$nClusterOptimal
    LandMarks = computeLandmarks(ForeGround = ForeGroundMatrix, 
                                 BackGround = BackGroundMatrix,
                                 nCluster = GapStat$nClusterOptimal, 
                                 lambda = lambda, nTop = nTop)
  }
  else {
    nClusters = 1:10
    GapStat = getGapStat(ForeGroundMatrix, BackGroundMedian, 
                         nClusters = nClusters, nPerm=nPerm, 
                         nTop = nTopGapStat, quiet = QUIET)
    if(PLOT){
      plotGapStat(GapStat, nClusters=nClusters, main = "Gap Stat")
    }
    nClusters = GapStat$nClusterOptimal
    LandMarks = computeLandmarks(ForeGround = ForeGroundMatrix, 
                                 BackGround = BackGroundMatrix,
                                 nCluster = nClusters, 
                                 lambda = lambda, nTop = nTop)
  }
  LandMarkAssignments = assign2landmarks(ForeGround = ForeGroundMatrix, LandMarks)
  PeakPvals = getClusterSpecificPvalue(ForeGround = ForeGroundMatrix, cluster = LandMarkAssignments, 
                                       background_medians = BackGroundMedian, 
                                       maxiter=maxiter, thresMLE=thresMLE, thresMAP=10^-5, 
                                       quiet=QUIET)
  which_peaks = which(!is.na(PeakPvals$pvalue[,1]))
  
  if(!QUIET){
    cat("dim(peaks) = ", dim(peaks), "\n")
  }
  return(list(ForeGroundMatrix = ForeGroundMatrix[which_peaks, ], peaks = peaks[which_peaks, ],
              nClusters = nClusters, cluster_assignments = LandMarkAssignments,
               LandMarks = LandMarks[which_peaks, ], PeakPVals = PeakPvals$pvalue[which_peaks, ]))
}