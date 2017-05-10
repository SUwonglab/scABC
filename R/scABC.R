#' scABC
#' @author Timothy Daley \email{tdaley@stanford.edu}, Zhixiang Lin \email{zl235@stanford.edu},
#'         Mahdi Zamanighomi \email{mzamani@stanford.edu}
#' @param bamfiles location of bam files
#' @param peakfile location of peak file
#' @param plot boolean variable to make plots
#' @param nClusters number of clusters, if not set manually then the number of clusters is determined by gap stat analysis
#' @import WeightedCluster 
#' @import Rsamtools
#' @import GenomicRanges
#' @output list of filtered data and peaks, cluster assigments, landmarks, and peak p-values
#' @examples bamfiles and peakfile contain the location of the bams and peaks (in bed format)
#'          scABCcluster = scABC(bamfiles, peakfile, nClusters = 1:5) chooses the optimal 
#'          cluster from 1 to 5 and returns the filtered data and peak along with cluster
#'          assignments, representative landmarks, and peak p-values
#' @export
scABC <- function(bamfiles, peakfile, PLOT = TRUE, QUIET = TRUE,
                  nClusters = c(), pValThresh = 2, 
                  nreadsThresh = 2, ncellsThresh = 10, readsFGthresh = NULL, 
                  lambda = 1, nTop = 2000, nPerm = 20, 
                  maxiter=1000, thresMLE=10^-3, thresMAP=10^-5){
  peaks = select_peaks(peakfile,thresh = pValThresh)
  if(!QUIET){cat("\nreading in foreground\n")}
    
  ForeGround = get_counts_matrix(bamfiles, peaks)
  ForeGroundFiltered = filter_peaks(ForeGround$ForeGroundMatrix, ForeGround$peaks,
                                    nreads_thresh = nreadsThresh, ncells_thresh = ncellsThresh)
  peaks = ForeGroundFiltered$peaks
  if(!QUIET){cat("\nreading in background\n")}
  BackGround = get_background(bamfiles, peaks)
  ForeGroundBackGroundFiltered = filter_samples(ForeGround = ForeGroundFiltered$ForeGroundMatrix, 
                                                   BackGround = BackGround$BackGroundMatrix, 
                                                   readsFGthresh = readsFGthresh)
  ForeGroundMatrix = ForeGroundBackGroundFiltered$ForeGroundMatrix
  BackGroundMatrix = ForeGroundBackGroundFiltered$BackGroundMatrix
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
    nClusters = 1:10
    GapStat = getGapStat(ForeGroundMatrix, BackGroundMedian, 
                         nClusters = nClusters, nPerm=nPerm, quiet = QUIET)
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
  PeakPvals = getClusterSpecificPvalue(ForeGround = ForeGroundMatrix, cluster = LandMarkAssignments, 
                                       background_medians = BackGroundMedian, landmark=LandMarks, 
                                       maxiter=maxiter, thresMLE=thresMLE, thresMAP=10^-5, 
                                       quiet=QUIET)
  which_peaks = which(!is.na(PeakPvals$pvalue[,1]))
  
  return(list(ForeGroundMatrix = ForeGroundMatrix[which_peaks, ], peaks = peaks[which_peaks, ],
              nClusters = nClusters, cluster_assignments = LandMarkAssignments,
               LandMarks = LandMarks[which_peaks, ], PeakPVals = PeakPvals$pvalue[which_peaks, ]))
}