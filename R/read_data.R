sort_peaks <- function(peaks){
  return(peaks[order(peaks$chrom, peaks$start), ])
}

#' select peaks
#' @param filename of a bed12+3 gapped peaks file obtain from peaking calling using MACS2
#' @return significant peaks obtained by filtering by p-value
#' @keywords peaks
#' @export select_peaks
select_peaks <- function(filename, thresh = 2){
  column_names = c("chrom", "start", "end", "name", "score", "strand",
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
                   "blockStarts", "signalValue", "pValue", "qValue");
  gapped_peaks = read.table(file = filename, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, col.names = column_names);
  wanted_peaks = which(gapped_peaks$pValue > thresh); # pValue is -log10(p), p < 0.1 => pValue > 2
  wanted_peaks = sort_peaks(gapped_peaks[wanted_peaks, ])
  return(wanted_peaks)
}

peaks2GRanges <- function(peaks, upstream = 0, downstream = 0){
  peaks.gr = with(peaks, GRanges(chrom, IRanges(sapply(start, function(x) max(0, x - upstream)), end + downstream), id = name, pVal = pValue))
}

# peaks should be in GenomicRanges 
get_counts_from_bam <- function(bamfile, peaks){
  param = ScanBamParam(which = peaks, what = c("rname", "pos", "strand", "qwidth"))
  counts = countBam(bamfile, param = param, flag = scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE))
  return(counts[,c("space", "start", "end", "file", "records")])
}

#' get_counts_matrix
#' @param bamfiles a vector of filenames of input bam files
#' @param peaks a bed15 format file returned from select_peaks
#' @return a matrix of chrom, start, end of peaks followed by counts of each bam file in bamfiles
#' @import Rsamtools 
#' @import GenomicRanges
#' @keywords peaks
#' @keywords counts
#' @export get_counts_matrix
get_counts_matrix <- function(bamfiles, peaks){
  peaks = peaks2GRanges(peaks)
  counts_list = lapply(bamfiles, function(x) get_counts_from_bam(x, peaks))
  sample_names = c(do.call(rbind, lapply(counts_list, function(x) head(toString(x$file[1])))))
  counts_mat = do.call(cbind, lapply(counts_list, function(x) x$records))
  colnames(counts_mat) = sample_names
  counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end, name = peaks$id, pValue = peaks$pVal)
  return(list(peaks = counts_info, ForeGroundMatrix = counts_mat))
}

#' compute_background
#' @param bamfiles a vector of filenames of input bamfiles
#' @param peaks a bed15 format file returned from select peaks
#' @param upstream number of bases upstream of peak to consider for computing background
#' @param downstream number of bases downstream of peak to consider for computing background
#' @return a matrix of chrom, start, end of peaks followed by background counts for each bam file in bamfiles
#' @import Rsamtools
#' @import GenomicRanges
#' @export get_background
get_background <- function(bamfiles, peaks, upstream = 500000, downstream = 500000){
  background_peaks = peaks2GRanges(peaks, upstream, downstream)
  counts_list = lapply(bamfiles, function(x) get_counts_from_bam(x, background_peaks))
  sample_names = c(do.call(rbind, lapply(counts_list, function(x) head(toString(x$file[1])))))
  background_counts = do.call(cbind, lapply(counts_list, function(x) x$records))
  colnames(background_counts) = sample_names
  counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end, name = peaks$name, pValue = peaks$pValue)
  return(list(peaks = counts_info, BackGroundMatrix = background_counts))
}

#' filter_samples
#' @param ForeGround matrix or data frame of ForeGround values
#' @param BackGround matrix or data frame of BackGround values
#' @param readsFGthresh threshold for the total reads per cell in ForeGround. Default is min(500, number of peaks/50)
#' @return filtered ForeGround and BackGround
#' @export filter_samples
filter_samples <- function(ForeGround, BackGround, readsFGthresh=NULL){
  stopifnot(dim(ForeGround) == dim(BackGround))
  if (is.null(readsFGthresh)){
    readsFGthresh <- min(500, nrow(ForeGround)/50)  
  } 
  which_samples_pass = which(colSums(ForeGround) > readsFGthresh) 
  return(list(ForeGroundMatrix = ForeGround[,which_samples_pass], 
              BackGroundMatrix = BackGround[,which_samples_pass]))
}

#' filter_peaks
#' @param ForeGround matrix or data frame of Foreground values
#' @param nreads_thresh threshold of the number of reads
#' @param ncells_thresh threshold of the number of cells
#' @return filtered ForeGround and peaks
#' @export filter_peaks
filter_peaks <- function(ForeGround, peaks, nreads_thresh = 2, ncells_thresh = 10){
  which_peaks_pass = which(rowSums(ForeGround[,4:dim(ForeGround)[2]] >= nreads_thresh) >= ncells_thresh)
  return(list(ForeGroundMatrix = ForeGround[which_peaks_pass,], peaks = peaks[which_peaks_pass, ]))
}

