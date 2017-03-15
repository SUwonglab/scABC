sort_peaks <- function(peaks){
  return(peaks[order(peaks$chrom, peaks$chromStart), ])
}

#' select peaks
#' @param filename of a bed12+3 gapped peaks file obtain from peaking calling using MACS2
#' @return significant peaks obtained by filtering by p-value
#' @keywords peaks
#' @export select_peaks
select_peaks <- function(filename, thresh = 1){
  column_names = c("chrom", "start", "end", "name", "score", "strand",
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
                   "blockStarts", "signalValue", "pValue", "qValue");
  gapped_peaks = read.table(file = filename, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, col.names = column_names);
  wanted_peaks = which(gapped_peaks$pValue > thresh); # pValue is -log10(p), p < 0.1 => pValue > 1
  wanted_peaks = sort_peaks(gapped_peaks[wanted_peaks, ])
  return(wanted_peaks)
}

peaks2GRanges <- function(peaks){
  peaks.gr = with(peaks, GRanges(chrom, IRanges(start, end), id = name, pVal = pValue))
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
  counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end)
  return(data.frame(counts_info, counts_mat))
}

