sort_peaks <- function(peaks){
  return(peaks[order(peaks$chrom, peaks$chromStart), ])
}

#' select peaks
#' @param filename of a bed12+3 gapped peaks file obtain from peaking calling using MACS2
#' @return significant peaks obtained by filtering by p-value
#' @keywords peaks
#' @export
#' select_peaks()
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
get_counts_from_bam <- function(bamfile, peaks, ISPAIRED = FALSE){
  param = ScanBamParam(which = peaks, what = c("rname", "pos", "strand", "qwidth"))
  counts = countBam(bamfile, param = param, flag = scanBamFlag(isDuplicate = FALSE, isUnmappedQuery = FALSE))
  return(counts[,c("space", "start", "end", "records")])
}

