#' select peaks
#' @param filename of a bed12+3 gapped peaks file obtain from peaking calling using MACS2
#' @return significant peaks obtained by filtering by p-value
#' @keywords peaks
#' @export
#' select_peaks()
select_peaks <- function(filename){
  column_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
                   "blockStarts", "signalValue", "pValue", "qValue");
  gapped_peaks = read.table(file = filename, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, col.names = column_names);
  wanted_peaks = which(gapped_peaks$pValue > 1); # pValue is -log10(p), p < 0.1 => pValue > 1
  return(gapped_peaks[wanted_peaks, ])
}

