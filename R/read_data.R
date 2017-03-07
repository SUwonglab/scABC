
select_peaks <- function(filename){
  column_names = c("chrom", "chromStart", "chromEnd", "name", "score", "strand",
                   "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
                   "blockStarts", "signalValue", "pValue", "qValue");
  gapped_peaks = read.table(file = filename, header = FALSE, sep = "\t",
                            stringsAsFactors = FALSE, col.names = column_names);
  wanted_peaks = which(gapped_peaks$pValue > 1); # pValue is -log10(p), p < 0.1 => pValue > 1
  return(gapped_peaks[wanted_peaks, ])
}
