sort_peaks <- function(peaks){
  return(peaks[order(peaks$chrom, peaks$start), ])
}

#' select peaks
#' 
#' read in peaks into a table and select only those with MACS2 pValue > 1 (p < 0.1)
#' @param filename of a bed12+3 gapped peaks file obtain from peaking calling using MACS2
#' @return significant peaks obtained by filtering by p-value
#' @keywords peaks
#' @export selectPeaks
selectPeaks <- function(filename, thresh = 2){
  peaks = read.table(file = filename, header = FALSE, sep = "\t",
                     stringsAsFactors = FALSE);
  if(dim(peaks)[2] == 15){
    # gapped peaks
    column_names = c("chrom", "start", "end", "name", "score", "strand",
                     "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes",
                     "blockStarts", "signalValue", "pValue", "qValue");
    colnames(peaks) = column_names
  }
  if(dim(peaks)[2] == 10){
    # narrow peaks
    column_names = c("chrom", "start", "end", "name", "score", "strand",
                     "foldChange", "pValue", "qValue", "summit2PeakDist")
    colnames(peaks) = column_names
  }
  
  wanted_peaks = which(peaks$pValue > thresh); # pValue is -log10(p), p < 0.1 => pValue > 2
  peaks = sort_peaks(peaks[wanted_peaks, ])
  return(peaks)
}

peaks2GRanges <- function(peaks, upstream = 0, downstream = 0){
  peaks.gr = with(peaks, GenomicRanges::GRanges(chrom, IRanges::IRanges(sapply(start, function(x) max(0, x - upstream)), end + downstream), id = name, pVal = pValue))
}

# peaks should be in GenomicRanges 
get_counts_from_bam <- function(bamfile, peaks){
  param = Rsamtools::ScanBamParam(which = peaks, what = c("rname", "pos", "strand", "qwidth"))
  counts = Rsamtools::countBam(bamfile, param = param, 
                               flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, 
                                                             isUnmappedQuery = FALSE))
  return(counts[,c("space", "start", "end", "file", "records")])
}


getTagCounts <- function(RGtag, bamfile, peaks){
  RGparam = ScanBamParam(which = peaks, 
                         what = c("rname", "pos", "strand", "qwidth"), 
                         tagFilter = list(RG= c("RG", RGtag)))
  counts = countBam(bamfile, 
                    param = RGparam, 
                    flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, 
                                                  isUnmappedQuery = FALSE))
  return(counts$records)
}

getCountsByReadGroup <- function(bamfile, peaks){
  scanned <- scanBam(bamfile, 
                     param = ScanBamParam(what = c("rname", "pos"),
                                                   tag = "RG"))[[1]]
  RGtags = unique(scanned$tag$RG)
  
  counts = do.call(cbind, lapply(RGtags, function(x) getTagCounts(x, bamfile, peaks)))
  colnames(counts) = RGtags
  
  return(counts)
}

getCountsByReadGroup2 <- function(bamfile, peaks, PAIRED = FALSE, VERBOSE = FALSE){
  scanned <- scanBam(bamfile, 
                     param = ScanBamParam(what = c("rname", "pos", "strand", "qwidth"),
                                          tag = "RG"))[[1]]
  RGtags = unique(scanned$tag$RG)
  
  counts_mat = matrix(nrow = length(peaks), ncol = length(RGtags))
  for(i in 1:length(RGtags)){
    tag = RGtags[i]
    if(VERBOSE){
      message("Processing tag ", tag)
    }
    match_RG <- which(scanned$tag$RG == tag)
    # convert bamfiles to Genomic Ranges
    bam.gr = GRanges(seqnames = scanned$rname[match_RG], 
                     IRanges(start = sapply(match_RG, function(i) ifelse(scanned$strand[i] == "-", 
                                                                         scanned$pos[i] + scanned$qwidth[i] - 1, 
                                                                         scanned$pos[i])),
                             width = scanned$qwidth[match_RG]))
    counts_mat[,i] = countOverlaps(peaks, bam.gr, type = "any", ignore.strand = TRUE)
  }
 # counts = do.call(cbind, lapply(RGtags, function(x) getTagCounts(x, bamfile, peaks)))
  colnames(counts_mat) = RGtags
  
  return(counts_mat)
}


#' get counts matrix
#' 
#' Obtain counts matrix from bamfiles and peaks
#' 
#' @param bamfiles a vector of filenames of input bam files
#' @param peaks a bed15 format file returned from select_peaks
#' @param byReadGroup boolean variable indicating whether or not individual experiments are separated by read group
#' @return a matrix of chrom, start, end of peaks followed by counts of each bam file in bamfiles
#' @import Rsamtools 
#' @import GenomicRanges
#' @keywords peaks
#' @keywords counts
#' @export getCountsMatrix
getCountsMatrix <- function(bamfiles, peaks, byReadGroup = FALSE){
  peaks.gr = peaks2GRanges(peaks)
  if(byReadGroup){
    counts_mat = getCountsByReadGroup(bamfile, peaks.gr);
    rownames(counts_mat) = peaks$name
  }
  else{
    counts_list = lapply(bamfiles, function(x) get_counts_from_bam(x, peaks.gr)) # this is the roadblock
    sample_names = c(do.call(rbind, lapply(counts_list, function(x) head(toString(x$file[1])))))
    counts_mat = do.call(cbind, lapply(counts_list, function(x) x$records))
    colnames(counts_mat) = sample_names
    rownames(counts_mat) = peaks$name
    #counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end, name = peaks$id, pValue = peaks$pVal)
  }
  peaks = peaks[,c("chrom", "start", "end", "name", "pValue")]
  return(list(peaks = peaks, ForeGroundMatrix = counts_mat))
}

#' compute background counts matrix
#' 
#' @param bamfiles a vector of filenames of input bamfiles
#' @param peaks a bed15 format file returned from select peaks
#' @param upstream number of bases upstream of peak to consider for computing background
#' @param downstream number of bases downstream of peak to consider for computing background
#' @return a matrix of chrom, start, end of peaks followed by background counts for each bam file in bamfiles
#' @import Rsamtools
#' @import GenomicRanges
#' @export getBackground
getBackground <- function(bamfiles, peaks, upstream = 500000, 
                          downstream = 500000, byReadGroup = FALSE){
  background_peaks.gr = peaks2GRanges(peaks, upstream, downstream)
  if(byReadGroup){
    counts_mat = getCountsByReadGroup(bamfile, background_peaks.gr);
    rownames(counts_mat) = peaks$name
  }
  else{
    counts_list = lapply(bamfiles, function(x) get_counts_from_bam(x, background_peaks.gr))
    sample_names = c(do.call(rbind, lapply(counts_list, function(x) head(toString(x$file[1])))))
    counts_mat = do.call(cbind, lapply(counts_list, function(x) x$records))
    colnames(counts_mat) = sample_names
    rownames(counts_mat) = peaks$name
    #counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end, name = peaks$id, pValue = peaks$pVal)
  }
  peaks = peaks[,c("chrom", "start", "end", "name", "pValue")]
  return(list(peaks = peaks, BackGroundMatrix = counts_mat))
}

#' get counts matrix
#' 
#' Obtain counts matrix from bamfiles and peaks
#' 
#' @param bamfiles a vector of filenames of input bam files
#' @param peaks a bed15 format file returned from select_peaks
#' @param paired a boolean variable indicating whether the experiment was paired end or not
#' @param byReadGroup boolean variable indicating whether or not individual experiments are separated by read group
#' @param VERBOSE a boolean variable to indicate whether to print update messages
#' @return a matrix of chrom, start, end of peaks followed by counts of each bam file in bamfiles
#' @import Rsamtools 
#' @import GenomicRanges
#' @import IRanges
#' @keywords peaks
#' @keywords counts
#' @export getCountsMatrix2
getCountsMatrix2 <- function(bamfiles, peaks, PAIRED = FALSE, 
                             byReadGroup = FALSE, VERBOSE = FALSE){
  peaks.gr = peaks2GRanges(peaks)
  if(byReadGroup){
    counts_mat = getCountsByReadGroup2(bamfile, peaks.gr);
    rownames(counts_mat) = peaks$name
  }
  else{
    nCells = length(bamfiles)
    counts_mat = matrix(nrow = length(peaks.gr), ncol = nCells)
    for(i in 1:nCells){
      if(VERBOSE){
        message("Processing file ", bamfiles[i])
      }
      # convert bamfiles to Genomic Ranges
      bam.gr = bam2gr(bamfiles[i], PAIRED = PAIRED)
      counts_mat[,i] = countOverlaps(peaks.gr, bam.gr, type = "any", ignore.strand = TRUE)
    }
    colnames(counts_mat) = bamfiles
    rownames(counts_mat) = peaks$name
    #counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end, name = peaks$id, pValue = peaks$pVal)
  }
  peaks = peaks[,c("chrom", "start", "end", "name", "pValue")]
  return(list(peaks = peaks, ForeGroundMatrix = counts_mat))
}

#' compute background counts matrix
#' 
#' @param bamfiles a vector of filenames of input bamfiles
#' @param peaks a bed15 format file returned from select peaks
#' @param upstream number of bases upstream of peak to consider for computing background
#' @param downstream number of bases downstream of peak to consider for computing background
#' @return a matrix of chrom, start, end of peaks followed by background counts for each bam file in bamfiles
#' @import Rsamtools
#' @import GenomicRanges
#' @import IRanges
#' @export getBackground2
getBackground2 <- function(bamfiles, peaks, upstream = 500000, 
                          downstream = 500000, byReadGroup = FALSE, 
                          VERBOSE = FALSE, PAIRED = FALSE){
  nCells = length(bamfiles)
  background_peaks.gr = peaks2GRanges(peaks, upstream, downstream)
  if(byReadGroup){
    counts_mat = getCountsByReadGroup2(bamfile, background_peaks.gr);
    rownames(counts_mat) = peaks$name
  }
  else{
    counts_mat = matrix(nrow = length(background_peaks.gr), ncol = nCells)
    for(i in 1:nCells){
      if(VERBOSE){
        message("Processing file ", bamfiles[i])
      }
      # convert bamfiles to Genomic Ranges
      bam.gr = bam2gr(bamfiles[i], PAIRED = PAIRED)
      counts_mat[,i] = countOverlaps(background_peaks.gr, bam.gr, type = "any", ignore.strand = TRUE)
    }
    colnames(counts_mat) = bamfiles
    rownames(counts_mat) = peaks$name
    #counts_info = data.frame(chrom = counts_list[[1]]$space, start = counts_list[[1]]$start, end = counts_list[[1]]$end, name = peaks$id, pValue = peaks$pVal)
  }
  peaks = peaks[,c("chrom", "start", "end", "name", "pValue")]
  return(list(peaks = peaks, BackGroundMatrix = counts_mat))
}


# from chromVAR
#' @importFrom IRanges PartitioningByEnd
#' @importFrom BiocGenerics relist
bam2gr <- function(bamfile, PAIRED = FALSE) {
  if (PAIRED) {
    scanned = scanBam(bamfile, param = ScanBamParam(flag =scanBamFlag(isMinusStrand = FALSE, 
                                                                      isProperPair = TRUE),
                                                    what = c("rname", "pos", "isize")))[[1]]
    out = GRanges(seqnames = scanned$rname, IRanges(start = scanned$pos, width = scanned$isize))
  } else {
    scanned = scanBam(bamfile, param = ScanBamParam(what = c("rname", "pos", 
                                                     "strand", 
                                                     "qwidth")))[[1]]
    out = GRanges(seqnames = scanned$rname, 
                  IRanges(start = ifelse(scanned$strand == "-", scanned$pos + scanned$qwidth - 1, scanned$pos),
                          width = scanned$qwidth))
  }
  return(out)
}


#' filter samples by comparing foreground against background
#' 
#' @param ForeGround matrix or data frame of ForeGround values
#' @param BackGround matrix or data frame of BackGround values
#' @param readsFGthresh threshold for the total reads per cell in ForeGround. Default is min(500, number of peaks/50)
#' @return filtered ForeGround and BackGround
#' @export filterSamples
filterSamples <- function(ForeGround, BackGround, readsFGthresh=NULL){
  stopifnot(dim(ForeGround) == dim(BackGround))
  if (is.null(readsFGthresh)){
    readsFGthresh <- min(500, nrow(ForeGround)/50)  
  } 
  which_samples_pass = which(colSums(ForeGround) > readsFGthresh) 
  return(list(ForeGroundMatrix = ForeGround[,which_samples_pass], 
              BackGroundMatrix = BackGround[,which_samples_pass]))
}

#' filter peaks
#' 
#' @param ForeGround matrix or data frame of Foreground values
#' @param peaks a bed format file of peaks
#' @param nreads_thresh threshold of the number of reads
#' @param ncells_thresh threshold of the number of cells
#' @return filtered ForeGround and peaks
#' @export filterPeaks
filterPeaks <- function(ForeGround, peaks, nreads_thresh = 1, ncells_thresh = 10){
  which_peaks_pass = which(rowSums(ForeGround[,4:dim(ForeGround)[2]] >= nreads_thresh) >= ncells_thresh)
  return(list(ForeGroundMatrix = ForeGround[which_peaks_pass,], peaks = peaks[which_peaks_pass, ]))
}

