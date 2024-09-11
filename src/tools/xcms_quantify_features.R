library('optparse')
library('xcms')
library('MsExperiment')
library('configr')

option_list <- list(
  make_option(c("--mzml"), type="character", 
              help="The .mzML file to search"),
  make_option(c("--output"), type="character", 
              help="The output file for identified features"),
  make_option(c("--xcms_params"), type="character", 
              help="The parameters for XCMS"),
  make_option(c("--peakmerge_params"), type="character", 
              help="The parameters for peak refinement"),
  make_option(c("--algorithm"), type="character", 
              help="what feature identification algorithm to use")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#read spectral data
mzml <- readMsExperiment(spectraFiles = opt$mzml)

#determine algorithm
if (opt$algorithm == 'xcms_cw') {
  alg_params <- CentWaveParam
} else if (opt$algorithm == 'xcms_cwip'){
  alg_params <- CentWavePredIsoParam
} else if (opt$algorithm == 'xcms_mf') {
  alg_params <- MatchedFilterParam
} else if (opt$algorithm == 'xcms_kalman') {
  alg_params <- MassifquantParam
}

#find peaks
xcms_params <- read.config(file = opt$xcms_params)
xcms_params <- do.call(alg_params, xcms_params)
peaks <- findChromPeaks(mzml, xcms_params)

#merge peaks
merge_params <- read.config(file = opt$peakmerge_params)
merge_params <- do.call(MergeNeighboringPeaksParam, merge_params)
peaks <- refineChromPeaks(peaks, merge_params)

#export data
write.table(chromPeaks(peaks), opt$output, sep = '\t', row.names = FALSE)