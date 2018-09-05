#
# MStractor- A molecular feature extractor for mass spectrometry data
#
#######################################################
#                                                     #
#              *           |                  |       #
#           *              |         |        |       #
#         *  _____         |      |  |        |       #
#        |  /    |         |      |  |    |   |       #
#      __|_/_____|         || |  || ||  | |   ||      #
#     |___    /   \..... =========================    #
#    ]=/  \__|     |         / \           / \        #
#      \__/   \___/          \_/           \_/        #
#.||.|.....||*.....|.|||..:|....*.||:.|||.|*:.|||.||..#
#                                                     #
#######################################################

# Author: Jeremy Hack

# Intro -------------------------------------------------------------------

# This R-script is suitable for non-targeted profiling experiments.
# Molecular features are extracted from MS data files, retention time aligned,
# grouped and filtered to generate a data matrix ready for downstream
# statistical analysis.

# It is recommended to create a project directory ("./MyProject" or whatever
# you'd like to call the project) and copy this script there.

# MS data files must be in open file formats such as:
#   - CDF (tested, but doesn't support MSn for future releases)
#   - mzData (tested, but larger file size relative to mzXML)
#   - mzML (not-tested, but should work)
#   - mzXML (tested and preferable)
# Export proprietary formats directly from vendor software or use ProteoWizard
# http://proteowizard.sourceforge.net/index.shtml

# To automatically determine sample classes, group files in folders according to
# their class e.g. "./MyProject/MSfiles/Control", "./MyProject/MSfiles/Treated".

# Outputs intended for downstream analysis are stored in the working directory
# alongside this script. A subdirectory "QC" is created for various quality
# control outputs which are generally only required for troubleshooting and
# development. If "FastMode" is enabled some of these outputs won't be generated
# but the final (BP) matrix will be the same.

# The first time a user runs this script software dependencies are installed
# which will take some minutes.

# The main sections of the script are listed below. It is intended that user
# input is generally only required in the "Project Settings" and
# "Instrument Parameters" sections.
#
# "Project Settings" - Project specific input
# "Acquisition Parameters" - Parameters relating to the acquisition conditions
# "Prepare Environment" - Installed required software etc
# "Create dataset" - Extract data from MS files
# "RT Alignment" - Aligns the extracted molecular features
# "Reconstruct Spectra" - Construct deconvoluted spectra
# "Export Matrix" - Filter peaks and export basepeak matrix

# To clean things up and start afresh you may wish to:
# Clear the console pane. Objects in the environment pane unaffected.
cat("\014")

# Clear the environment pane (the objects). Console pane unaffected.
rm(list = ls(all.names = TRUE))

# Invoke garbage collection to reclaim RAM
gc()

# "Project Settings" - Project specific input -----------------------------

# Prior to running this script create a project directory and specify it's
# location. If the project directory is named "MStractorDemo" the MS files will
# be loaded automatically. Otherwise, place MS files in a sub-directory named
# "MSfiles" and group files according to their class into sub-directories of
# MSfiles.
# e.g.
# D:/MyProject
#       |--MSfiles
#       |     |--Control
#       |     |     |--MSfileC1.mzXML
#       |     |     |--MSfileC2.mzXML
#       |     |--Treatment
#       |     |     |--MSfileT1.mzXML
#       |     |     |--MSfileT2.mzXML
#       |-Demo.R (this script)

# Specify the path to the project folder. If left as "" a folder named
# "MStractorDemo" will be created in the user's home directory.
# e.g.  "C:/users/%UserName% on Windows
#       "/home/$HOME" on Linux

#pathToProject <- "D:/MStractorDemo"
pathToProject <- "D:/MyProject"
setwd(pathToProject)
getwd()

# Specify a "typical" MS file which can be considered a reference. Leave blank
# for MStractorDemo.
ref <- "./MSfiles/Control/MSfileC1.mzXML"

#QCdataAll <- TRUE # IF TRUE, EICs for every molecular feature will be plotted
CPUs <- "max"      # Concurrent CPU threads to run. "max" to autodetect
updatePackages <- "FALSE" # If TRUE, all packages will be automatically updated

# "Acquisition Parameters" - Parameters relating to the acquisition -------

# Chromatography parameters
# Tuned to faahKO data (Agilent 1100 400 bar HPLC)
rtStart <- 1    # Start region of interest (in seconds)
rtEnd <- "max"  # End region of interest (in seconds). "max" for RT full range
FWHM_min <- 10  # FWHM in seconds of narrowest peak (10)
FWHM_max <- 90  # FWHM in seconds of broadest peak (40)
rtDelta <- 3705-3673   # Max observed difference in retention time (s) across all
# samples (peak apex to peak apex delta).

# MS parameters
# Tuned for faahKO data (Agilent MSD SL ion trap)
mzPol <- "negative" # Set to "positive" or "negative" ion mode
mzStart <- 100  # Start of m/z region of interest
mzEnd <- 1650   # End of m/z region of interest
mzErrAbs <- 0.01 # Max m/z delta expected for the same feature across all samples
mzZmax <- 3     # Max charge state expected
EICsMax <- 30   # Max number of chrom. peaks expected for a single EIC
sens <- 1      # Factor (between 0 and 1) for peak extraction sensitivity
# Impacts peak picking thresholds, RAM & CPU utilisation.
# Start with ~0.5.
fileType <- "mzXML" # MS data file type e.g. "mzData", "mzML", "mzXML", "CDF"

# The script may be run through to the section "Peak Curation" near the end.
# Then follow the instructions there.

# End "Instrument parameters"

# ToDo: Check "Set Parameters"

# pathToProject exists?
# pathToProject/MSfiles exists?
# MSfiles     valid format? mzData, mzML, mzXML, netCDF

# QCdataAll  logical
# CPUs        "max" or CPUs > 0
# updatePackages logica
# rtStart     "min" OR numeric 1 < rtStart < rtEnd
# rtEnd       "max" OR numeric rtEnd > rtStart
# FWHM_min    numeric < FWHM_max & < 0
# FWHM_max    numeric > FWHM_min
# rtDelta     numeric < rtEnd-rtStart & < 0

# mzPol       "positive" OR "negative"
# mzStart     numeric < mzEnd & < 0
# mzEnd       numeric > mzStart
# mzErrAbs    numeric < mzStart-mzEnd & < 0
# mzZmax      numeric 0 < mzZmax < 10
# EICsMax     numeric 0 < EICsMax
# sens        numeric 0 < sens <= 1

# "Prepare Environment" - Installed required software etc -----------------
# Install software dependencies
libUser <- Sys.getenv("R_LIBS_USER")

update.packages(ask= updatePackages,lib= libUser)

# Install packages from CRAN
# Add required packages to this list
#pkgsCRAN <- c("snow","tools","Rcpp","rgl","parallel")
pkgsCRAN <- c("snow", "RANN")

for (i in pkgsCRAN) {
  if(i %in% rownames(installed.packages()) == FALSE)
  {install.packages(i, dependencies= TRUE, lib = libUser)}
}

# Install packages from Bioconductor
# Add required packages to this list
pkgsBioC <- c("xcms","CAMERA","modeest","mzR")

if ((basename(pathToProject) == "MStractorDemo")
    | basename(pathToProject) == ""){
  pkgsBioC <- c(pkgsCRAN,"faahKO")
}

source("http://bioconductor.org/biocLite.R")

# Ensure base packages are installed and auto-update others
if (updatePackages == TRUE){
  biocLite(suppressUpdates= FALSE, suppressAutoUpdate= FALSE, ask= FALSE, lib.loc= libUser)
} else {
  biocLite(suppressUpdates= TRUE, suppressAutoUpdate= FALSE, lib= libUser)
}

for (i in pkgsBioC) {
  if(i %in% rownames(installed.packages()) == FALSE){
    biocLite(i, suppressUpdates= TRUE, dependencies = TRUE, lib= libUser)}
}

# Load libraries
library(xcms)
library(RANN)
library(modeest)
library(CAMERA)
library(tools)
library(parallel)

# "Create dataset" - Extract data from MS files ---------------------------

# Record the time so that the total processing time can be determined later.
runtimeStart<- Sys.time()
runtimeStart

if (CPUs == "max"){
  CPUs <- detectCores(all.tests= TRUE, logical= TRUE)
}

# Set the projects working directory
if (pathToProject == ""){
  project <- "MStractorDemo"
  dir.create(paste(Sys.getenv("HOME"), project ,sep= "/"))
  pathToProject <- paste(Sys.getenv("HOME"), project ,sep= "/")
}
setwd(pathToProject)
cat("Working directory set to:\n",getwd())

# Load all the MS files
if (basename(pathToProject) == "MStractorDemo"){
  datapath <- file.path(find.package("faahKO"), "cdf")
  rawfiles <- dir(datapath, full.names=TRUE, pattern="\\.CDF", recursive=TRUE)
  ref <- paste(datapath,"/WT/wt19.CDF", sep= "")
} else if (fileType == "CDF"){
  datapath <- paste(pathToProject,"MSfiles" , sep="/")
  rawfiles <- dir(datapath, full.names=TRUE,
                  pattern="\\.CDF", recursive=TRUE)
} else if (fileType == "mzData"){
  datapath <- paste(pathToProject,"MSfiles" , sep="/")
  rawfiles <- dir(datapath, full.names=TRUE,
                  pattern="\\.mzData", recursive=TRUE)
} else if (fileType == "mzML"){
  datapath <- paste(pathToProject,"MSfiles" , sep="/")
  rawfiles <- dir(datapath, full.names=TRUE,
                  pattern="\\.mzML", recursive=TRUE)
} else if (fileType == "mzXML"){
  datapath <- paste(pathToProject,"MSfiles" , sep="/")
  rawfiles <- dir(datapath, full.names=TRUE,
                  pattern="\\.mzXML", recursive=TRUE)
} else (cat("No raw files found!\n\nSupported filetypes include:
            \n CDF \n mzData \n mzML \n mzXML", fill= FALSE)
)

cat("The following MS files have been loaded:\n", rawfiles, fill= TRUE, sep="")

# Load a reference file & define the scan range
# Set mz step size for seeking new EIC traces
profmethod <- "bin"
profStep <- mzErrAbs*4
refRaw <- xcmsRaw(ref, profstep= profStep, includeMSn= FALSE, mslevel= NULL,
                  scanrange= NULL)
refRaw
# Plot TIC of reference sample to file
dir.create("./QC")
graphics.off()
png("./QC/Ref_TIC.png", width = 1024, height = 768, units = "px")
plotTIC(refRaw, ident= FALSE, msident= FALSE)
dev.off()

# Determine scan range
scanStart <- head(which(refRaw@scantime > rtStart & refRaw@scantime < rtEnd),
                  n= 1)

if (identical(rtEnd, "max")) {
  scanEnd <- max(refRaw@scanindex)
  rtEnd <- refRaw@scantime[which(refRaw@scanindex == scanEnd)]
} else {
  scanEnd <- tail(which(refRaw@scantime > rtStart & refRaw@scantime < rtEnd),
                  n= 1)
}

# Set peak picking parameters

pwMin <- FWHM_min*1.3
pwMax <- pwMin*1.3
mzErrPpmMin <- mzErrAbs/2/mzEnd*1000000
mzErrPpmMax <- mzErrAbs/2/mzStart*1000000
mzErrPpmMean <- mean(c(mzErrPpmMin,mzErrPpmMax))
mzdiff <- mzErrAbs/5
intThresh <- as.integer(quantile(refRaw@env$intensity,1-sens)*10)
snThresh <- 30/sens
integ <- 1
fitGauss <- FALSE
sleep <- 0

# Plot picked peaks to file
#sleep <- 0.001
#png(file.path("./QC/Pks/%003d.png"), h=768, w=1024)
refPks <- findPeaks(refRaw, method= 'centWave', ppm= mzErrPpmMin*2,
                    peakwidth= c(pwMin, pwMax), snthresh= snThresh,
                    prefilter= c(5,intThresh), mzCenterFun= "mean",
                    integrate= integ, mzdiff= mzdiff, verbose.columns= TRUE,
                    fitgauss= fitGauss, noise= intThresh, sleep= sleep)
#dev.off()
#Pks <- refPks[,c("rt","mz","maxo","into","intb","sn","egauss")]
#write.table(Pks, file= "./Pks.tsv", sep="\t")

png("./QC/Ref_EICs_100.png", width = 1024, height = 768, units = "px")
plotPeaks(refRaw, refPks,  c(10,10), width = FWHM_min*10)
dev.off()

# Create xcms data set
# Using Centwave
# N.B. Whilst Centwave does not perform binning (it uses the raw data directly)
# downstream methods such as getEIC will use the step paramter

xset <- xcmsSet(rawfiles, profmethod= "bin", profparam= list(step= profStep),
                method='centWave', ppm= mzErrPpmMin*2,
                peakwidth= c(pwMin, pwMax), snthresh= snThresh,
                prefilter= c(3,intThresh), mzCenterFun= "wMean",
                integrate= integ, mzdiff= mzErrAbs/2, verbose.columns= TRUE,
                fitgauss= fitGauss, scanrange= c(scanStart,scanEnd),
                nSlaves= CPUs)

xset # To checkout slots use slotNames(xset); ls(xset@groups)

# Determine the lowest level of replication within a class? ---------------
classes <- list.dirs(datapath, full.names= FALSE, recursive= FALSE)
classSize <- vector(length= length(classes))
dirs <- list.dirs(datapath, full.names= TRUE, recursive= FALSE)
lst <- vector("list", length(dirs))
names(lst) <- dirs
for (i in 1:length(dirs)){
  classSize[i] <- length(dir(dirs[i]))
}
minClassSize <- min(classSize)
#medClassSize <- median(classSize)
#aveClassSize <- ave(classSize)
#modeClassSize <- mfv(classSize)

#bw <- FWHM_min/2
mzWid <- mzErrAbs
bw <- mzWid/2
minsamp <- ceiling(minClassSize)-1

# group-methods {xcms}
# A peak group is considered valid when it contains either, the minimum number
# of samples set by minSamp, OR at least one sample class satisfies minFrac.

# xset <- group(
#   xset, method= "density", bw= bw, minfrac= 0.3, minsamp= minsamp,
#   mzwid= mzWid, max= EICsMax, sleep= 0)
# xset

xset <- group(
  xset, method= "nearest", mzVsRTbalance= 10, mzCheck= mzErrAbs,
  rtCheck= rtDelta, kNN=10)
xset

# Function: TICsOverlaid --------------------------------------------------

getTIC <- function(file,rtcor=NULL) {
  object <- xcmsRaw(file)
  cbind(if (is.null(rtcor)) object@scantime else rtcor,
        rawEIC(object,mzrange=range(object@env$mz))$intensity)
}

# Overlay TIC from all files in current folder or from xcmsSet, create PNG

getTICs <- function(xcmsSet=NULL,files=NULL, pngName="TICs.png",
                    rt=c("raw","corrected")) {
  if (is.null(xcmsSet)) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""),
                         collapse = "|")
    if (is.null(files))
      files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
  }
  else {
    files <- filepaths(xcmsSet)
  }

  N <- length(files)
  TIC <- vector("list",N)

  for (i in 1:N) {
    cat(files[i],"\n")
    if (!is.null(xcmsSet) && rt == "corrected")
      rtcor <- xcmsSet@rt$corrected[[i]] else
        rtcor <- NULL
    TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
  }

  png(pngName,h=768, w=1024)
  #cols <- rainbow(N)
  cols <- as.integer(sampclass(xcmsSet)) + 1
  lty = 1:N
  pch = 1:N
  xlim = range(sapply(TIC, function(x) range(x[,1])))
  ylim = range(sapply(TIC, function(x) range(x[,2])))
  plot(0, 0, type="n", xlim = xlim, ylim = ylim,
       main = "Total Ion Chromatograms", xlab = "Retention Time (s)",
       ylab = "TIC")
  for (i in 1:N) {
    tic <- TIC[[i]]
    points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
  }
  legend("topright",paste(basename(files)), col = cols, lty = lty, pch = pch)
  dev.off()

  invisible(TIC)
}
# End Function: TICsOverlaid

# QC Plots ----------------------------------------------------------------
# Plot TICs as PNG
getTICs(xcmsSet= xset, pngName= "./QC/TICs_Raw.png", rt= "raw")

# if (QCdataAll == TRUE ){
#   QCplots <- c("mzdevhist", "rtdevhist", "mzdevmass", "mzdevtime",
#                "mzdevsample", "rtdevsample")
#
#   for (i in QCplots){
#     png(filename= paste("./QC/","QCplot_","align0_", i, ".png", sep= ""),
#         w=1280, h=1024)
#     plotQC(xset, what= i)
#     dev.off()
#   }

# Get EICs
xset_grps <- xset@groups
eicRange <- 2*mean(c(FWHM_min,FWHM_max))
eicsRaw <- getEIC(xset, mzrange=xset_grps, rtrange= eicRange ,
                  groupidx = 1:nrow(xset_grps), rt= "raw")

# Plot as individual PNGs
dir.create("./QC/EICs_Raw/")
do.call(file.remove,list(list.files("./QC/EICs_Raw", full.names= TRUE)))
graphics.off()
png(file.path("./QC/EICs_Raw/%003d.png"), h=768, w=1024)
plot(eicsRaw, xset)
dev.off()

# Export matrix prior to RT Alignment -------------------------------------
xsTable <- peakTable(xset, filebase= "./QC/xsTable", method= "medret",
                     value= "maxo")

write.table(xsTable[with(xsTable, order(rt, mz)), ],
            file= "./QC/xsTable.tsv", sep= "\t", quote= FALSE, col.names= NA)

# "RT Alignment" - Aligns the extracted molecular features ----------------
align_ref <- match(basename(ref),basename(rawfiles[]))

png(filename= "./QC/rtAlignLoess.png", w=1280, h=1024)
xsAlign <- retcor(xset, method= "loess", missing= 3, extra= 0, span= 0.3,
                  family= "gaussian", plottype= "deviation")
dev.off()
xsAlign

xsAlign <- group(
  xsAlign, method= "nearest", mzVsRTbalance= 10, mzCheck= mzErrAbs,
  rtCheck= rtDelta, kNN=10)
xsAlign

# Plot TICs
getTICs(xcmsSet= xsAlign, pngName= "./QC/TICs_Aligned.png", rt= "corrected")

# Retrieve missing data
xsFilled <- fillPeaks(xsAlign, method="chrom", nSlaves=CPUs)

# Plot EICs ---------------------------------------------------------------
xsFilledGrps <- xsFilled@groups
eicsFilled <- getEIC(xsFilled, mzrange= xsFilledGrps, rtrange= FWHM_max*2,
                     groupidx= 1:nrow(xsFilledGrps), rt= "corrected")

dir.create("./QC/EICs_Aligned/")
do.call(file.remove,list(list.files("./QC/EICs_Aligned", full.names= TRUE)))
graphics.off()
png(file.path("./QC/EICs_Aligned/%003d.png"), h=768, w=1024)

plot(eicsFilled, xsFilled)
dev.off()

# "Reconstruct Spectra" - Construct deconvoluted spectra ------------------
ppm <- mzErrPpmMean/2
mzabs <- mzErrAbs
minfrac <- minClassSize/length(rawfiles)
xs <- xsFilled # Choose an xcmxSet object

# Using CAMERA:
xs_an <- xsAnnotate(xs, polarity= mzPol, nSlaves= 1) # nSlaves>1 is broken
xs_an <- groupFWHM(xs_an, sigma= 6, perfwhm= 1, intval= "maxo")

xs_an <- findIsotopes(xs_an, maxcharge=mzZmax, maxiso=4, ppm= ppm,
                      mzabs= mzabs, intval="maxo", minfrac=minfrac,
                      filter= TRUE)
xs_an <- groupCorr(xs_an, cor_eic_th= 0.7, pval=0.1,
                   graphMethod="hcs", calcIso= TRUE, calcCiS = TRUE,
                   calcCaS= TRUE, cor_exp_th= 0.7)
# Not appropriate for EI spectra
#xs_an <- findAdducts(xs_an, ppm= ppm, mzabs= mzabs, multiplier= 3,
#                     polarity= mzPol)

PksAn <- getPeaklist(xs_an, intval="maxo")

#PksAn <- PksAn[order(as.numeric(PksAn$rt)),]

write.table(PksAn, file=paste("./QC/Pks_An", "tsv", sep="."), sep= "\t",
            col.names= NA, row.names= TRUE)

# "Export Matrix" - Filter peaks and export basepeak matrix ---------------

# 1. Filter lone groups.
# Count number of entities for each factor
pcgroups <- table(PksAn$pcgroup)
min_ions <- 2
PksAnFilt <- droplevels(PksAn[PksAn$pcgroup %in% names(pcgroups)
                              [pcgroups >= min_ions],,drop=FALSE])

# 2. Optionally, select signals that are part of an isotope group
#PksAnFilt <- PksAnFilt[PksAnFilt$isotopes != "", ]
#PksAnFilt <- PksAn[PksAn$isotopes != "", ]

# 3. Determine median signal intensity
sNames <- sampnames(xs)
RespMed <- apply(PksAnFilt[c(sNames)], 1, median, na.rm= TRUE)
PksAnFilt <-cbind(PksAnFilt, RespMed)

# 4. Keep only the most intense signal from each pcgroup
# Other values may be preferred
#RespMed <- PksAnFilt$RespMed
# May be required if factors aren't recognised
# PeakList$pcgroup <- as.factor(PeakList$pcgroup)
BasePks <- PksAnFilt[ RespMed == ave(RespMed, PksAnFilt$pcgroup,
                                     FUN= function(RespMed) max(RespMed))
                      , ]
# BasePks[order(as.numeric(BasePks$pcgroup)),]

write.table(BasePks[with(BasePks, order(rt, mz)), ],
            file=paste("Pks_BPs", "tsv", sep="."), sep= "\t",
            col.names= NA, row.names= TRUE)

# Fetch basePeak EICs ---------------------------------------------------------------
BP_EICs <- paste(sprintf("%03d",as.numeric(rownames(BasePks))),"png", sep= ".")

dir.create("./EICs_BasePeaks/")
do.call(file.remove,list(list.files("./EICs_BasePeaks/", full.names= TRUE)))
for (i in 1:length(BP_EICs)){
  sourceF <- paste("./QC/EICs_Aligned/", BP_EICs[i], sep= "")
  destF <- paste("./EICs_BasePeaks/", BP_EICs[i], sep= "")
  file.copy(from= sourceF, to= destF, overwrite= FALSE)
}

dir.create("./EICs_BasePeaks_Curated/")
do.call(file.remove,list(list.files("./EICs_BasePeaks_Curated/",
                                    full.names= TRUE)))
for (i in 1:length(BP_EICs)){
  sourceF <- paste("./QC/EICs_Aligned/", BP_EICs[i], sep= "")
  destF <- paste("./EICs_BasePeaks_Curated/", BP_EICs[i], sep= "")
  file.copy(from= sourceF, to= destF, overwrite= TRUE)
}

# Stop watch --------------------------------------------------------------
sessionInfo <- sessionInfo()

now <- Sys.time()
time <- unclass(as.POSIXlt(now))
timeStamp <-paste(time$year,time$mon, time$mday, time$hour, time$min, sep="_")
imageName <- paste(timeStamp, "RData", sep= ".")
save.image(imageName)
historyName <- paste(timeStamp,"RHistory", sep= ".")
savehistory(historyName)

runtimeEnd <- Sys.time()
runtime <- runtimeEnd - runtimeStart
runtime

gc()

# Peak Curation -----------------------------------------------------------
# The EIC's in ./EIC_BasePeaks should be reviewed. Create a folder named
# "./EICs_BasePeaks_Curated/" and only copy into it the .PNG files that
# demonstrate good chromatographic peak shape. Then run the following

files <- dir("./EICs_BasePeaks_Curated/")
pks <- basename(file_path_sans_ext(dir("./EICs_BasePeaks_Curated/")))
BasePksCur <- BasePks[sprintf("%03d",as.numeric(rownames(BasePks))) %in% pks,]
write.table(BasePksCur[with(BasePksCur, order(rt,mz)), ],
            file= paste("PksBPsCurated", "tsv", sep="."),
            sep= "\t", col.names= NA, row.names= TRUE)
