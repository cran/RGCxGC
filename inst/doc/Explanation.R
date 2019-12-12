## ----Fig 1, echo=FALSE, fig.cap= "Figure 1. Overview of general data processing pipeline in chromatography that is presented in the RGCxGC package.", out.width= "80%",message=FALSE, warning=FALSE, paged.print=FALSE, fig.align='center'----
knitr::include_graphics("images/dataProcessing.png")

## ----Fig 2, echo=FALSE, fig.align="center", fig.cap="Figure 2. The basic workflow of RGCxGC package. The functions for each step are in parenthesis. The double line between smooth andbaseline correction refers to the interchangeable pathway.", message=FALSE, warning=FALSE, out.width="80%", paged.print=FALSE----
knitr::include_graphics("images/basicWorkflow.jpg")

## ----Cran install, eval=FALSE, include=TRUE-----------------------------------
#  install.packages("RGCxGC")

## ----github install, eval=FALSE, include=TRUE---------------------------------
#  library(devtools)
#  install_github("DanielQuiroz97/RGCxGC")

## ----library call-------------------------------------------------------------
library(RGCxGC)

## ----Chrom import-------------------------------------------------------------
chrom_08 <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
MTBLS08 <- read_chrom(chrom_08, mod_time = 5)
slotNames(MTBLS08)

## ----Chrom plot, message=FALSE, warning=FALSE, paged.print=FALSE, out.width= "60%"----
# nlevels: Number of levels
# color.palette: The color palette to employ
library(colorRamps)
plot(MTBLS08, nlevels = 100, color.palette = matlab.like2)

## ----baseline correction, out.width= "60%"------------------------------------
MTBLS08_bc <- baseline_corr(MTBLS08)
plot(MTBLS08_bc, nlevels = 100, color.palette = matlab.like2)

## ----smoothing, out.width= "60%"----------------------------------------------
# Linear penalty with lambda equal to 10
MTBLS08_sm1 <- wsmooth(MTBLS08_bc, penalty = 1, lambda = 1e1)
plot(MTBLS08_sm1, nlevels = 100, color.palette = matlab.like,
     main = expression(paste(lambda, "= 10, penalty = 1")) )

## ----smoothing 2, out.width= "60%"--------------------------------------------
# Cuadratic penalty with lambda equal to 10
MTBLS08_sm2 <- wsmooth(MTBLS08_bc, penalty = 2, lambda = 1e1)
plot(MTBLS08_sm2, nlevels = 100, color.palette = matlab.like,
     main = expression(paste(lambda, "= 10, penalty = 2")) )

## ----2DCOW, out.width= "60%"--------------------------------------------------
# Reference chromatogram
chrom_09 <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
MTBLS09 <- read_chrom(chrom_09, mod_time = 5L)
# Baseline correction
MTBL09_bc <- baseline_corr(MTBLS09)
# Smoothing
MTBL09_sm2 <- wsmooth(MTBL09_bc, penalty = 2, lambda = 1e1)
# Alignment 
aligned <- twod_cow(sample_chrom = MTBLS08_sm2, ref_chrom =  MTBL09_sm2,
                   segments = c(10, 40), max_warp = c(1, 8))
plot(aligned, nlevels = 100, color.palette = matlab.like,
     main = "Aligned chromatogram")

## ----Batch import-------------------------------------------------------------
# Read Sample chromatogram
GB08_fl <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
MTBLS08 <- read_chrom(GB08_fl, mod_time = 5, verbose = F)
 
# Read reference chromatogram
GB09_fl <- system.file("extdata", "09GB.cdf", package = "RGCxGC")
MTBLS09 <- read_chrom(GB09_fl, mod_time = 5, verbose = F)

## ----Batch list---------------------------------------------------------------
batch_samples <- list(Chrom1 = MTBLS08, Chrom2 = MTBLS08, Chrom3 = MTBLS08)

## ----Batch alignment----------------------------------------------------------
batch_alignment <- batch_2DCOW(MTBLS09, batch_samples, c(10, 40), c(1, 10))
names(batch_alignment@Batch_2DCOW)

## ----Join---------------------------------------------------------------------
allChrom <- join_chromatograms(MTBLS09, MTBLS08)

## ----Join complex, eval=FALSE, include=TRUE-----------------------------------
#  join_complex <- join_chromatograms(batch_samp1, batch_samp2,#Two batch samples
#                                     Ref_chrom1 = reference_1,#User named argument
#                                     Ref_chrom2 = reference_2,#User named argument
#                                     groups = metadata_exp)   #Metadata

## ----metadata, message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE----
metadata <- data.frame(Names = c("08GB", "09GB", "14GB", "29GB", "34GB", "24GB"), stringsAsFactors = F)
metadata$Type = factor(c(rep("S. typhy Carriege", 3), rep("Control", 3)))

## ----print metadata, echo=FALSE-----------------------------------------------
knitr::kable(metadata)

## ----batch alignment 2--------------------------------------------------------
data(MTBLS579)

## ----MPCA---------------------------------------------------------------------
exp_MPCA <- m_prcomp(MTBLS579, center = T, scale = F)

## ----Scores, out.width= "60%"-------------------------------------------------
scores(exp_MPCA)

## ----negative loadings, out.width= "60%"--------------------------------------
# Negative loadings
plot_loading(exp_MPCA, type = "n", main = "Negative loadings",
             color.palette = matlab.like)

## ----positive loadings, out.width= "60%"--------------------------------------
# Positive loadings
plot_loading(exp_MPCA, type = "p", main = "Positive loadings",
             color.palette = matlab.like)

