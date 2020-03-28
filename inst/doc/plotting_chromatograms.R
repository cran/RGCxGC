## ----slots, message=FALSE, warning=FALSE, paged.print=FALSE-------------------
library(RGCxGC)
chrom_08 <- system.file("extdata", "08GB.cdf", package = "RGCxGC")
MTBLS08 <- read_chrom(chrom_08, mod_time = 5, verbose = F)
slotNames(MTBLS08)

## ----default, fig.align='center'----------------------------------------------
plot(MTBLS08)

## ----color_palette, fig.align='center'----------------------------------------
library(colorRamps)
plot(MTBLS08, type = "f", color.palette = matlab.like,
     main =  "matlab.like")

plot(MTBLS08, type = "f", color.palette = matlab.like2,
     main =  "matlab.like2")

## ----message=FALSE, warning=FALSE, paged.print=FALSE, fig.align='center'------
myl_d5 <- system.file("extdata", "mylbd5.CDF", package = "RGCxGC")
myl <- read_chrom(myl_d5, mod_time = 5, sam_rate = 25, verbose = F)
plot(myl, color.palette = matlab.like2 )

## ----Fig 1, echo=FALSE, out.width= "60%",message=FALSE, warning=FALSE, paged.print=FALSE, fig.align='center'----
knitr::include_graphics("images/may_chrom.jpg")

## -----------------------------------------------------------------------------
myl_deph <- dephase_chrom(myl, rel_dephase = 65)
plot(myl_deph, color.palette = matlab.like2 )

## -----------------------------------------------------------------------------
plot(myl_deph, type = "c",col = matlab.like2(10) )

## -----------------------------------------------------------------------------
my_palette <- colorRampPalette(rev(c("red","yellow","springgreen",
                                     "blue", "white")))
plot(myl_deph, type = "c", col = my_palette(10) )

## -----------------------------------------------------------------------------
plot(myl_deph, type = "c", col = my_palette(35), nlevels = 100 )

