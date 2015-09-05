##c010_compare_filtered_unfiltered.r
##2015-09-01 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script compares the filtered (Godard's method) with the unfiltered results of the GSA analysis at miRNA level.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.2 (2015-08-14)"
#library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.0.0"
#library (openxlsx); packageDescription ("openxlsx", fields = "Version") #"3.0.0"
#help (package = mdgsa)

try (source (".job.r")); try (.job)

options (width = 170)

#corte <- 0.05

### DATA
setwd (file.path (.job$dir$proces))
load ("gsa_index_godard_unpaired_nofilter.RData")
load ("gsa_index_godard_unpaired.RData")
load ("gsa_index_godard_paired_nofilter.RData")
load ("gsa_index_godard_paired.RData")

ls ()

dim (imat.godard.pair)
dim (imat.godard.pair.nf)

dim (imat.godard.unpa)
dim (imat.godard.unpa.nf)

table (rownames (imat.godard.pair) == rownames (imat.godard.pair.nf))
table (rownames (imat.godard.unpa) == rownames (imat.godard.unpa.nf))

comunes <- intersect (rownames (imat.godard.pair), rownames (imat.godard.unpa))
length (comunes)


## EXPLORE CORRELATION
setwd (.job$dir$plots)
graphics.off ()
tags <- colnames (imat.godard.unpa)
tags

tag <- "blca"
for (tag in tags) {
    print (tag)
    png (filename = paste0 ("filtered_vs_unfiltered_", tag, ".png"))
    ##x11 ()
    par (mfrow = c (2,2))
    ##
    try ({plot (imat.godard.unpa[,tag], imat.godard.unpa.nf[,tag], main = tag); abline (0, 1, col = "red")})
    try ({plot (imat.godard.pair[,tag], imat.godard.pair.nf[,tag], main = tag); abline (0, 1, col = "red")})
    ##
    try ({plot (imat.godard.unpa   [comunes, tag], imat.godard.pair   [comunes, tag], main = tag); abline (0, 1, col = "red")})
    try ({plot (imat.godard.unpa.nf[comunes, tag], imat.godard.pair.nf[comunes, tag], main = tag); abline (0, 1, col = "red")})
    ##
    dev.off ()
}


###EXIT
warnings ()
sessionInfo ()
q ("no")
