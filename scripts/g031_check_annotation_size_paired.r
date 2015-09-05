##g031_get_annotation_size_paired.r
##2015-08-28 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script counts the number of go terms analyzed in the paired analysis

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.3"

try (source (".job.r")); try (.job)

## ANNOTATION
load (file.path (.job$dir$proces, "annot_for_paired_data.RData"))
ls ()

ontologias <- names (annot)
ontologias

## INDEX
load (file.path (.job$dir$proces, "rindex_paired.RData"))
ls ()

tags <- names (rindex)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

## GSA
lons <- list ()
for (tag in tags) {
    for (onto in ontologias) {
        cat ("\n=============== ", TAGS[tag], ":", onto, " ===============\n")
        anotacion <- annotFilter (annot[[onto]], index = rindex[[tag]], minBlockSize = 10, maxBlockSize = 300)
        lons[[onto]][tag] <- length (anotacion)
    }
}

as.data.frame (lons)
rowSums (as.data.frame (lons))
sort (unique (rowSums (as.data.frame (lons))))


###EXIT
warnings ()
sessionInfo ()
q ("no")
