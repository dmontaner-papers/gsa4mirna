##a020_gsa_at_mirna_level_paired.r
##2015-01-09 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script performs a GSA analysis at miRNA level following Godard's paradigm.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.2.2 (2015-08-14)"
##library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.0.0"
#help (package = mdgsa)

try (source (".job.r")); try (.job)

################################################################################

## DATA
#load (file.path (.job$dir$data, "data_processed", "res_dif_exp_paired.RData"))
load (file.path (.job$dir$data, "data_processed", "rindex0_paired.RData"))

## ANNOTATION
setwd (file.path (.job$dir$proces))
load ("annotation_paired.RData")
ls ()

## annotated miRNAs
annot.mirnas <- unique (unlist (annot))
length (annot.mirnas)

tag <- "blca"
res.godard.pair <- list ()
for (tag in names (rindex0)) {
    cat ("\n========== ", tag, " ==========\n")

    ## keep just annotated miRNAs as in Godard 2015
    rindex <- rindex0[[tag]]
    rindex <- rindex[names (rindex) %in% annot.mirnas]  ## FILTERING: using just the annotated miRNA universe
    
    ## transform (normalize) index
    rindex <- indexTransform (index = rindex, method = "normalize")
    ## boxplot (rindex)
    ## boxplot (rindex0[[tag]])
    ## plot (rindex0[[tag]], rindex)
    
    ## filter annotation
    annotF <- annotFilter (annot = annot, index = rindex, minBlockSize = 10, maxBlockSize = 300)
    if (.job$testmode) annotF <- annotF[1:3]
    
    ## uvGsa
    res <- uvGsa (index = rindex, annot = annotF)
    print (dim (uvSignif (res)))

    ## save
    res.godard.pair[[tag]] <- res
}

### SAVE
save (res.godard.pair, file = "res0_godard_paired.RData")


###EXIT
warnings ()
sessionInfo ()
q ("no")
