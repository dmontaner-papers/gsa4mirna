##a030_format_res_uvgsa_paired.r
##2015-09-01 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script collects all GSA results generated under Godard's paradigm.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.2 (2015-08-14)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.0.0"
library (openxlsx); packageDescription ("openxlsx", fields = "Version") #"3.0.0"
#help (package = mdgsa)

try (source (".job.r")); try (.job)

options (width = 170)

corte <- 0.05

### DATA
setwd (file.path (.job$dir$proces))
load ("res0_godard_paired_nofilter.RData")
ls ()

tag <- "blca"
for (tag in names (res.godard.pair.nf)) {
    cat ("\n========== ", tag, " ==========\n")
    
    res <- res.godard.pair.nf[[tag]]
    ##
    res[,"pat"] <- uvPat (res, cutoff = corte)
    res[,"index"] <- pval2index (pval = res[,"pval"], sign = res[,"lor"])
    res[,"Name"] <- getGOnames (res, verbose = FALSE)
    ##
    orden <- order (res[,"index"], decreasing = TRUE)
    res <- res[orden,]
    ##
    res.godard.pair.nf[[tag]] <- res
}

## SAVE
save (res.godard.pair.nf, file = "res_godard_paired_nofilter.RData")

## xls
write.xlsx (res.godard.pair.nf, file = file.path (.job$dir$res, "res_godard_paired_nofilter.xlsx"))

################################################################################

### Matrix Format
gos <- sort (unique (unlist (lapply (res.godard.pair.nf, rownames))))
length (gos)

imat.godard.pair.nf <- matrix (NA, nrow = length (gos), ncol = length (res.godard.pair.nf))
rownames (imat.godard.pair.nf) <- gos
colnames (imat.godard.pair.nf) <- names (res.godard.pair.nf)
head (imat.godard.pair.nf)

for (tag in names (res.godard.pair.nf)) {
    cat ("\n========== ", tag, " ==========\n")
    imat.godard.pair.nf[,tag] <- res.godard.pair.nf[[tag]][gos, "index"]
}
head (imat.godard.pair.nf)

round (cor (imat.godard.pair.nf, use = "pairwise.complete.obs"), 2)

###SALVAMOS
save (imat.godard.pair.nf, file = "gsa_index_godard_paired_nofilter.RData")


###EXIT
warnings ()
sessionInfo ()
q ("no")
