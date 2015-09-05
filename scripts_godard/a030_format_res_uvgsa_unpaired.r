##a030_format_res_uvgsa_unpaired.r
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
load ("res0_godard_unpaired.RData")
ls ()

tag <- "blca"
for (tag in names (res.godard.unpa)) {
    cat ("\n========== ", tag, " ==========\n")
    
    res <- res.godard.unpa[[tag]]
    ##
    res[,"pat"] <- uvPat (res, cutoff = corte)
    res[,"index"] <- pval2index (pval = res[,"pval"], sign = res[,"lor"])
    res[,"Name"] <- getGOnames (res, verbose = FALSE)
    ##
    orden <- order (res[,"index"], decreasing = TRUE)
    res <- res[orden,]
    ##
    res.godard.unpa[[tag]] <- res
}

## SAVE
save (res.godard.unpa, file = "res_godard_unpaired.RData")

## xls
write.xlsx (res.godard.unpa, file = file.path (.job$dir$res, "res_godard_unpaired.xlsx"))

################################################################################

### Matrix Format
gos <- sort (unique (unlist (lapply (res.godard.unpa, rownames))))
length (gos)

imat.godard.unpa <- matrix (NA, nrow = length (gos), ncol = length (res.godard.unpa))
rownames (imat.godard.unpa) <- gos
colnames (imat.godard.unpa) <- names (res.godard.unpa)
head (imat.godard.unpa)

for (tag in names (res.godard.unpa)) {
    cat ("\n========== ", tag, " ==========\n")
    imat.godard.unpa[,tag] <- res.godard.unpa[[tag]][gos, "index"]
}
head (imat.godard.unpa)

round (cor (imat.godard.unpa, use = "pairwise.complete.obs"), 2)

###SALVAMOS
save (imat.godard.unpa, file = "gsa_index_godard_unpaired.RData")


###EXIT
warnings ()
sessionInfo ()
q ("no")
