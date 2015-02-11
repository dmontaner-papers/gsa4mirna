##s010_sample_description_paired.r
##2014-07-07 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script describes the analyzed samples

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
#library (); packageDescription ("", fields = "Version") #

try (source (".job.r")); try (.job)

## DATA
setwd (file.path (.job$dir$proces))
load ("datos_paired.RData")
ls ()

tags <- names (datos)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

sapply (datos, names)
sapply (datos, function (x) dim (x$contmat))

sapply (datos, function (x) table (colnames (x$contmat) == rownames (x$sinfo)))

################################################################################


## Count cases and controls
tag <- "blca"

case.cont <- matrix (NA, nrow = length (tags), ncol = 3)
rownames (case.cont) <- tags
colnames (case.cont) <- c ("N.cases", "N.controls", "N.total")
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    sinfo   <- datos[[tag]]$sinfo
    case.cont[tag,] <- c (sum (sinfo$tumor == 1), sum (sinfo$tumor == 0), nrow (sinfo))
}
case.cont <- as.data.frame (case.cont)
case.cont[,"case.cont.ratio"] <- case.cont[,"N.cases"]  / case.cont[,"N.controls"]
case.cont


table (case.cont[,"N.cases"] + case.cont[,"N.controls"] == case.cont[,"N.total"])

#plot (case.cont[,"N.controls"] , case.cont[,"N.cases"])
#abline (0, 1, col = "red")

summary (case.cont)


## SAVE
save (list = c("case.cont"), file = file.path (.job$dir$proces, "report", "case_control_counts_paired.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
