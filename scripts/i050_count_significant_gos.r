##i050_count_significant_gos.r
##2014-12-04 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script counts the number of significant GO terms in all cancers

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.2 (2014-10-31)"
#library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.7"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.2"
#library (ellipse); packageDescription ("ellipse", fields = "Version") #"0.3-8"

try (source (".job.r")); try (.job)

options (width = 170)
#options (width = 250)

corte <- 0.05
por <- 1
ontologias <- c ("bp", "cc", "mf")

## DATOS
load (file.path (.job$dir$proces, "res_gsa_paired.RData"))
load (file.path (.job$dir$proces, "res_gsa_unpaired.RData"))
#load (file.path (.job$dir$proces, "res_dif_exp_unpaired.RData"))
ls ()
names (mat.gsa.unpa)
names (mat.gsa.pair)

mat.gsa.unpa[["desc"]]

setwd (.job$dir$res)

################################################################################

##function
conteo <- function (x) {
    up = sum (x == 1)
    do = sum (x ==-1)
    no = sum (x == 0)
    return (c(deregulated = do, noDif = no, inhibited = up))
}


## unpa
gsa.res.count.bp.unpa <- t (apply (mat.gsa.unpa[["bp"]][["pat"]], 2, conteo))
gsa.res.count.cc.unpa <- t (apply (mat.gsa.unpa[["cc"]][["pat"]], 2, conteo))
gsa.res.count.mf.unpa <- t (apply (mat.gsa.unpa[["mf"]][["pat"]], 2, conteo))

gsa.res.count.bp.unpa
gsa.res.count.cc.unpa
gsa.res.count.mf.unpa

total.unpa <- as.data.frame (gsa.res.count.bp.unpa + gsa.res.count.cc.unpa + gsa.res.count.mf.unpa)
total.unpa

## pair
gsa.res.count.bp.pair <- t (apply (mat.gsa.pair[["bp"]][["pat"]], 2, conteo))
gsa.res.count.cc.pair <- t (apply (mat.gsa.pair[["cc"]][["pat"]], 2, conteo))
gsa.res.count.mf.pair <- t (apply (mat.gsa.pair[["mf"]][["pat"]], 2, conteo))

gsa.res.count.bp.pair
gsa.res.count.cc.pair
gsa.res.count.mf.pair

total.pair <- as.data.frame (gsa.res.count.bp.pair + gsa.res.count.cc.pair + gsa.res.count.mf.pair)
total.pair

###

touse <- rownames (total.unpa) %in% rownames (total.pair)
total.unpa[touse,]
total.unpa[!touse,]

cbind (total.unpa, total.pair[rownames (total.unpa),])


### SAVE for the wiki
save (list = c (
          "gsa.res.count.bp.unpa", "gsa.res.count.cc.unpa", "gsa.res.count.mf.unpa", "total.unpa",
          "gsa.res.count.bp.pair", "gsa.res.count.cc.pair", "gsa.res.count.mf.pair", "total.pair"),
      file = file.path (.job$dir$proces, "res_gsa_conteos.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
