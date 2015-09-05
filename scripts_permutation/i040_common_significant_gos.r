##i040_common_significatn_gos.r
##2014-12-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script explores the GSA significant terms

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
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

pat.unpa <- as.data.frame (rbind (mat.gsa.unpa[["bp"]][["pat"]],
                                  mat.gsa.unpa[["cc"]][["pat"]],
                                  mat.gsa.unpa[["mf"]][["pat"]]))
pat.unpa[,"SUM"] <- rowSums (abs (pat.unpa))
pat.unpa[,"onto"] <- rep (ontologias, sapply (ontologias, function (x) nrow (mat.gsa.unpa[[x]][["pat"]])))
pat.unpa[,"name"] <- getGOnames (pat.unpa)
pat.unpa[1:3,]
table (pat.unpa[,"SUM"])
table (pat.unpa[,"SUM"], pat.unpa[,"onto"])


pat.pair <- as.data.frame (rbind (mat.gsa.pair[["bp"]][["pat"]],
                                  mat.gsa.pair[["cc"]][["pat"]],
                                  mat.gsa.pair[["mf"]][["pat"]]))
pat.pair[,"SUM"] <- rowSums (abs (pat.pair))
pat.pair[,"onto"] <- rep (ontologias, sapply (ontologias, function (x) nrow (mat.gsa.pair[[x]][["pat"]])))
pat.pair[,"name"] <- getGOnames (pat.pair)
pat.pair[1:3,]
table (pat.pair[,"SUM"])
table (pat.pair[,"SUM"], pat.pair[,"onto"])

###

table (rownames (pat.unpa) == rownames (pat.pair)) ## OK same order
table (unpa = pat.unpa[,"SUM"], pair = pat.pair[,"SUM"])

################################################################################

file.unpa <- "common_enrichment_unpaired.xlsx"
unlink (file.unpa)
order.unpa <- order (pat.unpa[,"SUM"], pat.unpa[,"onto"], pat.unpa[,"name"], decreasing = TRUE)
order.unpa <- pat.unpa[order.unpa,]
order.unpa[1:3,]

suma.unpa <- sort (setdiff (unique (pat.unpa[,"SUM"]), 0), decreasing = TRUE)
suma.unpa
comen.unpa <- list ()  ## common enrichment
for (su in suma.unpa) {
    print (su)
    touse <- pat.unpa[,"SUM"] == su
    mat <- pat.unpa[touse,]
    comen.unpa[[su]] <- mat
    print (mat)
    write.xlsx2 (mat, file = file.unpa, sheetName = paste ("common", su),
                 col.names = TRUE, row.names = TRUE, append = TRUE)
}
sapply (comen.unpa, dim )

###


file.pair <- "common_enrichment_paired.xlsx"
unlink (file.pair)
order.pair <- order (pat.pair[,"SUM"], pat.pair[,"onto"], pat.pair[,"name"], decreasing = TRUE)
order.pair <- pat.pair[order.pair,]
order.pair[1:3,]

suma.pair <- sort (setdiff (unique (pat.pair[,"SUM"]), 0), decreasing = TRUE)
suma.pair
comen.pair <- list ()  ## common enrichment
for (su in suma.pair) {
    print (su)
    touse <- pat.pair[,"SUM"] == su
    mat <- pat.pair[touse,]
    comen.pair[[su]] <- mat
    print (mat)
    write.xlsx2 (mat, file = file.pair, sheetName = paste ("common", su),
                 col.names = TRUE, row.names = TRUE, append = TRUE)
}
sapply (comen.pair, dim )


### SAVE for the wiki
save (list = c("comen.pair", "comen.unpa"), file = file.path (.job$dir$proces, "common_enrichment.RData"))

## PLOTS
graphics.off ()
setwd (.job$dir$plots)




















###EXIT
warnings ()
sessionInfo ()
q ("no")
