##i010_explore_res_uvgsa_unpaired.r
##2014-11-26 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script explores the GSA ressults

###NOTES:
## Gene Set "evidence" is not related to Block Size

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
#library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.2"
library (ellipse); packageDescription ("ellipse", fields = "Version") #"0.3-8"
#library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mdgsa)
#help (package = mirbaseID)

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


## PLOTS
graphics.off ()
setwd (.job$dir$plots)



### Inhibition effects scatter
for (tag in mat.gsa.pair[["tags"]]) {
    print (tag)
    ##
    fichero.png <- paste0 ("inhibition_effect_paired_vs_unpaired_", tag, ".png")
    png (filename = fichero.png, width = por * 480 * 3, height = por * 480)
    par (mfcol = c(1, 3))
    for (onto in ontologias) {    
        plot (mat.gsa.unpa[[onto]][["index"]][,tag],
              mat.gsa.pair[[onto]][["index"]][,tag],
              xlab = "unpaired", ylab = "paired", main = paste0 (toupper (tag), " (", onto, ")"))
        abline (0, 1, col = "blue")
    }
    dev.off ()
}


##All correlation plots
por <- 2
ids <- mat.gsa.pair[["tags"]]
pair.unpa.gsa.res.cor <- list ()
for (onto in ontologias) {
    fichero.png <- paste0 ("inhibition_effect_cor_paired_vs_unpaired_", onto, ".png")
    png ("paired_cor_rindex0.png", width = por * 480, height = por * 480, pointsize = 12, bg = "white")
    mico <- cor (mat.gsa.unpa[[onto]][["index"]][,ids],
                 mat.gsa.pair[[onto]][["index"]])
    plotcorr (mico)
    abline (length (ids) + 1, -1)
    dev.off ()
    ##
    ##keep correlation
    pair.unpa.gsa.res.cor[[onto]] <- mico
}

sapply (pair.unpa.gsa.res.cor, diag)

###SAVE
save (list = "pair.unpa.gsa.res.cor", file.path (.job$dir$proces, "gsa_res_corelation_paired_unpaired.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
