##i010_explore_res_uvgsa_paired.r
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

ontologias <- c ("bp", "cc", "mf")

## DATOS
load (file.path (.job$dir$proces, "res_gsa_paired.RData"))
ls ()
names (mat.gsa.pair)

por <- 1

graphics.off ()
setwd (.job$dir$plots)

## Size dependencies
for (tag in mat.gsa.pair[["tags"]]) {
    print (tag)
    fichero.png <- paste0 ("paired_size_effect_", tag, ".png")
    png (filename = fichero.png, width = por * 480 * 3, height = por * 480 *2)
    par (mfcol = c(2, 3))
    for (onto in ontologias) {
        ## PLOT
        plot (log (mat.gsa.pair[[onto]][["N"]]), mat.gsa.pair[[onto]][["pval"]][,tag],  xlab = "block size (log scale)", ylab = "p-value",              main = paste0 (toupper (tag), " (", onto, ")"), cex.lab = 2, cex.main = 2)
        plot (log (mat.gsa.pair[[onto]][["N"]]), mat.gsa.pair[[onto]][["index"]][,tag], xlab = "block size (log scale)", ylab = "-log pval * sign lor", main = paste0 (toupper (tag), " (", onto, ")"), cex.lab = 2, cex.main = 2)
        abline (h = 0, col = "blue")
        ## PLOT
        ##plot (mat.gsa.pair[[onto]][["N"]], mat.gsa.pair[[onto]][["lor"]][,tag], xlab = "block size", ylab = "lor", main = paste (onto, tag, sep = " : "))
        ##plot (log (mat.gsa.pair[[onto]][["N"]]), mat.gsa.pair[[onto]][["lor"]][,tag], xlab = "block size", ylab = "lor", main = paste (onto, tag, sep = " : "))
        ##plot (log (mat.gsa.pair[[onto]][["N"]]), mat.gsa.pair[[onto]][["pval"]][,tag], xlab = "block size", ylab = "pval", main = paste (onto, tag, sep = " : "))
        ##plot (log (mat.gsa.pair[[onto]][["N"]]), - log (mat.gsa.pair[[onto]][["pval"]][,tag]), xlab = "block size", ylab = "- log pval", main = paste (onto, tag, sep = " : "))
        ##plot (mat.gsa.pair[[onto]][["N"]], mat.gsa.pair[[onto]][["index"]][,tag], xlab = "block size", ylab = "index", main = paste (onto, tag, sep = " : "))
    }
    dev.off ()
}


## GO INHIBITION CORRELATION
por <- 2
for (onto in ontologias) {
    print (onto)
    ##
    fichero.png <- paste0 ("paired_cor_", onto, ".png")
    png (fichero.png, width = por * 480, height = por * 480, pointsize = 12, bg = "white")
    plotcorr (cor (mat.gsa.pair[[onto]][["index"]]), main = paste (toupper (onto), "inhibition effect correlation"))
    dev.off ()
    ##
    fichero.png <- paste0 ("paired_dist_of_cor", onto, ".png")
    png (fichero.png, width = por * 480, height = por * 480, pointsize = 12 * 2, bg = "white")
    boxplot (cor (mat.gsa.pair[[onto]][["index"]]), las = 3, main = paste (toupper (onto), "inhibition effect correlation distribution"))
    abline (h = 0, col = "blue")
    dev.off ()
}


png ("paired_inhibition_effect_correlation_across_ontologies", width = por * 480 * 3, height = por * 480, pointsize = 12, bg = "white")
par (mfrow = c(1, 3))
plot (cor (mat.gsa.pair[["bp"]][["index"]]), cor (mat.gsa.pair[["cc"]][["index"]]))
abline (0, 1, col = "blue")
plot (cor (mat.gsa.pair[["bp"]][["index"]]), cor (mat.gsa.pair[["mf"]][["index"]]))
abline (0, 1, col = "blue")
plot (cor (mat.gsa.pair[["mf"]][["index"]]), cor (mat.gsa.pair[["cc"]][["index"]]))
abline (0, 1, col = "blue")
dev.off ()


###EXIT
warnings ()
sessionInfo ()
q ("no")
