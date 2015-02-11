##s040_save_res_uvgsa_paired.r
##2014-12-07 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script collects saves all GSA results for the supplementary materials

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.7"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.2"
#library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mdgsa)
#help (package = mirbaseID)

try (source (".job.r")); try (.job)

corte <- 0.05
por <- 1

ontologias <- c ("bp", "cc", "mf")

source (file.path (.job$dir$scripts, "000_function_signifCount_1.r"))

###DATOS
load (file.path (.job$dir$proces, "res_gsa_paired.RData"))
ls ()
names (res.gsa.pair)

tags <- names (res.gsa.pair[["bp"]])
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS


################################################################################

### Count significant
n.bp.pair <- signifCount (res.gsa.pair[["bp"]], stat = "lor", pval = "padj", cutoff = corte, colNames = c("bp.Down.pair", "bp.noDif.pair", "bp.UP.pair"))
n.cc.pair <- signifCount (res.gsa.pair[["cc"]], stat = "lor", pval = "padj", cutoff = corte, colNames = c("cc.Down.pair", "cc.noDif.pair", "cc.UP.pair"))
n.mf.pair <- signifCount (res.gsa.pair[["mf"]], stat = "lor", pval = "padj", cutoff = corte, colNames = c("mf.Down.pair", "mf.noDif.pair", "mf.UP.pair"))
##
n.go.pair <- n.bp.pair + n.cc.pair + n.mf.pair
colnames (n.go.pair) <- c ("go.Down.pair", "go.noDif.pair", "go.UP.pair")

n.bp.pair
n.cc.pair
n.mf.pair
n.go.pair

## SAVE
save (list = c ("n.bp.pair", "n.cc.pair", "n.mf.pair", "n.go.pair"), file = file.path (.job$dir$proces, "report", "uvgsa_counts_pair.RData"))

################################################################################


###SAVE xls
setwd (.job$dir$res)

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    for (onto in ontologias) {
        cat ("\n===== ", onto, " =====\n")
        fichero <- paste0 ("res_gsa_paired_", tag, "_", onto, ".xlsx")
        ##
        datos <- res.gsa.pair[[onto]][[tag]]
        orden <- order (datos$index, decreasing = TRUE)
        datos <- datos[orden, c ("N", "lor", "pval", "padj", "Name")]
        ##
        write.xlsx2 (datos, file = fichero)
    }
}


###EXIT
warnings ()
sessionInfo ()
q ("no")
