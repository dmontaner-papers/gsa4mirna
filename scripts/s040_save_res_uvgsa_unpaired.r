##s040_save_res_uvgsa_unpaired.r
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
load (file.path (.job$dir$proces, "res_gsa_unpaired.RData"))
ls ()
names (res.gsa.unpa)

tags <- names (res.gsa.unpa[["bp"]])
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS


################################################################################

### Count significant
n.bp.unpa <- signifCount (res.gsa.unpa[["bp"]], stat = "lor", pval = "padj", cutoff = corte, colNames = c("bp.Down.unpa", "bp.noDif.unpa", "bp.UP.unpa"))
n.cc.unpa <- signifCount (res.gsa.unpa[["cc"]], stat = "lor", pval = "padj", cutoff = corte, colNames = c("cc.Down.unpa", "cc.noDif.unpa", "cc.UP.unpa"))
n.mf.unpa <- signifCount (res.gsa.unpa[["mf"]], stat = "lor", pval = "padj", cutoff = corte, colNames = c("mf.Down.unpa", "mf.noDif.unpa", "mf.UP.unpa"))
##
n.go.unpa <- n.bp.unpa + n.cc.unpa + n.mf.unpa
colnames (n.go.unpa) <- c ("go.Down.unpa", "go.noDif.unpa", "go.UP.unpa")

n.bp.unpa
n.cc.unpa
n.mf.unpa
n.go.unpa

## SAVE
save (list = c ("n.bp.unpa", "n.cc.unpa", "n.mf.unpa", "n.go.unpa"), file = file.path (.job$dir$proces, "report", "uvgsa_counts_unpa.RData"))

################################################################################


###SAVE xls
setwd (.job$dir$res)

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    for (onto in ontologias) {
        cat ("\n===== ", onto, " =====\n")
        fichero <- paste0 ("res_gsa_unpaired_", tag, "_", onto, ".xlsx")
        ##
        datos <- res.gsa.unpa[[onto]][[tag]]
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
