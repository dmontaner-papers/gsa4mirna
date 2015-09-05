##c010_compare_minra_level_vs_gene_level.r
##2015-09-01 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script compares the GSA analysis done at miRNA level (Godard's method) and the GSA analysis done at gene level using our transference index.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.2 (2015-08-14)"
#library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.0.0"
#library (openxlsx); packageDescription ("openxlsx", fields = "Version") #"3.0.0"
#help (package = mdgsa)

try (source (".job.r")); try (.job)

options (width = 170)

por <- 2
#corte <- 0.05

### DATA
setwd (file.path (.job$dir$proces))
load ("gsa_index_godard_unpaired_nofilter.RData")
load ("gsa_index_godard_unpaired.RData")
load ("gsa_index_godard_paired_nofilter.RData")
load ("gsa_index_godard_paired.RData")

load (file.path (.job$dir$data, "data_processed", "res_gsa_unpaired.RData"))
load (file.path (.job$dir$data, "data_processed", "res_gsa_paired.RData"))

table (colnames (mat.gsa.unpa[["bp"]][["index"]]) == colnames (mat.gsa.unpa[["cc"]][["index"]]))
table (colnames (mat.gsa.unpa[["bp"]][["index"]]) == colnames (mat.gsa.unpa[["mf"]][["index"]]))
table (colnames (mat.gsa.unpa[["bp"]][["index"]]) == colnames (imat.godard.unpa))
table (colnames (mat.gsa.unpa[["bp"]][["index"]]) == colnames (imat.godard.unpa.nf))

table (colnames (mat.gsa.pair[["bp"]][["index"]]) == colnames (mat.gsa.pair[["cc"]][["index"]]))
table (colnames (mat.gsa.pair[["bp"]][["index"]]) == colnames (mat.gsa.pair[["mf"]][["index"]]))
table (colnames (mat.gsa.pair[["bp"]][["index"]]) == colnames (imat.godard.pair))
table (colnames (mat.gsa.pair[["bp"]][["index"]]) == colnames (imat.godard.pair.nf))

################################################################################

## format unpa
gen.unpa    <- rbind (mat.gsa.unpa[["bp"]][["index"]], mat.gsa.unpa[["cc"]][["index"]], mat.gsa.unpa[["mf"]][["index"]])
mir.unpa    <- imat.godard.unpa
mir.unpa.nf <- imat.godard.unpa.nf

table (colnames (gen.unpa) == colnames (mir.unpa))
table (colnames (gen.unpa) == colnames (mir.unpa.nf))
table (rownames (mir.unpa) == rownames (mir.unpa.nf))

dim (gen.unpa)
dim (mir.unpa)

comunes.unpa <- intersect (rownames (gen.unpa), rownames (mir.unpa))
length (comunes.unpa)


## format pair
gen.pair    <- rbind (mat.gsa.pair[["bp"]][["index"]], mat.gsa.pair[["cc"]][["index"]], mat.gsa.pair[["mf"]][["index"]])
mir.pair    <- imat.godard.pair
mir.pair.nf <- imat.godard.pair.nf

table (colnames (gen.pair) == colnames (mir.pair))
table (colnames (gen.pair) == colnames (mir.pair.nf))
table (rownames (mir.pair) == rownames (mir.pair.nf))

dim (gen.pair)
dim (mir.pair)

comunes.pair <- intersect (rownames (gen.pair), rownames (mir.pair))
length (comunes.pair)

################################################################################


## EXPLORE CORRELATION
setwd (.job$dir$plots)
graphics.off ()
tags <- colnames (gen.unpa)
tags

tag <- "blca"
for (tag in tags) {
    print (tag)
    png (filename = paste0 ("gene_vs_mirna_level_gsa", tag, ".png"), width = 480 * por, height = 480 * por)
    par (mfrow = c (2, 2))
    ##
    try ({plot (gen.unpa[comunes.unpa, tag], mir.unpa   [comunes.unpa, tag], xlab = "gene level", ylab = "miRNA level",             main = paste (tag, "(unpa)")); title (sub = paste ("cor:", round (cor (gen.unpa[comunes.unpa, tag], mir.unpa   [comunes.unpa, tag], use = "pairwise.complete.obs"), 2))); abline (h = 0, v = 0, col = "blue")})
    try ({plot (gen.unpa[comunes.unpa, tag], mir.unpa.nf[comunes.unpa, tag], xlab = "gene level", ylab = "miRNA level (no filter)", main = paste (tag, "(unpa)")); title (sub = paste ("cor:", round (cor (gen.unpa[comunes.unpa, tag], mir.unpa.nf[comunes.unpa, tag], use = "pairwise.complete.obs"), 2))); abline (h = 0, v = 0, col = "blue")})
    ##
    try ({plot (gen.pair[comunes.pair, tag], mir.pair   [comunes.pair, tag], xlab = "gene level", ylab = "miRNA level",             main = paste (tag, "(pair)")); title (sub = paste ("cor:", round (cor (gen.pair[comunes.pair, tag], mir.pair   [comunes.pair, tag], use = "pairwise.complete.obs"), 2))); abline (h = 0, v = 0, col = "blue")})
    try ({plot (gen.pair[comunes.pair, tag], mir.pair.nf[comunes.pair, tag], xlab = "gene level", ylab = "miRNA level (no filter)", main = paste (tag, "(pair)")); title (sub = paste ("cor:", round (cor (gen.pair[comunes.pair, tag], mir.pair.nf[comunes.pair, tag], use = "pairwise.complete.obs"), 2))); abline (h = 0, v = 0, col = "blue")})
    ##
    dev.off ()
}

###EXIT
warnings ()
sessionInfo ()
q ("no")
