##i010_format_res_uvgsa_unpaired.r
##2014-06-03 dmontaner@cipf.es
##2014-11-26 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script collects all GSA results

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.99.2"
#library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mdgsa)
#help (package = mirbaseID)

try (source (".job.r")); try (.job)

options (width = 170)

corte <- 0.05

ontologias <- c ("bp", "cc", "mf")


## INDEX
load (file.path (.job$dir$proces, "rindex_unpaired.RData"))
ls ()

tags <- names (rindex)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS


### DATA DESCRIPTION
estudios <- read.xlsx (file = file.path (.job$dir$docs, "lista_estudios.xls"), sheetIndex = 1, stringsAsFactors = FALSE)
rownames (estudios) <- tolower (estudios$ID)
estudios

tags.desc <- estudios[tags, "Available.Cancer.Types"]
tags.desc <- sapply (strsplit (tags.desc, split = " \\["), function (x) x[[1]])
names (tags.desc) <- tags
tags.desc

################################################################################

## GSA list format
setwd (file.path (.job$dir$proces, "res_uvgsa_unpaired"))
dir ()

res.gsa.unpa <- list ()
for (tag in tags) {
    cat ("\n", tag, fill = TRUE, sep = "")
    cat ("====================================================================================================", fill = TRUE)
    for (onto in ontologias) {
        fichero <- paste (tag, "_", onto, ".RData", sep = "")
        load (fichero)
        res[,"pat"] <- uvPat (res, cutoff = corte)
        res[,"index"] <- pval2index (pval = res[,"pval"], sign = res[,"lor"])
        res[,"Name"] <- getGOnames (res, verbose = FALSE)
        res.gsa.unpa[[onto]][[tag]] <- res
    }
}

res[1:3,]
res.gsa.unpa[[onto]][[tag]][1:3,]

## SOME CHECKINGS
sapply (res.gsa.unpa, length)

for (onto in ontologias) {  
    cat ("\n")
    print (onto)
    print (sapply (res.gsa.unpa[[onto]], dim))  ##OK all equal dimensions
    print (unique (sapply (res.gsa.unpa[[onto]], function (x) sum (rownames (x) == rownames (res.gsa.unpa[[onto]][[1]]))))) ##OK all equal rownames
    print (unique (sapply (res.gsa.unpa[[onto]], function (x) sum (x[,"N"] == res.gsa.unpa[[onto]][[1]][,"N"])))) ##OK all equal go sizes
    print (unique (sapply (res.gsa.unpa[[onto]], function (x) sum (x[,"conv"])))) ## OK all converged
}


################################################################################

### Matrix Format

mat.gsa.unpa <- list ()

for (onto in ontologias) {
    N <- res.gsa.unpa[[onto]][[1]][,"N"]
    names (N) <- rownames (res.gsa.unpa[[onto]][[1]])
    mat.gsa.unpa[[onto]][["N"]] <- N
    
    mat.gsa.unpa[[onto]][["lor"]]   <- sapply (res.gsa.unpa[[onto]], function (x) {s <- x[,"lor"];   names (s) <- rownames (x); s})
    mat.gsa.unpa[[onto]][["pval"]]  <- sapply (res.gsa.unpa[[onto]], function (x) {s <- x[,"pval"];  names (s) <- rownames (x); s})
    mat.gsa.unpa[[onto]][["padj"]]  <- sapply (res.gsa.unpa[[onto]], function (x) {s <- x[,"padj"];  names (s) <- rownames (x); s})
    mat.gsa.unpa[[onto]][["pat"]]   <- sapply (res.gsa.unpa[[onto]], function (x) {s <- x[,"pat"];   names (s) <- rownames (x); s})
    mat.gsa.unpa[[onto]][["index"]] <- sapply (res.gsa.unpa[[onto]], function (x) {s <- x[,"index"]; names (s) <- rownames (x); s})
}

mat.gsa.unpa[["tags"]] <- tolower (names (res.gsa.unpa[[onto]]))
mat.gsa.unpa[["TAGS"]] <- toupper (names (res.gsa.unpa[[onto]]))
mat.gsa.unpa[["desc"]] <- tags.desc[mat.gsa.unpa[["tags"]]]


###SALVAMOS
save (list = c ("res.gsa.unpa", "mat.gsa.unpa"), file = file.path (.job$dir$proces, "res_gsa_unpaired.RData"))



###EXIT
warnings ()
sessionInfo ()
q ("no")
