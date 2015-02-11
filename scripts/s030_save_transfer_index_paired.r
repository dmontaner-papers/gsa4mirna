##s030_save_transfer_index_paired.r
##2014-07-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script saves the transferred index from miRNAs to GENE

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"

try (source (".job.r")); try (.job)

## DATA
load (file.path (.job$dir$proces, "rindexT_paired.RData")) #Transferred
ls ()
load (file.path (.job$dir$proces, "rindex_paired.RData"))  #Transferred and normalized
ls ()

class (rindex)
class (rindexT)

sapply (rindex,  class)
sapply (rindexT, class)

#boxplot (rindexT, las = 3)
#boxplot (rindex)

table (names (rindex) == names (rindexT))

tags <- names (rindex)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

tag <- "blca"

setwd (.job$dir$res)

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")

    if (any (names (rindex[[tag]]) != names (rindexT[[tag]]))) stop ("problem")

    mat <- as.data.frame (cbind (transfer.index = rindexT[[tag]], normalized.index = rindex[[tag]]))

    ## order
    orden <- order (mat[,"transfer.index"], decreasing = TRUE)
    mat <- mat[orden,]
    mat[1:10,]

    ## GENE ids
    mat <- cbind (gene = rownames (mat), mat, stringsAsFactors = FALSE)
    
    ## SAVE xls
    fichero <- paste0 ("transfer_index_paired_", tag, ".xlsx")
    write.xlsx2 (mat, file = fichero, row.names = FALSE)
}

###EXIT
warnings ()
sessionInfo ()
q ("no")
