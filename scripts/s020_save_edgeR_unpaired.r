##s020_save_edgeR_unpaired.r
##2014-07-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script exports the differential gene expression results at gene level

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"
library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mirbaseID)
#help (package = mdgsa)

try (source (".job.r")); try (.job)

source (file.path (.job$dir$scripts, "000_function_signifCount_1.r"))

corte <- 0.05

## DATA
load (file.path (.job$dir$proces, "res_dif_exp_unpaired.RData"))
ls ()

tags <- names (res.edger)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

## ANNOTATION
load (file.path (.job$dir$proces, "annot_for_unpaired_data.RData"))
names (annot)

anotacion <- c(annot[["bp"]], annot[["cc"]], annot[["mf"]])
length (anotacion)
sum (sapply (annot, length))
table (duplicated (names (anotacion)))

anotacion <- annotFilter (anotacion, minBlockSize = 10, maxBlockSize = 300)
gen2go <- revList (anotacion)
length (gen2go)
gen2go[1:2]

################################################################################

tag <- "blca"

res.edger.formatted <- list ()
n.target.genes <- matrix (NA, nrow = length (tags), ncol = 3)
colnames (n.target.genes) <- c ("Down", "Intersect", "Up")
rownames (n.target.genes) <- tags
n.target.genes

n.target.gos <- n.target.genes  ## directly targeted gos

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    
    mat <- res.edger[[tag]]$table

    ## gene level p-value and statistic
    pvalue    <- mat[,"PValue"]
    statistic <- mat[,"logFC"]
    names (pvalue) <- names (statistic) <- rownames (mat)
        
    ## p-value adjustment
    mat[,"PAdjusted"] <- p.adjust (pvalue, method = "BH")
    
    ## ranking index at miRNA level
    mat[,"RankingIndex"] <- pval2index (pval = pvalue, sign = statistic)    
    
    ## differential expression pattern
    mat[,"DiffExp"] <- uvPat (mat, cutoff = corte, pvalue = "PAdjusted", statistic = "logFC")

    ## mirbase ids
    mat[,"miRBase20"] <- mirIDmat[rownames (mat), "mirBase20"]
    
    ## order by ranking index
    orden <- order (mat[,"RankingIndex"], decreasing = TRUE)
    mat <- mat[orden,]
    
    ## up and down regulated genes
    up <- unique (unlist (acc2geneTSp[rownames (mat)[mat$DiffExp ==  1]]))
    do <- unique (unlist (acc2geneTSp[rownames (mat)[mat$DiffExp == -1]]))

    ## up and down regulated GOs
    upGOs <- unique (unlist (gen2go[up]))
    doGOs <- unique (unlist (gen2go[do]))
    
    ## SAVE
    res.edger.formatted[[tag]] <- mat
    n.target.genes[tag,] <- c (length (do), length (intersect (up, do)), length (up))
    ##
    n.target.gos[tag,] <- c (length (doGOs), length (intersect (upGOs, doGOs)), length (upGOs))
}

sapply (res.edger.formatted, dim)
t (sapply (res.edger.formatted, function (x) table (x[,"DiffExp"])))
signif.counts <- signifCount (res.edger.formatted, stat = "logFC", pval = "PAdjusted", cutoff = corte)
signif.counts
n.target.genes
n.target.gos

## SAVE
save (list = c("res.edger.formatted", "signif.counts", "n.target.genes", "n.target.gos"),
      file = file.path (.job$dir$proces, "report", "res_dif_exp_unpaired_formatted.RData"))

## SAVE XLS
setwd (.job$dir$res)
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    fichero <- paste0 ("res_edger_unpaired_", tag, ".xlsx")
    write.xlsx2 (res.edger.formatted[[tag]], file = fichero)
}


###EXIT
warnings ()
sessionInfo ()
q ("no")
