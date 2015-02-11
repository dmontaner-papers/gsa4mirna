##g010_prepare_index_paired.r
##2014-06-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script computes the gene ranking index from the results of the differential gene expression

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.1.0 (2014-04-10)"
library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"
library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
help (package = mirbaseID)

try (source (".job.r")); try (.job)

################################################################################

## miRNA to gene transfer information
class (acc2geneTSp)
length (acc2geneTSp)
lapply (acc2geneTSp[1:2], head)
dim (annotList2mat (acc2geneTSp))
##dim (annotList2mat (acc2geneTSc)) ##five times bigger... we do not use it. Introduces some noise...
length (unique (unlist (acc2geneTSp)))

################################################################################

## DATA edgeR results
setwd (file.path (.job$dir$proces))
load ("res_dif_exp_paired.RData")
ls ()

tags <- names (res.edger)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

tag <- "blca"

rindex0 <- rindexT <- rindex <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")

    mat <- res.edger[[tag]]$table
    
    ## gene level p-value and statistic
    pvalue    <- mat[,"PValue"]
    statistic <- mat[,"logFC"]
    names (pvalue) <- names (statistic) <- rownames (mat)
    
    ## Index
    rindex0[[tag]] <- pval2index (pval = pvalue, sign = statistic)
    rindex0[[tag]][1:3]

    ## Transform (to normal distribution) at miRNA level BETTER NOT USE
    ##rindex0[[tag]] <- indexTransform (index = rindex0[[tag]], method = "normalize")
    
    ## Transfer
    rindexT[[tag]] <- transferIndex (index = rindex0[[tag]], targets = acc2geneTSp, method = "sum") ## transfer before transformation (normalization)
    rindexT[[tag]][1:3]
    
    ## Transform (to normal distribution)
    rindex[[tag]] <- indexTransform (index = rindexT[[tag]], method = "normalize")
    rindex[[tag]][1:3]
}

t (sapply (rindex0, summary))
t (sapply (rindexT, summary))
t (sapply (rindex, summary))


## SAVE
save (list = "rindex0", file = "rindex0_paired.RData")
save (list = "rindexT", file = "rindexT_paired.RData")
save (list = "rindex",  file = "rindex_paired.RData")

###EXIT
warnings ()
sessionInfo ()
q ("no")
