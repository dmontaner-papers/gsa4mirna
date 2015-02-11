##u010_two_canceled_out_genes_example.r
##2015-01-23 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script searches for genes targeted for two miRNAs that cancel out

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.1.0 (2014-04-10)"
#library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"
library (Matrix)  ## needed for the rowSums
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

tags <- "kich"

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

    tind <- transferIndex (index = rindex0[[tag]], targets = acc2geneTSp, method = "sum") ## transfer before transformation (normalization)
    tmat <- transferIndex (index = rindex0[[tag]], targets = acc2geneTSp, method = "sum", transferMatrix = TRUE) ## transfer before transformation (normalization)
    dim (tmat)
    length (tind)

    table (rownames (tmat) == names (tind))
    tind[1:3]
    tmat[1:3, 1:5]
    rowSums (tmat[1:3,])
    rowSums (tmat[1:3, 1:5])
    rowSums (tmat[1:3, 1:5] == 0)
    
    Nmirna <- rowSums (tmat != 0)
    summary (Nmirna)
    Nmirna[1:3]
    misd <- function (x) sd (x[x!=0])
    sds <- apply (tmat, 1, misd)
    length (sds)
    summary (sds)
    sds[1:3]
    table (rownames (tmat) == names (sds))
    ##table (rownames (tmat) == names (Nmirna))

    ## par (mfrow = c(2,1))
    ## plot (Nmirna, tind)
    ## abline (h = 0, col = "blue")
    ## abline (v = 0:5)
    ## ##
    ## plot (Nmirna, sds)
    ## abline (h = 0, col = "blue")
    ## abline (v = 0:5)
   
    ## par (mfrow = c(2,1))
    ## plot (Nmirna, tind, log = "x")
    ## abline (h = 0, col = "blue")
    ## abline (v = 1:5)
    ## ##
    ## plot (Nmirna, sds, log = "x")
    ## abline (h = 0, col = "blue")
    ## abline (v = 1:5)

    touse <- Nmirna == 2
    table (touse)
    
    dosp <- t (apply (tmat[touse,], 1, function (x) x[x != 0]))
    ## x11()
    ## plot (dosp, main = tag)
    ## abline (0, -1, col = "blue")
    ## abline (0,  1, col = "blue")
    ## abline (h = 0, v = 0, col = "blue")

    touse <- dosp[,1] > 50 &dosp[,2] < -50
    table (touse)

    print (dosp[touse,])
    print (tind[rownames (dosp[touse,])])
    
    ## ## Transfer
    ## rindexT[[tag]] <- transferIndex (index = rindex0[[tag]], targets = acc2geneTSp, method = "sum") ## transfer before transformation (normalization)
    ## rindexT[[tag]][1:3]
    
    ## ## Transform (to normal distribution)
    ## rindex[[tag]] <- indexTransform (index = rindexT[[tag]], method = "normalize")
    ## rindex[[tag]][1:3]
}

touse <- tmat["GPR162",] != 0
table (touse)
tmat[c ("GPR162", "UBALD1"), touse]

intersect (acc2geneTSp[["MIMAT0000077"]], acc2geneTSp[["MIMAT0000271"]])

mirIDmat[c ("MIMAT0000077", "MIMAT0000271"), "mirBase20"]

rindex0[["kich"]][c ("MIMAT0000077", "MIMAT0000271")]

tind[c ("GPR162", "UBALD1")]

mat[c ("MIMAT0000077", "MIMAT0000271"),]


###EXIT
warnings ()
sessionInfo ()
q ("no")
