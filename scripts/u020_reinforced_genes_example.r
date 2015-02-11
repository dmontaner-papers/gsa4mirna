##u020_reinforced_genes_example.r
##2015-01-23 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script searches for miRNAs with reinforced effect over a gene

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

tag <- "esca"
tags <- "esca"

graphics.off ()

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
    names (Nmirna) <- rownames (tmat)
    summary (Nmirna)
    Nmirna[1:3]
    misd <- function (x) sd (x[x!=0])
    sds <- apply (tmat, 1, misd)
    length (sds)
    summary (sds)
    sds[1:3]
    table (rownames (tmat) == names (sds))
    ##table (rownames (tmat) == names (Nmirna))
    mimax <- apply (tmat, 1, function (x) max (abs (x)))
    length (mimax)
    mimax[1:3]

    touse <- 9 < Nmirna & Nmirna < 20
    print (table (touse))
    
    ## x11 ()
    ## plot (mimax[touse], tind[touse], log = "x", main = tag)
    ## abline (h = c(0, -50, 50, 200, -200), col = "blue")
    ## abline (v = 0:10)
    ## #abline (h = -45, col = "blue")
    
    ## ## Transfer
    ## rindexT[[tag]] <- transferIndex (index = rindex0[[tag]], targets = acc2geneTSp, method = "sum") ## transfer before transformation (normalization)
    ## rindexT[[tag]][1:3]
    
    ## ## Transform (to normal distribution)
    ## rindex[[tag]] <- indexTransform (index = rindexT[[tag]], method = "normalize")
    ## rindex[[tag]][1:3]
}

touse <- 9 < Nmirna & Nmirna < 20 & 7 < mimax & mimax < 10 & tind < -45
table (touse)
rownames (tmat)[touse]
Nmirna["TRIB1"]
Nmirna["GREB1"]

u <- tmat["TRIB1",]
u <- u[u != 0]
sort (u)
sum (u)
tind["TRIB1"]

u <- tmat["GREB1",]
u <- u[u != 0]
sort (u)
sum (u)
tind["GREB1"]


###EXIT
warnings ()
sessionInfo ()
q ("no")
