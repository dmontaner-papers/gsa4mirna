##f020_edgeR_unpaired.r
##2014-06-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script performs the differential gene expression analysis for the paired data

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"

try (source (".job.r")); try (.job)

## DATA
setwd (file.path (.job$dir$proces))
load ("datos_unpaired.RData")
ls ()

tags <- names (datos)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

sapply (datos, names)
sapply (datos, function (x) dim (x$contmat))

sapply (datos, function (x) table (colnames (x$contmat) == rownames (x$sinfo)))

################################################################################

tag <- "blca"
tag <- "skcm"

res.edger <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    sinfo   <- datos[[tag]]$sinfo
    contmat <- datos[[tag]]$contmat
    print (dim (contmat))
    
    ## compute FOLD CHANGES
    ## just for checking purposes
    ## no normalization has been done
    ## at least RPKM or similar will be needed before the fold change is computed
    media0 <- rowMeans (contmat[,sinfo$tumor == 0, drop = FALSE])
    media1 <- rowMeans (contmat[,sinfo$tumor == 1, drop = FALSE])
    mifold <- log2 ((media1 + 0.001) / (media0 + 0.001))
    
    ## edgeR classic analysis
    clase <- sinfo$tumor   ##does not need to be a factor
    edger.cla <- DGEList (counts = contmat, group = clase)
    edger.cla <- calcNormFactors     (edger.cla)       ## HAS TO BE COMPUTED IN THIS ORDER
    edger.cla <- estimateCommonDisp  (edger.cla)
    edger.cla <- estimateTagwiseDisp (edger.cla)
    edger.cla.res <- exactTest       (edger.cla)
    ##topTags (edger.cla.res)
    print (cor (mifold, edger.cla.res$table[,"logFC"], use = "pairwise.complete.obs"))
    
    ## ## PAIRED edgeR glm analysis
    ## indiv <- as.factor (sinfo$patient)
    ## clase <- sinfo$tumor
    ## design <- model.matrix (~ clase + indiv)
    ## rownames (design) <- colnames (contmat)
    ## edger.par <- DGEList (counts = contmat)
    ## edger.par <- calcNormFactors (edger.par)
    ## system.time (edger.par <- estimateGLMCommonDisp  (edger.par, design))  #takes some time
    ## system.time (edger.par <- estimateGLMTrendedDisp (edger.par, design))
    ## system.time (edger.par <- estimateGLMTagwiseDisp (edger.par, design))
    ## edger.par.fit <- glmFit (edger.par, design)
    ## edger.par.res <- glmLRT (edger.par.fit, coef = "clase")
    ## ##topTags (edger.par.res)
    ## print (cor (mifold, edger.par.res$table[,"logFC"], use = "pairwise.complete.obs"))
    
    ##STORE RESULTS
    ##res.edger[[tag]] <- edger.par.res
    res.edger[[tag]] <- edger.cla.res
}

names (res.edger)


## SAVING
save (list = "res.edger", file = file.path (.job$dir$proces, "res_dif_exp_unpaired.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
