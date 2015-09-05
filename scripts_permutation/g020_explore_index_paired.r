##g020_explore_index_paired.r
##2014-06-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script explores the computed gene ranking index from the results of the differential gene expression

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.1.0 (2014-04-10)"
library (ellipse); packageDescription ("ellipse", fields = "Version") #"0.3-8"
##library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"

try (source (".job.r")); try (.job)

options (width = 200)

################################################################################

## DATOS
setwd (file.path (.job$dir$proces))

load ("rindex0_paired.RData")
load ("rindexT_paired.RData")
load ("rindex_paired.RData")

ls ()

tags <- names (rindex)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

## MATRIX FORMAT

ids0 <- idsT <- idsN <- NULL
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")

    ids0 <- c (ids0, names (rindex0[[tag]]))
    idsT <- c (idsT, names (rindexT[[tag]]))
    idsN <- c (idsN, names (rindex [[tag]]))
}
ids0 <- sort (unique (ids0))
idsT <- sort (unique (idsT))
idsN <- sort (unique (idsN))

length (ids0)
length (idsT)
length (idsN)

table (idsT == idsN) ##ok


mat0         <- matrix (NA, nrow = length (ids0), ncol = length (rindex), dimnames = list (ids0, names (rindex)))
matT <- matN <- matrix (NA, nrow = length (idsT), ncol = length (rindex), dimnames = list (idsT, names (rindex)))
##
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    ##
    mat0[,tag] <- rindex0[[tag]][rownames (mat0)]
    matT[,tag] <- rindexT[[tag]][rownames (matT)]
    matN[,tag] <- rindex [[tag]][rownames (matN)]
}

summary (mat0)
summary (matT)
summary (matN)

################################################################################

## PLOTS
graphics.off ()
setwd (.job$dir$plots)
dir ()

## boxplot
por <- 2
cex <- 2
png ("paired_rindex_boxplot.png", width = por * 480, height = por * 480 * 2, pointsize = 12, bg = "white")
par (mfrow = c (3,1), cex.axis = cex, cex.lab = cex, cex.main = cex)
boxplot (mat0, las = 3, main = "rindex0"); abline (h = 0, col = "blue")
boxplot (matT, las = 3, main = "rindexT"); abline (h = 0, col = "blue")
boxplot (matN, las = 3, main = "rindex"); abline (h = 0, col = "blue")
dev.off ()


## correlation
por <- 2

corre <- cor (matT, matN)
diag (corre)  ##OK

corre0 <- cor (mat0, use = "pairwise.complete.obs")
round (corre0, 2)
png ("paired_cor_rindex0.png", width = por * 480, height = por * 480, pointsize = 12, bg = "white")
plotcorr (corre0, main = "Ranking index correlation (miRNAs)")
dev.off ()

correT <- cor (matT, use = "pairwise.complete.obs")
round (correT, 2)
png ("paired_cor_rindexT.png", width = por * 480, height = por * 480, pointsize = 12, bg = "white")
plotcorr (correT, main = "Transferred ranking index correlation")
dev.off ()

correN <- cor (matN, use = "pairwise.complete.obs")
round (correN, 2)
png ("paired_cor_rindexN.png", width = por * 480, height = por * 480, pointsize = 12, bg = "white")
plotcorr (correN, main = "Ranking index correlation (genes)")
dev.off ()


## correlation vs correlationx
tri0 <- corre0[upper.tri (corre0)]
triT <- correT[upper.tri (correT)]
triN <- correN[upper.tri (correN)]

por <- 1
png ("paired_rindex_cor_vs_cor.png", width = por * 480 * 3, height = por * 480, pointsize = 12, bg = "white")
par (mfrow = c (1, 3))
plot (tri0, triT); abline (0, 1, col = "red")
plot (tri0, triN); abline (0, 1, col = "red")
plot (triT, triT); abline (0, 1, col = "red")
dev.off ()


###EXIT
warnings ()
sessionInfo ()
q ("no")
