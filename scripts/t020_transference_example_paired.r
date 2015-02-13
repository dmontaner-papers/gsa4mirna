##t020_transference_example_paired.r
##2014-12-18 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##We revise how the transference to gene has worked

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.1.0 (2014-04-10)"
library (ellipse); packageDescription ("ellipse", fields = "Version") #"0.3-8"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"
library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mirbaseID)

try (source (".job.r")); try (.job)

options (width = 200)

por <- 1

################################################################################

## DATOS
setwd (file.path (.job$dir$proces))

load ("annot_for_unpaired_data.RData")

load ("rindex0_paired.RData")
load ("rindexT_paired.RData")
load ("rindex_paired.RData")

ls ()

tags <- names (rindex)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

sapply (rindex0, length)
sapply (rindexT, length)

universo.mirnas <- unique (unlist (lapply (rindex0, names)))
length (universo.mirnas)


################################################################################
## Transfered Index & number of targeting miRNAS
################################################################################

graphics.off ()
setwd (.job$dir$plots)

for (tag in tags) {
    print (tag)

    ind0 <- rindex0[[tag]]
    indT <- rindexT[[tag]]
    mirna2gene <- acc2geneTSp[names (ind0)]
    gene2mirna <- revList (mirna2gene)
    nmirnas <- sapply (gene2mirna, length)
    nmirnas <- nmirnas[names (indT)]
    table (is.na (nmirnas))
    
    fichero.png <- paste0 ("paired_explore_transfer_", tag, ".png")
    png (filename = fichero.png, width = por * 480 * 3, height = por * 480 *2)
    par (mfrow = c(1,2))
    ##
    boxplot (indT, ylab = "indexT (gene level)", main = TAGS[tag])
    abline (h = 0, col = "blue")
    ##
    plot (nmirnas, indT, xlab = "N miRNAs (log scale)", ylab = "indexT (gene level)", log = "x" , main = TAGS[tag])
    ##abline (h = 0, v = c(1, 2, 10, 20, 25), col = "blue")
    abline (h = 0, col = "blue")
    abline (h = summary (indT)[c("1st Qu.", "Median", "3rd Qu.")], col = "grey")
    ##
    dev.off ()
}

################################################################################
### DIAGRAM
################################################################################

## selected GO
tag <- "brca"
go <- "GO:0060053"
getGOnames (go)

ind0 <- rindex0[[tag]]
indT <- rindexT[[tag]]
mirna2gene <- acc2geneTSp[names (ind0)]
gene2mirna <- revList (mirna2gene)
nmirnas <- sapply (gene2mirna, length)
nmirnas <- nmirnas[names (indT)]
table (is.na (nmirnas))

genes <- annot[["cc"]][[go]]
table (genes %in% names (indT))

ordenados <- sort (rindexT[[tag]][genes], decreasing = TRUE)
ordenados
genes <- names (ordenados)
genes

colores <- rainbow (length (genes))

## list for the boxplot
lista <- list ("All" = rindex0[[tag]])
for (gen in genes) {
    mirnas <- gene2mirna[[gen]]
    lista[[gen]] <- rindex0[[tag]][mirnas]
}


graphics.off ()
por <- 2
#png (filename = file.path (.job$dir$code, "paper", "images", "diagram.png"), width = por * 480 * 1.5, height = por * 480)
pdf (file = file.path (.job$dir$code, "paper", "images", "diagram.pdf"), width = por * 7 * 1.5, height = por * 7)
##par (mfcol = c(2, 1))
layout (mat = matrix (c (1,1,1,1, 2,2,2,2,2,2,  3,3,3,3,3,3,  4,4,4,4), ncol = 2), widths = c(2,1), heights = 1)
par (mar = c (5.1, 6, 4.1, 2.1)) #A vector of the form ‘c(bottom, left, top, right)
##
xlimite <- c (0.1, length (genes) + 1.5)
plot (1:(length (genes) + 1), c (0, indT[genes]), col = c ("black", colores), pch = 15, xlim = xlimite, ylim = c (-20, 20 + max (indT)), ylab = "Transferred   Index   (gene level)", xlab = "", xaxt = "n",
      main = paste0 (getGOnames (go), " (", go, ")"),
      cex = 2, cex.main = 3, cex.axis = 1.7, cex.lab = 2)
axis (side = 1, at = 1:(length (genes)+1), labels = c ("All", genes), tick = TRUE, cex.axis = 1.7)
abline (h = 0, col = "blue", lty = 2)
text (1.5 + length (genes), y = 1300, labels = "B", col = "brown", cex = 4)
#abline (h = 1000)
##
boxplot (lista, ylab = "Ranking Index (miRNA level) computed from p-values", col = c("white", colores), xlim = xlimite, xlab = "genes", cex.axis = 1.7, cex.lab = 2)
abline (h = 0, col = "blue", lty = 2)
text (1:(length (genes)+1), y = -190, labels = c("Number of miRNAs:", nmirnas[genes]), cex = 1.5)
text (0, y = -120, labels = "Under expressed miRNAs", offset = 0, srt = 90, col = "darkgrey", cex = 1.5)
text (0, y =  140, labels = "Over expressed miRNAs",  offset = 0, srt = 90, col = "darkgrey", cex = 1.5)
text (1.5 + length (genes), y = 190, labels = "A", col = "brown", cex = 4)
##
boxplot (indT, indT[genes], names = c("All genes", go), ylab = "Transferred   Index   (gene level)", cex.axis = 1.7, cex.lab = 2)
abline (h = 0, col = "blue", lty = 2)
points (seq (1.25, 1.75, length.out = length (genes)), indT[genes], col = colores, pch = 15, cex = 2)
text (2.4, y = 1300, labels = "C", col = "brown", cex = 4)
#abline (h = 1000)
dev.off ()


###EXIT
warnings ()
sessionInfo ()
q ("no")
