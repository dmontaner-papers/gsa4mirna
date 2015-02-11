##s050_combine_study_info.r
##2014-07-07 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script combines information for all cancer studies

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
library (knitr); packageDescription ("knitr", fields = "Version") #"1.6"
#library (xtable); packageDescription ("xtable", fields = "Version") #"1.7-3"
#library (tables); packageDescription ("tables", fields = "Version") #"0.7.79"
#help (package = xtable)
#help (package = tables)

try (source (".job.r")); try (.job)

options (width = 170)


### DATA
setwd (file.path (.job$dir$proces, "report"))

load ("case_control_counts_paired.RData")
case.cont.paired <- case.cont

load ("case_control_counts_unpaired.RData")
datos <- case.cont

## combine
datos[,"N.paired.samples"] <- case.cont.paired[rownames (datos), "N.cases"]
datos[is.na (datos$N.paired.samples), "N.paired.samples"] <- 0L
dim (datos)
datos

################################################################################

tags <- rownames (datos)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS


################################################################################
## DATA DESCRIPTION: cancer type
################################################################################

estudios <- read.xlsx (file = file.path (.job$dir$docs, "lista_estudios.xls"), sheetIndex = 1, stringsAsFactors = FALSE)
rownames (estudios) <- tolower (estudios$ID)

tags.desc <- estudios[tags, "Available.Cancer.Types"]
names (tags.desc) <- tags
tags.desc.simple <- sub (" \\[.*","", tags.desc)
tags.desc.simple <- sub ("Cervical squamous cell carcinoma and endocervical adenocarcinoma", "Cervical squamous cell carcinoma", tags.desc.simple)

tags.desc
tags.desc.simple


## cancer type
datos[,"cancer.type"] <- tags.desc.simple[rownames (datos)]
datos

## id
table (rownames (datos) == tags)
#datos[,"ID"] <- TAGS[rownames (datos)]
datos <- cbind (ID = TAGS[rownames (datos)], datos, stringsAsFactors = FALSE)
sapply (datos, class)
datos

kable (datos)
#xls

################################################################################
## RESULTS AT miRNA LEVEL
################################################################################

## miRNA level results UNpaired
load (file.path (.job$dir$proces, "report", "res_dif_exp_unpaired_formatted.RData"))
signif.counts.unpa <- signif.counts
n.target.genes.unpa <- data.frame (n.target.genes)
n.target.gos.unpa   <- data.frame (n.target.gos)

table (rownames (datos) == rownames (signif.counts.unpa))
table (rownames (datos) == rownames (n.target.genes.unpa))
table (rownames (datos) == rownames (n.target.gos.unpa))


## miRNA level results Paired
load (file.path (.job$dir$proces, "report", "res_dif_exp_paired_formatted.RData"))
signif.counts.pair <- signif.counts[rownames (datos),]
n.target.genes.pair <- as.data.frame (n.target.genes)[rownames (datos),]
n.target.gos.pair   <- as.data.frame (n.target.gos)  [rownames (datos),]

rownames (signif.counts.pair) <- rownames (n.target.genes.pair) <- rownames (n.target.gos.pair) <- rownames (datos)

###

colnames (signif.counts.unpa) <- c ("miRNA.Down.unpa", "miRNA.noDif.unpa", "miRNA.Up.unpa")
colnames (signif.counts.pair) <- c ("miRNA.Down.pair", "miRNA.noDif.pair", "miRNA.Up.pair")

colnames (n.target.genes.unpa) <- paste0 ("targets.", colnames (n.target.genes.unpa), ".unpa")
colnames (n.target.genes.pair) <- paste0 ("targets.", colnames (n.target.genes.pair), ".pair")

colnames (n.target.gos.unpa) <- paste0 ("targetGOs.", colnames (n.target.gos.unpa), ".unpa")
colnames (n.target.gos.pair) <- paste0 ("targetGOs.", colnames (n.target.gos.pair), ".pair")

datos <- cbind (datos, signif.counts.unpa[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, signif.counts.pair[rownames (datos),], stringsAsFactors = FALSE)
##
datos <- cbind (datos, n.target.genes.unpa[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.target.genes.pair[rownames (datos),], stringsAsFactors = FALSE)
##
datos <- cbind (datos, n.target.gos.unpa[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.target.gos.pair[rownames (datos),], stringsAsFactors = FALSE)


################################################################################
## RESULTS AT GO LEVEL
################################################################################

load (file.path (.job$dir$proces, "report", "uvgsa_counts_unpa.RData"))
datos <- cbind (datos, n.bp.unpa[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.cc.unpa[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.mf.unpa[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.go.unpa[rownames (datos),], stringsAsFactors = FALSE)

load (file.path (.job$dir$proces, "report", "uvgsa_counts_pair.RData"))
datos <- cbind (datos, n.bp.pair[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.cc.pair[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.mf.pair[rownames (datos),], stringsAsFactors = FALSE)
datos <- cbind (datos, n.go.pair[rownames (datos),], stringsAsFactors = FALSE)

t (t (colnames (datos)))

datos


###SAVE
save (list = "datos", file = file.path (.job$dir$proces, "report", "allStudyInfo.RData"))

###SAVE xls
setwd (.job$dir$res)
write.xlsx2 (datos, file = "allStudyInfo.xlsx")


###EXIT
warnings ()
sessionInfo ()
q ("no")
