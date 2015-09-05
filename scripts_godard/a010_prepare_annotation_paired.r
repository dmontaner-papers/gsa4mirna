##a010_prepare_annotation_paired.r
##2015-08-28 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script prepares GO annotation at miRNA level to carry out the logistic analysis in the miRNA dimension as in the paper of Godard.

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.1.0 (2014-04-10)"
library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"
library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mirbaseID)

try (source (".job.r")); try (.job)

################################################################################

## DATA edgeR results
setwd (file.path (.job$dir$proces))
##load ("res_dif_exp_paired.RData")
load (file.path (.job$dir$data, "data_processed", "res_dif_exp_paired.RData"))

ids <- sapply (res.edger, function (x) rownames (x$table))
class (ids)
sapply (ids, length)

unimirna <- unique(unlist (ids))
length (unimirna)

## the annotated ones
intermirna <- intersect (unimirna, names (acc2geneTSp))
length (intermirna)  #just 250

mirna2gene <- acc2geneTSp[intermirna]
genes <- unique (unlist (mirna2gene))
length (genes)

length (unique (unlist (acc2geneTSp)))


################################################################################

## ANNOTATION
ensgo <- read.table (file.path (.job$dir$annotation, "mart_export.txt"), 
                     header = TRUE, sep = "\t", quote = "", na.strings = "", colClasses = "character")
dim (ensgo)
ensgo[1:3,]

table (duplicated (ensgo))

table (ensgo$GO.Term.Accession == "")
table (ensgo$HGNC.symbol       == "")

na.go   <- is.na (ensgo$GO.Term.Accession)
na.gene <- is.na (ensgo$HGNC.symbol)

table (na.go, na.gene) ## some missing

touse <- !na.go & !na.gene
table (touse)

ensgo <- ensgo[touse,]
dim (ensgo)


## keep just genes in the dataset
touse <- ensgo[,"HGNC.symbol"] %in% genes  ##most of them are present
table (touse)

ensgo <- ensgo[touse,]
dim (ensgo)
ensgo[1:3,]
length (unique (ensgo[,"HGNC.symbol"]))
length (unique (ensgo[,"GO.Term.Accession"]))


## list format
gene2go <- annotMat2list(ensgo[, c("GO.Term.Accession", "HGNC.symbol")])
gene2go[1:3]

################################################################################

## miRNA to go annotation
mirna2go <- list ()
for (mi in names (mirna2gene)) {
    mirna2go[[mi]] <- unique (unlist (gene2go[mirna2gene[[mi]]]))
}

length (mirna2go)
names (mirna2go)
mirna2go[1:3]

lon <- sapply (mirna2go, length)
summary (lon)

################################################################################

## go annotation at MIRNA LEVEL
annot <- revList (mirna2go)
length (annot)
summary (sapply (annot, length))

## PROPAGATE ONTOLOGY
system.time (annot <- propagateGO (annot))
annot[1:3]
length (annot)
summary (sapply (annot, length))

## ## Split Ontologies
## annot <- splitOntologies (annot, na.rm = TRUE, verbose = TRUE)
## sapply (annot, length)
## sum (sapply (annot, length))

## SAVE
save (list = "annot",  file =  "annotation_paired.RData")


###EXIT
warnings ()
sessionInfo ()
q ("no")
