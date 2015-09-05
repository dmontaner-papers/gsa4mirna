##h030_prepare_annotation_unpaired.r
##2014-06-03 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##Format GO annotation from Ensembl-Biomart

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"

try (source (".job.r")); try (.job)

##DATOS
setwd (file.path (.job$dir$proces))
##load ("rindex_unpaired.RData")
load (file.path (.job$dir$data, "data_processed", "rindex_unpaired.RData"))
ls ()

sapply (rindex, class)
sapply (rindex, head)
sapply (rindex, length)
sapply (rindex, summary)

sapply (rindex, function (x) table (x == 0))

genes <- unique (unlist (lapply (rindex, names)))
length (genes)


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
ensgo[1:3,]
length (unique (ensgo[,"HGNC.symbol"]))
length (unique (ensgo[,"GO.Term.Accession"]))



## LIST FORMAT
ensgo <- ensgo[, c("HGNC.symbol", "GO.Term.Accession")]
ensgo[1:3,]

system.time (annot <- annotMat2list (ensgo))
length (annot)


## PROPAGATE ONTOLOGY
system.time (annot <- propagateGO (annot))


## FILTERING: better to be done for each dataset
# annot <- annotFilter (annot, minBlockSize = 10, maxBlockSize = 500)
# length (annot)


## Split Ontologies
annot <- splitOntologies (annot, na.rm = TRUE, verbose = TRUE)
sapply (annot, length)

## SAVING
save (list = "annot", file = file.path (.job$dir$proces, "annot_for_unpaired_data.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
