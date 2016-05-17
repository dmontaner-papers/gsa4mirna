##d010_read_expression.r
##2016-05-04 david.montaner@gmail.com
##Read Expression Data

### NOTES about the normalized data from; <http://seqanswers.com/forums/showthread.php?t=42911>
## It took me a while to get my head around this, since the column names in the rsem.genes/isoforms.results files don't match the default output of RSEM, neither the version they claim to have used nor the most current version.
## The (first) RSEM paper explains that the program calculates two values. One represent the (estimated) number of reads that aligned to a transcript. This value is not an integer because RSEM only reports a guess of how many ambiguously mapping reads belong to a transcript/gene. This number is what the TCGA slightly misleadingly calls raw counts.
## The scaled estimate value on the other hand is the estimated frequency of the gene/transcript amongst the total number of transcripts that were sequenced. Newer versions of RSEM call this value (multiplied by 1e6) TPM - Transcripts Per Million. It's closely related to FPKM, as explained on the RSEM website. The important point is that TPM, like FPKM, is independent of transcript length, whereas "raw" counts are not!
## The *.normalized_results files on the other hand just contain a scaled version of the raw_counts column. The values are divided by the 75-percentile and multiplied by 1000. This should make the values a bit more comparable between experiments. The Perl code for this quantile normalisation can be found here.
## In conclusion, I would strongly recommend using the TPM/scaled_estimate values for all intents and purposes. It seems to me to be the more robust and mathematically sound value.
## Hope that helps, best wishes,
## Benjamin


date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.5 (2016-04-14)"

try (source (".job.r")); try (.job)

options (width = 240)

##ftag <- ".rsem.genes.results"
ftag <- ".rsem.genes.normalized_results"

##ecol <- "raw_count"
ecol <- "normalized_count"


### SAMPLE INFO
sinfo <- read.table (file = file.path (.job$dir$rawdat, "KICH", "METADATA", "UNC__IlluminaHiSeq_RNASeqV2", "unc.edu_KICH.IlluminaHiSeq_RNASeqV2.1.2.0.sdrf.txt"),
                     header = TRUE, sep = "\t", quote = "", as.is = TRUE)
table (sapply (sinfo, class))
dim (sinfo)
colnames (sinfo)

head (sinfo[,c ("Extract.Name", "Comment..TCGA.Barcode.", "Derived.Data.File")])
#sinfo[,"Extract.Name"]
#sinfo[,"Comment..TCGA.Barcode."]

touse <- grep (ftag, sinfo[,"Derived.Data.File"])
sinfo <- sinfo[touse,]
dim (sinfo)

table (duplicated (sinfo[,"Comment..TCGA.Barcode."])) ## OK no dups

sinfo[,"sample"] <- substring (sinfo[,"Comment..TCGA.Barcode."], 1, 15)



### COUNT DATA
setwd (file.path (.job$dir$rawdat, "KICH", "RNASeqV2", "UNC__IlluminaHiSeq_RNASeqV2", "Level_3"))
ficheros <- dir (pattern = ftag)
head (ficheros)
length (ficheros)

table (ficheros %in% sinfo[,"Derived.Data.File"])
table (ficheros == sinfo[,"Derived.Data.File"])   ### NOT IN THE SAME ORDER


## read first file
i <- 1

tabla <- read.table (file = sinfo[i,"Derived.Data.File"], header = TRUE, sep = "\t", quote = "", as.is = TRUE)
print (dim (tabla))
head (tabla)

gexp <- tabla[ecol]
rownames (gexp) <- tabla[,'gene_id']
colnames (gexp) <- sinfo[i,"sample"]

for (i in 2:nrow (sinfo)) {
    print (i)
    tabla <- read.table (file = sinfo[i,"Derived.Data.File"], header = TRUE, sep = "\t", quote = "", as.is = TRUE)
    if (any (rownames (gexp) != tabla[,'gene_id'])) {
        stop ('GENES DO NOT MATCH')
    }
    gexp[,sinfo[i,"sample"]] <- tabla[,ecol]
}

gexp[1:5, 1:3]
summary (unlist (gexp))
##boxplot (log (1 + gexp))


### REORDER AS MIRNA DATA
load (file.path (.job$dir$data, "data_processed", "datos_unpaired.RData")) ## datos
ls ()

dim (datos[["kich"]]$contmat)
table (colnames (datos[["kich"]]$contmat) == rownames (datos[["kich"]]$sinfo))
colnames (datos[["kich"]]$contmat)

table (colnames (datos[["kich"]]$contmat) %in% colnames (gexp))
table (colnames (datos[["kich"]]$contmat) == colnames (gexp))  ## not in the same order

## reorder expression data
gexp <- gexp[,colnames (datos[["kich"]]$contmat)]
table (colnames (datos[["kich"]]$contmat) == colnames (gexp))  ## OK


###SAVE
setwd (.job$dir$proces)
save (gexp, file = "kich_expression.RData")


###EXIT
warnings ()
sessionInfo ()
q ("no")
