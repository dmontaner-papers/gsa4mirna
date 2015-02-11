##d010_read_sample_info.r
##2014-05-23 dmontaner@cipf.es
##Collecting miRNA data from The Cancer Genome Atlas
##This script organizes all sample information

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
#library (); packageDescription ("", fields = "Version") #

try (source (".job.r")); try (.job)

options (width = 200)


## FUNCTION
barcode2 <- function (bc, to = "sample") {
    y <- sapply (bc, strsplit, split = "-")
    lon <- sapply (y, length)
    L <- unique (lon)
    if (length (L) != 1) stop ("unequal pa terns")
    mat <- matrix (unlist (y), ncol = L, byrow = TRUE)
    sample.number <- substr (mat[,4], start = 1, stop = 2)
    if (to == "sample") {
        res  <- paste (mat[,1], mat[,2], mat[,3], sample.number, sep = "-")
    } else { ##patient
        res <- paste (mat[,1], mat[,2], mat[,3], sep = "-")
    }
    return (res)
}


###DATOS
setwd (.job$dir$raw)
TAGS <- dir ()
tags <- tolower (TAGS)
names (TAGS) <- tags
tags
TAGS

tag <- "blca"

################################################################################
## File Info: information about the files available
################################################################################
sinfos <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    
    ## READ FILE NAMES
    finfo <- read.table (file.path (.job$dir$rawdat, TAGS[tag], "file_manifest.txt"), header = TRUE, sep = "\t", quote = "", as.is = TRUE)
    ## eliminate duplicated lines (just in case)
    finfo <- unique (finfo)
    ## keep just miRNA files
    finfo <- finfo[grep ("miRNASeq", finfo$Platform.Type),]  
    ## separate isoform and mirna files
    finfo.i <- finfo[grep ("isoform", finfo$File.Name),]
    finfo.m <- finfo[grep ("mirna",   finfo$File.Name),]
    ## common columns
    comcol <- colnames (finfo) != "File.Name"
    if (any (finfo.i[,comcol] != finfo.m[,comcol])) stop ("Differences in isoform and mirna files")
    ## reshape data to single lines
    finfo <- cbind (finfo.i[,comcol], File.Name.isoform = finfo.i[,"File.Name"], File.Name.mirna = finfo.m[,"File.Name"], stringsAsFactors = FALSE)
    finfo[1:3,]
    
    ## READ SAMPLE DATA
    fichero <- file.path (.job$dir$rawdat, TAGS[tag], "Clinical", "Biotab", paste0 ("nationwidechildrens.org_biospecimen_sample_", tag, ".txt"))
    ## format col names
    columnas <- read.table (fichero, header = FALSE, sep = "\t", quote = "", colClasses = "character", nrow = 1)
    columnas <- unlist (columnas, use.names = FALSE)
    columnas <- make.names (columnas)
    ## read data
    sinfo <- read.table (fichero, header = TRUE, sep = "\t", quote = "", as.is = TRUE, na.strings = "[Not Available]", skip = 1)
    ## rename columns
    if (colnames (sinfo)[1] == "CDE_ID.2673864") {
        colnames (sinfo) <- columnas
    } else {
        stop ("The file does not have a second row of colnames")
    }
    ## eliminate duplicated lines (just in case)
    sinfo <- unique (sinfo)
    ## include sample and patient
    sinfo[,"Sample"]  <- barcode2 (sinfo$bcr_sample_barcode, to = "sample")
    sinfo[,"patient"] <- barcode2 (sinfo$bcr_sample_barcode, to = "patient")
    ## keep just sample and type information
    dim (sinfo)
    sinfo <- unique (sinfo[,c("Sample", "patient", "sample_type")])
    dim (sinfo)
    
    ## MERGE
    if (!all (table (finfo[,"Sample"] %in% sinfo[,"Sample"]))) stop ("Some files are not reported in the biospecimen sample file")
    sinfo <- merge (sinfo, finfo, all = TRUE, sort = TRUE, stringsAsFactors = FALSE) ##sort: logical.  Should the result be sorted on the ‘by’ columns?

    ## path
    sinfo[,"path"] <- file.path (.job$dir$raw, TAGS[tag], "miRNASeq", paste0 ("BCGSC__", sinfo$Platform), "Level_3")
    sinfo[is.na (sinfo$Platform), "path"] <- NA
    
    ## STORE DATA
    sinfos[[tag]] <- sinfo
}


###SALVAMOS
save (list = "sinfos", file = file.path (.job$dir$proces, "sample_info_all.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
