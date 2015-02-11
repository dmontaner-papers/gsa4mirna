##d020_explore_sample_info.r
##2014-05-23 dmontaner@cipf.es
##Collecting miRNA data from The Cancer Genome Atlas
##Explore the sample info organization

## OBS: "Blood Derived Normal" samples do not have miRNA expression data

## We use just those patients with case and control measured by IlluminaHiSeq
## See: http://www.biostars.org/p/66612/

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
#library (); packageDescription ("", fields = "Version") #

try (source (".job.r")); try (.job)

options (width = 200)

###DATOS
load (file.path (.job$dir$proces, "sample_info_all.RData"))
ls ()

tags <- names (sinfos)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

tag <- "blca"

mat <- matrix (NA, nrow = length (tags), ncol = 3)
rownames (mat) <- tags
colnames (mat) <- c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    
    sinfo <- sinfos[[tag]]
    
    ## sinfo[,"path"] <- file.path (.job$dir$raw, TAGS[tag], "miRNASeq", paste0 ("BCGSC__", sinfo$Platform), "Level_3")
    ## sinfo[is.na (sinfo$Platform), "path"] <- NA
    
    sinfo[,"exists.mirna"]   <- file.exists (file.path (sinfo[,"path"], sinfo[,"File.Name.mirna"]))
    sinfo[,"exists.isoform"] <- file.exists (file.path (sinfo[,"path"], sinfo[,"File.Name.isoform"]))

    if (!identical (sinfo[,"exists.mirna"], sinfo[,"exists.isoform"])) stop ("mirna and isoform files do not match")
    
    table (sinfo[,"exists.mirna"], sinfo[,"exists.isoform"], sinfo$Platform, exclude = NULL)

    ## type
    ta <- table (sinfo$patient, sinfo$sample_type)
    suma <- colSums (ta)
    print (suma)
    
    mat[tag, c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")] <- suma[c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")]
}
mat[is.na (mat)] <- 0
mat

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    
    sinfo <- sinfos[[tag]]
    print (table (sinfo$Platform, exclude = NULL))
}


################################################################################


### SAMPLES WITH DATA
mat2 <- matrix (NA, nrow = length (tags), ncol = 3)
rownames (mat2) <- tags
colnames (mat2) <- c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")

for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    
    sinfo <- sinfos[[tag]]
    touse <- !is.na (sinfo$Platform.Type)
    sinfo <- sinfo[touse,]
    
    ## sinfo[,"path"] <- file.path (.job$dir$raw, TAGS[tag], "miRNASeq", paste0 ("BCGSC__", sinfo$Platform), "Level_3")
    ## sinfo[is.na (sinfo$Platform), "path"] <- NA
    
    sinfo[,"exists.mirna"]   <- file.exists (file.path (sinfo[,"path"], sinfo[,"File.Name.mirna"]))
    sinfo[,"exists.isoform"] <- file.exists (file.path (sinfo[,"path"], sinfo[,"File.Name.isoform"]))

    if (!identical (sinfo[,"exists.mirna"], sinfo[,"exists.isoform"])) stop ("mirna and isoform files do not match")
    
    table (sinfo[,"exists.mirna"], sinfo[,"exists.isoform"], sinfo$Platform, exclude = NULL)

    ## type
    ta <- table (sinfo$patient, sinfo$sample_type)
    suma <- colSums (ta)
    print (suma)
    
    mat2[tag, c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")] <- suma[c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")]
}
mat2[is.na (mat2)] <- 0
mat2

###EXIT
warnings ()
sessionInfo ()
q ("no")
