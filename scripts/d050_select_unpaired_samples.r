##d050_select_unpaired_samples_snt.r
##2014-05-23 dmontaner@cipf.es
##Collecting miRNA data from The Cancer Genome Atlas
##Select the samples to be analyzed in an unpaired analysis
## "Primary Tumor" samples are compared against "Solid Tissue Normal" as in the paired analysis
## Remember there is no miRNA expression data available for "Blood Derived Normal" samples

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

new.sinfos <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    sinfo <- sinfos[[tag]]
    ## Use just IlluminaHiSeq
    touse <- which (sinfo$Platform == "IlluminaHiSeq_miRNASeq")
    sinfo <- sinfo[touse,]
    ## Use just some types
    touse <- sinfo$sample_type %in% c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")
    sinfo <- sinfo[touse,]
    ## Keep just samples with data
    sinfo <- sinfo[!is.na (sinfo$File.Name.isoform),]
    ## Eliminate duplicated samples
    dup <- duplicated (sinfo$Sample)
    print (table (dup))
    sinfo <- sinfo[!dup,]
    ## Rownames
    rownames (sinfo) <- sinfo$Sample
    
    ## KEEP JUST DESIRED CONTROL TYPE
    touse <- sinfo$sample_type %in% c ("Primary Tumor", "Solid Tissue Normal")
    ##touse <- sinfo$sample_type %in% c ("Primary Tumor", "Blood Derived Normal")
    sinfo <- sinfo[touse,]
    ## Define tumor variale
    sinfo[,"tumor"] <- 1 * (sinfo$sample_type == "Primary Tumor")
    print (table (sinfo[,"tumor"], sinfo[,"sample_type"], exclude = NULL))

    ## STORE
    new.sinfos[[tag]] <- sinfo
}

sapply (new.sinfos, dim)

recuento <- sapply (new.sinfos, function (x) c (tumor = sum (x$tumor == 1), control = sum (x$tumor == 0)))
recuento
table (colnames (recuento) == names (new.sinfos))

touse <- (recuento[1,] != 0) &  (recuento[2,] != 0)
table (touse)

new.sinfos <- new.sinfos[touse]
names (new.sinfos)


##final counts
sapply (new.sinfos, function (x) c (tumor = sum (x$tumor == 1), control = sum (x$tumor == 0)))


###SALVAMOS
sinfos <- new.sinfos
save (list = "sinfos", file = file.path (.job$dir$proces, "sample_info_unpaired.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
