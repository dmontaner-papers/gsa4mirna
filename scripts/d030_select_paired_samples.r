##d030_select_paired_samples.r
##2014-05-23 dmontaner@cipf.es
##Collecting miRNA data from The Cancer Genome Atlas
##Select the samples to be analyzed

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

pinfos <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    sinfo <- sinfos[[tag]]
    ## Use just IlluminaHiSeq
    touse <- which (sinfo$Platform == "IlluminaHiSeq_miRNASeq")
    sinfo <- sinfo[touse,]
    ## Use just some types
    touse <- sinfo$sample_type %in% c ("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal")
    sinfo <- sinfo[touse,]
    ## Patients with more than one sample
    ta <- table (sinfo$patient)
    pacientes <- unique (names (ta)[ta > 1])
    touse <- sinfo$patient %in% pacientes
    sinfo <- sinfo[touse,]
    ## sample pairs
    pinfo <- NULL
    for (pa in pacientes) {
        pinfo <- rbind (pinfo,
                        c (patient = pa, 
                           PTs = sinfo[sinfo$patient == pa & sinfo$sample_type == "Primary Tumor",        "Sample"][1],
                           BDs = sinfo[sinfo$patient == pa & sinfo$sample_type == "Blood Derived Normal", "Sample"][1],
                           SNs = sinfo[sinfo$patient == pa & sinfo$sample_type == "Solid Tissue Normal",  "Sample"][1]))
    }
    pinfo <- pinfo[!is.na (pinfo[,"BDs"]) | !is.na (pinfo[,"SNs"]),]
    pinfo <- pinfo[!is.na (pinfo[,"PTs"]),]
    pinfo <- as.data.frame (pinfo, stringsAsFactors = FALSE)
    pinfos[[tag]] <- pinfo
}
pinfos

dims <- sapply (pinfos, dim)
dims

touse <- dims[1,] > 0
pinfos <- pinfos[touse]
names (pinfos)

###SALVAMOS
save (list = "pinfos", file = file.path (.job$dir$proces, "sample_info_paired.RData"))


###EXIT
warnings ()
sessionInfo ()
q ("no")
