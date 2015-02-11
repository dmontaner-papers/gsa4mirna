##d040_format_paired_samples.r
##2014-05-30 dmontaner@cipf.es
##Collecting miRNA data from The Cancer Genome Atlas
##Formatting paired samples

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
#library (); packageDescription ("", fields = "Version") #

try (source (".job.r")); try (.job)

options (width = 200)

###DATOS
load (file.path (.job$dir$proces, "sample_info_all.RData"))     #sinfos
load (file.path (.job$dir$proces, "sample_info_paired.RData"))  #pinfos
ls ()

#tags <- names (sinfos)
tags <- names (pinfos)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

tag <- "blca"

paired.sinfo.s <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    sinfo <- sinfos[[tag]]
    pinfo <- pinfos[[tag]]
    print (colSums (!is.na (pinfo[,-1])))

    ##reformat to 2 columns
    dat.t <- cbind (pinfo[,c("patient", "PTs")], tumor = 1)
    dat.n <- cbind (pinfo[,c("patient", "SNs")], tumor = 0)
    colnames (dat.t) <- colnames (dat.n) <- c ("patient", "sample", "tumor")
    dat <- rbind (dat.t, dat.n)

    ##cut samples
    touse <- sinfo$Sample %in% dat$sample
    misinfo <- sinfo[touse, c ("Sample", "File.Name.isoform", "path")]
    dups <- duplicated (misinfo[,"Sample"])
    misinfo <- misinfo[!dups,]
    rownames (misinfo) <- misinfo[,"Sample"]

    ##combine sample and file
    dat[,c("fichero", "path")] <- misinfo[dat$sample, c("File.Name.isoform", "path")]
    rownames (dat) <- dat$sample

    ##store
    paired.sinfo.s[[tag]] <- dat
}

sapply (paired.sinfo.s, dim)

dat[1:3,]

###SALVAMOS
save (list = "paired.sinfo.s", file = file.path (.job$dir$proces, "sample_info_paired_clean.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
