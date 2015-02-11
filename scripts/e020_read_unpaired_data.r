##e020_read_unpaired_data.r
##2014-06-04 dmontaner@cipf.es
##Collecting miRNA data from The Cancer Genome Atlas
##Reading data UNpaired data

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
#library (); packageDescription ("", fields = "Version") #

try (source (".job.r")); try (.job)

options (width = 200)

###DATOS
load (file.path (.job$dir$proces, "sample_info_unpaired.RData"))
ls ()

tags <- names (sinfos)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################

tag <- "blca"

datos <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")
    sinfo <- sinfos[[tag]]
    conteos <- list ()

    sinfo[,"fichero"] <- sinfo[,"File.Name.isoform"]  ## to reuse the paired samples code
    sinfo[,"sample"]  <- sinfo[,"Sample"]
    
    ## Reading files
    for (i in 1:nrow (sinfo)) {
        fichero <- file.path (sinfo[i,"path"], sinfo[i,"fichero"])
        dat <- read.table (fichero, header = TRUE, sep = "\t", quote = "", as.is = TRUE)
        dat[1:3,]
        
        touse <- grep (",", dat$miRNA_region)
        malo <- setdiff (1:nrow (dat), touse)
        ##print (table (dat[malo, "miRNA_region"]))
        
        dat <- dat[touse,]
        dat[,"ID"] <- sub ("mature,", "", dat$miRNA_region)
        dat[,"ID"] <- sub ("star,",   "", dat[,"ID"])
        res.by <- by (dat$read_count, dat$ID, sum)
        
        conteos[[sinfo[i, "sample"]]] <- unlist (as.list (res.by))
    }
    ## t (sapply (conteos, summary)) ## minimum number of counts is 1
    
    ## MATRIX format
    misIDs <- unique (unlist (lapply (conteos, names)))
    length (misIDs)
    misIDs[1:3]
    
    contmat <- matrix (NA, nrow = length (misIDs), ncol = length (conteos))
    rownames (contmat) <- misIDs
    colnames (contmat) <- names (conteos)
    ##
    for (id in names (conteos)) {
        ##print (id)
        contmat[,id] <- conteos[[id]][misIDs]
    }
    summary (contmat)
    contmat[is.na (contmat)] <- 0
    summary (contmat)

    if (any (colnames (contmat) != rownames (sinfo))) stop ("PROBLEMA 1")

    ## FILTERING:
    ## Use those miRNAs having more than 10 counts in any sample
    ## or a sum big enough:
    ## higher than 50 or 
    ## higher than half the number of samples
    ## (ie. on average a count in cases or in controls)
    maxi <- apply (contmat, 1, max)
    touse.maxi <- maxi > 10
    ##
    sumas <- rowSums (contmat)
    touse.sumas <- sumas > min (50, ncol (contmat) / 2)
    ##
    print (table (touse.sumas, touse.maxi))
    ##
    touse <- touse.sumas | touse.maxi
    table (touse)
    ##
    contmat <- contmat[touse,]
    ##
    table (colnames (contmat) == rownames (sinfo))

    if (any (colnames (contmat) != rownames (sinfo))) stop ("PROBLEMA 1")
    
    ## store
    datos[[tag]] <- list (contmat = contmat, sinfo = sinfo)
}

names (datos)
sapply (datos, names)

t (sapply (datos, function (x) dim (x$contmat)))


###SALVAMOS
save (list = "datos", file = file.path (.job$dir$proces, "datos_unpaired.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
