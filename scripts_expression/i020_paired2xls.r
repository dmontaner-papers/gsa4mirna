##i020_paired2xls.r
##2016-05-17 dmontaner@cipf.es
##Save to xls

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.3"
library (openxlsx); packageDescription ("openxlsx", fields = "Version") #"3.0.0"
#help (package = mdgsa)

try (source (".job.r")); try (.job)

options (widht = 240)

tag <- "kich"

lista <- list ()
for (onto in c ("bp", "cc", "mf")) {
    load (file.path (.job$dir$proces, "res_uvgsa_paired", paste0 (tag, "_", onto, ".RData"))) ## res
    res[,"pat"] <- uvPat (res)
    res[,"name"] <- getGOnames (res)    
    lista[[paste0 (onto, '_all')]] <- res
    lista[[paste0 (onto, '_non_redundatn')]] <- goLeaves (res)
}

lapply (lista, head)


### SAVE xlsx
setwd (.job$dir$res)

write.xlsx (lista, file = "filtered_analysis_kich_paired.xlsx")


###EXIT
warnings ()
sessionInfo ()
q ("no")
