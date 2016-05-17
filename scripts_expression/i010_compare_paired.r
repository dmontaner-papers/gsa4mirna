##i010_compare_paired.r
##2016-05-17 dmontaner@cipf.es
##Compare results with and without filtering

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.0.2 (2013-09-25)"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.3"
#library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"
#help (package = mdgsa)
#help (package = mirbaseID)

try (source (".job.r")); try (.job)

options (widht = 240)

#dir (file.path (.job$dir$data, "data_processed", "res_uvgsa_paired"))
#dir (file.path (.job$dir$proces, "res_uvgsa_paired"))


tag <- "kich"

res.u <- NULL
res.f <- NULL
##
for (onto in c ("bp", "cc", "mf")) {
    load (file.path (.job$dir$data, "data_processed", "res_uvgsa_paired", paste0 (tag, "_", onto, ".RData"))) ## res
    res[,"pat"] <- uvPat (res)
    res.u <- rbind (res.u, res)
    ##
    load (file.path (.job$dir$proces, "res_uvgsa_paired", paste0 (tag, "_", onto, ".RData"))) ## res
    res[,"pat"] <- uvPat (res)
    res.f <- rbind (res.f, res)
}
rm (res)

table (rownames (res.u) == rownames (res.f))

res.u[,"name"] <- getGOnames (res.u)
res.f[,"name"] <- getGOnames (res.f)

head (res.u)
head (res.f)

plot (res.u$lor, res.f$lor)
cor (res.u$lor, res.f$lor)

cor.test (res.u$lor, res.f$lor)
cor.test (pval2index (res.u$pval, res.u$lor),
          pval2index (res.f$pval, res.f$lor))

table (u = res.u$pat, f = res.f$pat)

############################################################

touse <- res.u$pat != 0
table (touse)

res.u[touse,]
res.f[touse,]

############################################################

touse <- res.f$pat != 0
table (touse)

res.u[touse,]
res.f[touse,]

############################################################

goLeaves (res.u)
goLeaves (res.f)


###EXIT
warnings ()
sessionInfo ()
q ("no")
