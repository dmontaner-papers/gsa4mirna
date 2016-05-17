##d020_select_expressed_genes.r
##2016-05-04 david.montaner@gmail.com
##Select Expressed Genes

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.5 (2016-04-14)"

try (source (".job.r")); try (.job)

options (width = 130)
options (width = 240)


## Ranking INDEX
load (file.path (.job$dir$data, "data_processed", "rindex_paired.RData")) ## rindex
ls ()
names (rindex)
rindex <- rindex[["kich"]]
head (rindex)



### EXPRESSION DATA
setwd (.job$dir$proces)
load ("kich_expression.RData")
dim (gexp)
gexp[100:103,1:5]

sp <- strsplit (rownames (gexp), split = "\\|")
table (sapply (sp, length))
sp[[100]]
rownames (gexp)[100]

gname <- sapply (sp, "[", 1)
extra <- sapply (sp, "[", 2)

length (intersect (names (rindex), gname)) ## use this one
length (intersect (names (rindex), extra))

table (names (rindex) %in% gname)
table (gname %in% names (rindex))




### SELECT EXPRESSED GENES
table (apply (gexp >= 1, 1, any))
table (apply (gexp >= 1, 1, all))
table (apply (gexp >  0, 1, any))
table (apply (gexp >  0, 1, all))

touse <- apply (gexp > 1, 1, all)
table (touse)

exp.genes <- gname[touse]
length (exp.genes)

table (names (rindex) %in% exp.genes)

### SAVE
save (exp.genes, file = "kich_expressed_genes.RData")

###EXIT
warnings ()
sessionInfo ()
q ("no")
