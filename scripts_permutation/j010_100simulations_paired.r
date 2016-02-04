##j010_100simulations_paired.r
##2016-01-29  fgarcia@cipf.es, dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script generates 100 simulations for the functional analysis in paired studies


date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.1.0 (2014-04-10)"
library (edgeR); packageDescription ("edgeR", fields = "Version") #"3.6.1"
library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.4"
library (mirbaseID); packageDescription ("mirbaseID", fields = "Version") #"0.0.2"


setwd("/home/paco/gsa4mirna/scripts_permutation")
try (source (".job.r")); try (.job)

## getting parameters 
args <- commandArgs(trailingOnly = TRUE)
nsim <- as.character(args[1])

set.seed (201508280 + as.numeric(nsim))


################################################################################

## miRNA to gene transfer information
class (acc2geneTSp)
length (acc2geneTSp)
lapply (acc2geneTSp[1:2], head)
dim (annotList2mat (acc2geneTSp))
##dim (annotList2mat (acc2geneTSc)) ##five times bigger... we do not use it. Introduces some noise...
length (unique (unlist (acc2geneTSp)))


### PERMUTE gene column
#head (acc2geneTSp)
mat <- annotList2mat (acc2geneTSp)
# head (mat)
dim(mat)
mat[,1] <- sample (mat[,1])
#head (mat)
acc2geneTSp <- annotMat2list (mat)
# head (acc2geneTSp)
# save (acc2geneTSp, file = file.path (.job$dir$proces, "permuted_targets_paired.Rdata"))

################################################################################

## DATA edgeR results
setwd (file.path (.job$dir$proces))
##load ("res_dif_exp_paired.RData")
load (file.path (.job$dir$data, "data_processed", "res_dif_exp_paired.RData"))
ls ()

tags <- names (res.edger)
TAGS <- toupper (tags)
names (TAGS) <- tags
tags
TAGS

################################################################################


rindex0 <- rindexT <- rindex <- list ()
for (tag in tags) {
    cat ("\n=============== ", TAGS[tag], " ===============\n")

    mat <- res.edger[[tag]]$table
    
    ## gene level p-value and statistic
    pvalue    <- mat[,"PValue"]
    statistic <- mat[,"logFC"]
    names (pvalue) <- names (statistic) <- rownames (mat)
    
    ## Index
    rindex0[[tag]] <- pval2index (pval = pvalue, sign = statistic)
    rindex0[[tag]][1:3]

    ## Transform (to normal distribution) at miRNA level BETTER NOT USE
    ##rindex0[[tag]] <- indexTransform (index = rindex0[[tag]], method = "normalize")
    
    ## Transfer
    rindexT[[tag]] <- transferIndex (index = rindex0[[tag]], targets = acc2geneTSp, method = "sum") ## transfer before transformation (normalization)
    rindexT[[tag]][1:3]
    
    ## Transform (to normal distribution)
    rindex[[tag]] <- indexTransform (index = rindexT[[tag]], method = "normalize")
    rindex[[tag]][1:3]
}

t (sapply (rindex0, summary))
t (sapply (rindexT, summary))
t (sapply (rindex, summary))


# ## SAVE
# save (list = "rindex0", file = "rindex0_paired.RData")
# save (list = "rindexT", file = "rindexT_paired.RData")
# save (list = "rindex",  file = "rindex_paired.RData")



################################################################################

##Format GO annotation from Ensembl-Biomart


##DATOS
# setwd (file.path (.job$dir$proces))
# ##load ("rindex_paired.RData")
# load (file.path (.job$dir$data, "data_processed", "rindex_paired.RData"))
# ls ()

sapply (rindex, class)
sapply (rindex, head)
sapply (rindex, length)
sapply (rindex, summary)

sapply (rindex, function (x) table (x == 0))

genes <- unique (unlist (lapply (rindex, names)))
length (genes)


## ANNOTATION
ensgo <- read.table (file.path (.job$dir$annotation, "mart_export.txt"), 
                     header = TRUE, sep = "\t", quote = "", na.strings = "", colClasses = "character")
dim (ensgo)
ensgo[1:3,]

table (duplicated (ensgo))

table (ensgo$GO.Term.Accession == "")
table (ensgo$HGNC.symbol       == "")

na.go   <- is.na (ensgo$GO.Term.Accession)
na.gene <- is.na (ensgo$HGNC.symbol)

table (na.go, na.gene) ## some missing

touse <- !na.go & !na.gene
table (touse)

ensgo <- ensgo[touse,]
dim (ensgo)


## keep just genes in the dataset
touse <- ensgo[,"HGNC.symbol"] %in% genes  ##most of them are present
table (touse)

ensgo <- ensgo[touse,]
ensgo[1:3,]
length (unique (ensgo[,"HGNC.symbol"]))
length (unique (ensgo[,"GO.Term.Accession"]))



## LIST FORMAT
ensgo <- ensgo[, c("HGNC.symbol", "GO.Term.Accession")]
ensgo[1:3,]

system.time (annot <- annotMat2list (ensgo))
length (annot)


## PROPAGATE ONTOLOGY
system.time (annot <- propagateGO (annot))


## FILTERING: better to be done for each dataset
# annot <- annotFilter (annot, minBlockSize = 10, maxBlockSize = 500)
# length (annot)


## Split Ontologies
annot <- splitOntologies (annot, na.rm = TRUE, verbose = TRUE)
sapply (annot, length)
sum (sapply (annot, length))

# ## SAVING
# save (list = "annot", file = file.path (.job$dir$proces, "annot_for_paired_data.RData"))


################################################################################
##We perform the GSA analysis of the transferred miRNA
# 
# ## ANNOTATION
# load (file.path (.job$dir$proces, "annot_for_paired_data.RData"))
# ls ()

ontologias <- names (annot)
ontologias

# ## INDEX
# load (file.path (.job$dir$proces, "rindex_paired.RData"))
# ls ()
# 
# tags <- names (rindex)
# TAGS <- toupper (tags)
# names (TAGS) <- tags
# tags
# TAGS

################################################################################

length (tags)

## tags <- tags[1:3]    ## 1
## tags <- tags[4:6]    ## 2
## tags <- tags[7:9]    ## 3
## tags <- tags[10:12]  ## 4
## tags <- tags[13:15]  ## 5
## tags <- tags[16:17]  ## 6

################################################################################


## GSA
setwd (file.path (.job$dir$proces, "res_uvgsa_paired"))


for (tag in tags) {
  for (onto in ontologias) {
    cat ("\n=============== ", TAGS[tag], ":", onto, " ===============\n")
    
    anotacion <- annotFilter (annot[[onto]], index = rindex[[tag]], minBlockSize = 10, maxBlockSize = 300)
    
    if (.job$testmode) anotacion <- anotacion[1:3]
    
    res <- try (uvGsa (rindex[[tag]], anotacion, fulltable = TRUE))
    
    fichero <- paste (tag, "_", onto, ".RData", sep = "")
    save (list = "res", file = fichero)
  }
}


################################################################################

## We collect all GSA results
options (width = 170)
corte <- 0.05
ontologias <- c ("bp", "cc", "mf")

## INDEX
tags <- names (rindex)
TAGS <- toupper (tags)   
names (TAGS) <- tags
tags
TAGS

################################################################################

## GSA list format
setwd (file.path (.job$dir$proces, "res_uvgsa_paired"))
dir ()

res.gsa.pair <- list ()
for (tag in tags) {
  cat ("\n", tag, fill = TRUE, sep = "")
  cat ("====================================================================================================", fill = TRUE)
  for (onto in ontologias) {
    fichero <- paste (tag, "_", onto, ".RData", sep = "")
    load (fichero)
    res[,"pat"] <- uvPat (res, cutoff = corte)
    res[,"index"] <- pval2index (pval = res[,"pval"], sign = res[,"lor"])
    res[,"Name"] <- getGOnames (res, verbose = FALSE)
    res.gsa.pair[[onto]][[tag]] <- res
  }
}

res[1:3,]
res.gsa.pair[[onto]][[tag]][1:3,]


## SOME CHECKINGS
sapply (res.gsa.pair, length)

for (onto in ontologias) {  
  cat ("\n")
  print (onto)
  print (sapply (res.gsa.pair[[onto]], dim))  ##OK all equal dimensions
  print (unique (sapply (res.gsa.pair[[onto]], function (x) sum (rownames (x) == rownames (res.gsa.pair[[onto]][[1]]))))) ##OK all equal rownames
  print (unique (sapply (res.gsa.pair[[onto]], function (x) sum (x[,"N"] == res.gsa.pair[[onto]][[1]][,"N"])))) ##OK all equal go sizes
  print (unique (sapply (res.gsa.pair[[onto]], function (x) sum (x[,"conv"])))) ## OK all converged
}


################################################################################

### Matrix Format

mat.gsa.pair <- list ()

for (onto in ontologias) {
  N <- res.gsa.pair[[onto]][[1]][,"N"]
  names (N) <- rownames (res.gsa.pair[[onto]][[1]])
  mat.gsa.pair[[onto]][["N"]] <- N
  
  mat.gsa.pair[[onto]][["lor"]]   <- sapply (res.gsa.pair[[onto]], function (x) {s <- x[,"lor"];   names (s) <- rownames (x); s})
  mat.gsa.pair[[onto]][["pval"]]  <- sapply (res.gsa.pair[[onto]], function (x) {s <- x[,"pval"];  names (s) <- rownames (x); s})
  mat.gsa.pair[[onto]][["padj"]]  <- sapply (res.gsa.pair[[onto]], function (x) {s <- x[,"padj"];  names (s) <- rownames (x); s})
  mat.gsa.pair[[onto]][["pat"]]   <- sapply (res.gsa.pair[[onto]], function (x) {s <- x[,"pat"];   names (s) <- rownames (x); s})
  mat.gsa.pair[[onto]][["index"]] <- sapply (res.gsa.pair[[onto]], function (x) {s <- x[,"index"]; names (s) <- rownames (x); s})
}

mat.gsa.pair[["tags"]] <- tolower (names (res.gsa.pair[[onto]]))
mat.gsa.pair[["TAGS"]] <- toupper (names (res.gsa.pair[[onto]]))
# mat.gsa.pair[["desc"]] <- tags.desc[mat.gsa.pair[["tags"]]]


###SAVE
save (list = c ("res.gsa.pair", "mat.gsa.pair"), file = file.path (.job$dir$proces, paste("res_gsa_paired", nsim, ".RData", 
                                                                                          sep = "")))



###EXIT
warnings ()
sessionInfo ()
q ("no")
