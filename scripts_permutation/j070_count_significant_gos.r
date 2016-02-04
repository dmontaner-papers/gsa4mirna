##j70_count_significant_gos.r
##2014-12-04 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script counts the number of significant GO terms in all cancers
############################################################


date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.2 (2015-06-18)"

fper5  <- function(x) {quantile(x, probs = 0.05)}
fper95 <- function(x) {quantile(x, probs = 0.95)}

try (source (".job.r")); try (.job)
options (width = 170)
corte <- 0.05




## PAIRED DATA
############################################################

setwd(file.path (.job$dir$proces))
lista <- dir()
pair    <- lista [grep ("res_gsa_paired", lista)]
length(pair)


# measuring total number of GO terms
load (file.path (.job$dir$proces, pair[1]))
n_bp  <- nrow(res.gsa.pair$bp[[1]])
n_mf  <- nrow(res.gsa.pair$mf[[1]])
n_cc  <- nrow(res.gsa.pair$cc[[1]])
n_gos <- n_bp + n_mf + n_cc
n_gos

# defining a matrix to save results
length(names(res.gsa.pair$bp))   #17
pair.tum   <- names(res.gsa.pair$bp)

mat.pair <- matrix(data =NA, nrow = length(names(res.gsa.pair$bp)), ncol= length(pair))
rownames(mat.pair) <-names(res.gsa.pair$bp)
colnames(mat.pair) <- c(1:100)
mat.pair[1:5,1:15]

for (i in 1:length(pair)) {
  load (file.path (.job$dir$proces, pair[i]))
  for (j in pair.tum){
    sig.bp <- sum(as.numeric(res.gsa.pair$bp[[j]][,"padj"] < corte))
    sig.cc <- sum(as.numeric(res.gsa.pair$cc[[j]][,"padj"] < corte))
    sig.mf <- sum(as.numeric(res.gsa.pair$mf[[j]][,"padj"] < corte))
    mat.pair[j,i] <-  sig.bp + sig.cc + sig.mf
  }
}
mat.pair[1:5,1:15]


#  % of significant GOs 
mat.pair.per <- (mat.pair / n_gos) * 100

summary(mat.pair.per)
sim.median            <- round(apply(mat.pair.per, 1, median),3) 
sim.percentile5       <- round(apply(mat.pair.per, 1, fper5),3)
sim.percentile95      <- round(apply(mat.pair.per, 1, fper95),3)  
Cancer                <- toupper(rownames(mat.pair.per))
sim.res.pair <- cbind(Cancer, sim.median, sim.percentile5,  sim.percentile95) 
colnames(sim.res.pair) <- c("Cancer", "median", "perc5", "perc95")
sim.res.pair 





  
## UNPAIRED DATA
############################################################


setwd(file.path (.job$dir$proces))
lista <- dir()
unpair    <- lista [grep ("res_gsa_unpaired", lista)]
length(unpair)


# measuring total number of GO terms
load (file.path (.job$dir$proces, unpair[1]))
n_bp  <- nrow(res.gsa.unpa$bp[[1]])
n_mf  <- nrow(res.gsa.unpa$mf[[1]])
n_cc  <- nrow(res.gsa.unpa$cc[[1]])
n_gos <- n_bp + n_mf + n_cc
n_gos


# defining a matrix to save results
length(names(res.gsa.unpa$bp))   
unpair.tum   <- names(res.gsa.unpa$bp)

mat.unpair <- matrix(data =NA, nrow = length(names(res.gsa.unpa$bp)), ncol= length(unpair))
rownames(mat.unpair) <-names(res.gsa.unpa$bp)
colnames(mat.unpair) <- c(1:100)
mat.unpair[1:5,1:15]

for (i in 1:length(unpair)) {
  load (file.path (.job$dir$proces, unpair[i]))
  for (j in unpair.tum){
    sig.bp <- sum(as.numeric(res.gsa.unpa$bp[[j]][,"padj"] < corte))
    sig.cc <- sum(as.numeric(res.gsa.unpa$cc[[j]][,"padj"] < corte))
    sig.mf <- sum(as.numeric(res.gsa.unpa$mf[[j]][,"padj"] < corte))
    mat.unpair[j,i] <-  sig.bp + sig.cc + sig.mf
  }
}

#  % of significant GOs
mat.unpair.per <- (mat.unpair / n_gos) * 100

summary(mat.unpair.per)
sim.median            <- round(apply(mat.unpair.per, 1, median),3) 
sim.percentile5       <- round(apply(mat.unpair.per, 1, fper5),3)
sim.percentile95      <- round(apply(mat.unpair.per, 1, fper95),3)  
Cancer                <- toupper(rownames(mat.unpair.per))
sim.res.unpair <- cbind(Cancer, sim.median, sim.percentile5, sim.percentile95) 
colnames(sim.res.unpair) <- c("Cancer", "median", "perc5", "perc95")
sim.res.unpair 



### save results
save (list = c ("sim.res.pair", "sim.res.unpair"),
      file = file.path (.job$dir$proces, "res_gsa_conteos.RData"))

write.table(sim.res.unpair, file = file.path (.job$dir$proces, "res_gsa_conteos_unpaired.txt"), 
            quote = F, sep = "\t", row.names = F)
write.table(sim.res.pair, file = file.path (.job$dir$proces, "res_gsa_conteos_paired.txt"), 
            quote = F, sep = "\t", row.names = F)

###EXIT
warnings ()
sessionInfo ()
q ("no")
