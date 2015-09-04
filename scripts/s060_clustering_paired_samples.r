########################################################################################
## s060_clustering_paired_samples.r
## 2015-08-25 fgarcia@cipf.es
## Clustering analysis for functional results in PAIRED samples
## We use this transformation: sign(LOR) * -1  * log(raw p-value) from pval2index 
########################################################################################


date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
try (source (".job.r")); try (.job)
options (width = 170)



### A. LOADING LIBRARIES AND FUNCTIONS
########################################################################################

library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
# library (knitr); packageDescription ("knitr", fields = "Version") #"1.6"
library(mdgsa); packageDescription ("mdgsa", fields = "Version")

# library(parallel); packageDescription ("parallel", fields = "Version") #"3.2.1"
# library(pvclust); packageDescription ("pvclust", fields = "Version")# "1.3-2"

source("000_function_arbol_2.r")
source("000_function_pcaGenes_2.r")
fun0 <- function (x) {
  y <- unlist (strsplit (as.vector(x), split = "paired_")) [2]
  z <- unlist (strsplit (as.vector(y), split = "_")) [1]
  return (z)
}



### B. DATA
########################################################################################

setwd (file.path (.job$dir$res))
lista <- dir()

listap     <- lista[grep("res_gsa_paired", lista)]
listap_bp  <- listap[grep("_bp", listap)]
listap_mf  <- listap[grep("_mf", listap)]
listap_cc  <- listap[grep("_cc", listap)]
length(listap_bp); length(listap_cc); length(listap_mf); length(listap)


bp <- read.xlsx(listap_bp[1], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
head(bp); dim(bp)  #3976
p_bp <- array(0, dim =c(nrow(bp), length(listap_bp)))
rownames(p_bp) <- rownames(bp)
colnames(p_bp) <- as.vector(sapply (listap_bp, fun0))
for (i in 1:length(listap_bp)){
 print(listap_bp[i])
  datos <- read.xlsx(listap_bp[i], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
  datos[,"rindex"]<- pval2index(pval = datos[,"pval"], sign = datos[, "lor"])
  datos <- datos[rownames(p_bp),]
  print(table(row.names(datos)==row.names(p_bp)))
  p_bp[,i] <- datos[,"rindex"]
}



cc <- read.xlsx(listap_cc[1], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
head(cc); dim(cc)  #446
p_cc <- array(0, dim =c(nrow(cc), length(listap_cc)))
rownames(p_cc) <- rownames(cc)
colnames(p_cc) <- as.vector(sapply (listap_cc, fun0))
for (i in 1:length(listap_cc)){
  print(listap_cc[i])
  datos <- read.xlsx(listap_cc[i], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
  datos[,"rindex"]<- pval2index(pval = datos[,"pval"], sign = datos[, "lor"])
  datos <- datos[rownames(p_cc),]
  print(table(row.names(datos)==row.names(p_cc)))
  p_cc[,i] <- datos[,"rindex"]
}


mf <- read.xlsx(listap_mf[1], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
head(mf); dim(mf)  #747
p_mf <- array(0, dim =c(nrow(mf), length(listap_mf)))
rownames(p_mf) <- rownames(mf)
colnames(p_mf) <- as.vector(sapply (listap_mf, fun0))
for (i in 1:length(listap_mf)){
  print(listap_mf[i])
  datos <- read.xlsx(listap_mf[i], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
  datos[,"rindex"]<- pval2index(pval = datos[,"pval"], sign = datos[, "lor"])
  datos <- datos[rownames(p_mf),]
  print(table(row.names(datos)==row.names(p_mf)))
  p_mf[,i] <- datos[,"rindex"]
}




#Convert colnames from lowercase to uppercase
colnames(p_bp) <- toupper(colnames(p_bp))
colnames(p_cc) <- toupper(colnames(p_cc))
colnames(p_mf) <- toupper(colnames(p_mf))


summary(p_bp)
summary(p_cc)
summary(p_mf)




### C. CLUSTERING
########################################################################################

setwd (file.path (.job$dir$plots))
x.por <- 1
y.por <- 1

### cluster with correlation distance, BP
correlacion <- cor (p_bp)
distancia <- as.dist ((1 - correlacion) / 2)
hc <- hclust (distancia)
png (filename ="clust_paired_correlationd_bp.png",   width = 480 * x.por, height = 480 * y.por)
arbol (cluster = hc, main = "Clustering. Correlation distance. Paired. BP ")
dev.off ()
### cluster with euclidean distance, BP
distancia <- dist (t (p_bp), method = "euclidean") #dist trabaja por filas
he <- hclust (distancia)
png (filename = "clust_paired_euclideand_bp.png", width = 480 * x.por, height = 480 * y.por)
arbol (cluster = he, main = "Clustering. Euclidean  distance. Paired. BP")
dev.off ()


### cluster with correlation distance, CC
correlacion <- cor (p_cc)
distancia <- as.dist ((1 - correlacion) / 2)
hc <- hclust (distancia)
png (filename ="clust_paired_correlationd_cc.png",   width = 480 * x.por, height = 480 * y.por)
arbol (cluster = hc, main = "Clustering. Correlation distance. Paired. CC ")
dev.off ()
### cluster with euclidean distance, CC 
distancia <- dist (t (p_cc), method = "euclidean") #dist trabaja por filas
he <- hclust (distancia)
png (filename = "clust_paired_euclideand_cc.png", width = 480 * x.por, height = 480 * y.por)
arbol (cluster = he, main = "Clustering. Euclidean  distance. Paired.CC")
dev.off ()


### cluster with correlation distance, MF
correlacion <- cor (p_mf)
distancia <- as.dist ((1 - correlacion) / 2)
hc <- hclust (distancia)
png (filename ="clust_paired_correlationd_mf.png",   width = 480 * x.por, height = 480 * y.por)
arbol (cluster = hc, main = "Clustering. Correlation distance. Paired. MF ")
dev.off ()
### cluster with euclidean distance, MF 
distancia <- dist (t (p_mf), method = "euclidean") #dist trabaja por filas
he <- hclust (distancia)
png (filename = "clust_paired_euclideand_mf.png", width = 480 * x.por, height = 480 * y.por)
arbol (cluster = he, main = "Clustering. Euclidean  distance. Paired. MF")
dev.off ()



# ### SIGNIFICANT cluster with correlation distance, BP 
# setwd (file.path (.job$dir$plots))
# cl<- makeCluster(2, type = "PSOCK")
# x.por <- 4; y.por <- 2
# mydat <- p_bp
# trans.pv <- parPvclust(cl, mydat, nboot=10000)
# ## highlight clusters with high au p-values
# png (filename = "sigcluster_corelationd_bp_paired.png", width = 480 * x.por, height = 480 * y.por, res = 200)
# plot(trans.pv, main = "Clustering. Correlation distance. Paired. BP")
# pvrect(trans.pv)
# dev.off ()
# 
# ### SIGNIFICANT cluster with correlation distance, CCP 
# setwd (file.path (.job$dir$plots))
# cl<- makeCluster(2, type = "PSOCK")
# x.por <- 4; y.por <- 2
# mydat <- p_cc
# trans.pv <- parPvclust(cl, mydat, nboot=10000)
# ## highlight clusters with high au p-values
# png (filename = "sigcluster_corelationd_cc_paired.png", width = 480 * x.por, height = 480 * y.por, res = 200)
# plot(trans.pv, main = "Clustering. Correlation distance. Paired. CC")
# pvrect(trans.pv)
# dev.off ()
# 
# ### SIGNIFICANT cluster with correlation distance, MF 
# setwd (file.path (.job$dir$plots))
# cl<- makeCluster(2, type = "PSOCK")
# x.por <- 4; y.por <- 2
# mydat <- p_mf
# trans.pv <- parPvclust(cl, mydat, nboot=10000)
# ## highlight clusters with high au p-values
# png (filename = "sigcluster_corelationd_mf_paired.png", width = 480 * x.por, height = 480 * y.por, res = 200)
# plot(trans.pv, main = "Clustering. Correlation distance. Paired. MF")
# pvrect(trans.pv)
# dev.off ()
# 




### D. PCAPLOT
########################################################################################

setwd (file.path (.job$dir$plots))
x.por <- 2
y.por <- 1

### BP
mi.pca <- pcaGenes (p_bp)
names (mi.pca)
sapply (mi.pca, dim)
png (filename = "pca_bp_paired.png",     width = 480 * x.por, height = 480 * y.por)
plot.pca.genes (mi.pca, main = "PCA plot. Paired. BP", pch = ".", cex=1.2)
dev.off ()

### CC
mi.pca <- pcaGenes (p_cc)
names (mi.pca)
sapply (mi.pca, dim)
png (filename = "pca_cc_paired.png",     width = 480 * x.por, height = 480 * y.por)
plot.pca.genes (mi.pca, main = "PCA plot. Paired. CC", pch = ".", cex=1.2)
dev.off ()

### MF
mi.pca <- pcaGenes (p_mf)
names (mi.pca)
sapply (mi.pca, dim)
png (filename = "pca_mf_paired.png",     width = 480 * x.por, height = 480 * y.por)
plot.pca.genes (mi.pca, main = "PCA plot. Paired. MF", pch = ".", cex=1.2)
dev.off ()



###SAVE
save(p_bp, p_cc, p_mf, file = file.path (.job$dir$proces,"clustering_paired_data.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
