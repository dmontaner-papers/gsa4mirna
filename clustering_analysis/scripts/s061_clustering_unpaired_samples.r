########################################################################################
## s061_clustering_unpaired_samples.r
## 2015-08-25 fgarcia@cipf.es
## Clustering analysis for functional results in UNPAIRED samples
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
library (knitr); packageDescription ("knitr", fields = "Version") #"1.6"
library(mdgsa); packageDescription ("mdgsa", fields = "Version")

library(parallel); packageDescription ("parallel", fields = "Version") #"3.2.1"
library(pvclust); packageDescription ("pvclust", fields = "Version")# "1.3-2"

source ("function_arbol_2.r")
source ("function_pcaGenes_2.r")
fun0 <- function (x) {
  y <- unlist (strsplit (as.vector(x), split = "paired_")) [2]
  z <- unlist (strsplit (as.vector(y), split = "_")) [1]
  return (z)
}



### B. DATA
########################################################################################

setwd (file.path (.job$dir$res))
lista <- dir()

listau     <- lista[grep("res_gsa_unpaired", lista)]
listau_bp  <- listau[grep("_bp", listau)]
listau_mf  <- listau[grep("_mf", listau)]
listau_cc  <- listau[grep("_cc", listau)]
length(listau_bp); length(listau_cc); length(listau_mf); length(listau) 


bp <- read.xlsx(listau_bp[1], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
head(bp); dim(bp)  #3976
u_bp <- array(0, dim =c(nrow(bp), length(listau_bp)))
rownames(u_bp) <- rownames(bp)
colnames(u_bp) <- as.vector(sapply (listau_bp, fun0))
for (i in 1:length(listau_bp)){
  print(listau_bp[i])
  datos <- read.xlsx(listau_bp[i], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
  datos[,"rindex"]<- pval2index(pval = datos[,"pval"], sign = datos[, "lor"])
  datos <- datos[rownames(u_bp),]
  print(table(row.names(datos)==row.names(u_bp)))
  u_bp[,i] <- datos[,"rindex"]
}


cc <- read.xlsx(listau_cc[1], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
head(cc); dim(cc)  #446
u_cc <- array(0, dim =c(nrow(cc), length(listau_cc)))
rownames(u_cc) <- rownames(cc)
colnames(u_cc) <- as.vector(sapply (listau_cc, fun0))
for (i in 1:length(listau_cc)){
  print(listau_cc[i])
  datos <- read.xlsx(listau_cc[i], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
  datos[,"rindex"]<- pval2index(pval = datos[,"pval"], sign = datos[, "lor"])
  datos <- datos[rownames(u_cc),]
  print(table(row.names(datos)==row.names(u_cc)))
  u_cc[,i] <- datos[,"rindex"]
}


mf <- read.xlsx(listau_mf[1], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
head(mf); dim(mf)  #747
u_mf <- array(0, dim =c(nrow(mf), length(listau_mf)))
rownames(u_mf) <- rownames(mf)
colnames(u_mf) <- as.vector(sapply (listau_mf, fun0))
for (i in 1:length(listau_mf)){
  print(listau_mf[i])
  datos <- read.xlsx(listau_mf[i], row.names = 1,stringsAsFactors = FALSE,sheetIndex = 1)
  datos[,"rindex"]<- pval2index(pval = datos[,"pval"], sign = datos[, "lor"])
  datos <- datos[rownames(u_mf),]
  print(table(row.names(datos)==row.names(u_mf)))
  u_mf[,i] <- datos[,"rindex"]
}



#Convert colnames from lowercase to uppercase
colnames(u_bp) <- toupper(colnames(u_bp))
colnames(u_cc) <- toupper(colnames(u_cc))
colnames(u_mf) <- toupper(colnames(u_mf))


summary(u_bp)
summary(u_cc)
summary(u_mf)




### C. CLUSTERING
########################################################################################

setwd (file.path (.job$dir$plots))
x.por <- 1
y.por <- 1

### cluster with correlation distance, BP
correlacion <- cor (u_bp)
distancia <- as.dist ((1 - correlacion) / 2)
hc <- hclust (distancia)
png (filename ="cluster_corelationd_bp_unpaired.png",   width = 480 * x.por, height = 480 * y.por)
arbol (cluster = hc, main = "Clustering. Correlation distance. Unpaired. BP ")
dev.off ()
### cluster with euclidean distance, BP
distancia <- dist (t (u_bp), method = "euclidean") #dist trabaja por filas
he <- hclust (distancia)
png (filename = "cluster_euclideand_bp_unpaired.png", width = 480 * x.por, height = 480 * y.por)
arbol (cluster = he, main = "Clustering. Euclidean  distance. Unpaired.BP")
dev.off ()


### cluster with correlation distance, CC
correlacion <- cor (u_cc)
distancia <- as.dist ((1 - correlacion) / 2)
hc <- hclust (distancia)
png (filename ="cluster_corelationd_cc_unpaired.png",   width = 480 * x.por, height = 480 * y.por)
arbol (cluster = hc, main = "Clustering. Correlation distance. CC ")
dev.off ()
### cluster with euclidean distance, CC 
distancia <- dist (t (u_cc), method = "euclidean") #dist trabaja por filas
he <- hclust (distancia)
png (filename = "cluster_euclideand_cc_unpaired.png", width = 480 * x.por, height = 480 * y.por)
arbol (cluster = he, main = "Clustering. Euclidean  distance. Unpaired. CC")
dev.off ()


### cluster with correlation distance, MF
correlacion <- cor (u_mf)
distancia <- as.dist ((1 - correlacion) / 2)
hc <- hclust (distancia)
png (filename ="cluster_corelationd_mf_unpaired.png",   width = 480 * x.por, height = 480 * y.por)
arbol (cluster = hc, main = "Clustering. Correlation distance. Unpaired. MF ")
dev.off ()
### cluster with euclidean distance, MF 
distancia <- dist (t (u_mf), method = "euclidean") #dist trabaja por filas
he <- hclust (distancia)
png (filename = "cluster_euclideand_mf_unpaired.png", width = 480 * x.por, height = 480 * y.por)
arbol (cluster = he, main = "Clustering. Euclidean  distance. Unpaired.MF")
dev.off ()



### SIGNIFICANT cluster with correlation distance, BP 
setwd (file.path (.job$dir$plots))
cl<- makeCluster(2, type = "PSOCK")
x.por <- 4; y.por <- 2
mydat <- u_bp
trans.pv <- parPvclust(cl, mydat, nboot=10000)
## highlight clusters with high au p-values
png (filename = "sigcluster_corelationd_bp_unpaired.png", width = 480 * x.por, height = 480 * y.por, res = 200)
plot(trans.pv, main = "Clustering. Correlation distance. Unpaired.BP")
pvrect(trans.pv)
dev.off ()

### SIGNIFICANT cluster with correlation distance, CCP 
setwd (file.path (.job$dir$plots))
cl<- makeCluster(2, type = "PSOCK")
x.por <- 4; y.por <- 2
mydat <- u_cc
trans.pv <- parPvclust(cl, mydat, nboot=10000)
## highlight clusters with high au p-values
png (filename = "sigcluster_corelationd_cc_unpaired.png", width = 480 * x.por, height = 480 * y.por, res = 200)
plot(trans.pv, main = "Clustering. Correlation distance. Unpaired. CC")
pvrect(trans.pv)
dev.off ()

### SIGNIFICANT cluster with correlation distance, MF 
setwd (file.path (.job$dir$plots))
cl<- makeCluster(2, type = "PSOCK")
x.por <- 4; y.por <- 2
mydat <- u_mf
trans.pv <- parPvclust(cl, mydat, nboot=10000)
## highlight clusters with high au p-values
png (filename = "sigcluster_corelationd_mf_unpaired.png", width = 480 * x.por, height = 480 * y.por, res = 200)
plot(trans.pv, main = "Clustering. Correlation distance. Unpaired.MF")
pvrect(trans.pv)
dev.off ()





### D. PCAPLOT
########################################################################################

setwd (file.path (.job$dir$plots))
x.por <- 2
y.por <- 1

### BP
mi.pca <- pcaGenes (u_bp)
names (mi.pca)
sapply (mi.pca, dim)
png (filename = "pca_bp_unpaired.png",     width = 480 * x.por, height = 480 * y.por)
plot.pca.genes (mi.pca, main = "PCA plot. Unpaired. BP", pch = ".", cex=1.2)
dev.off ()

### CC
mi.pca <- pcaGenes (u_cc)
names (mi.pca)
sapply (mi.pca, dim)
png (filename = "pca_cc_unpaired.png",     width = 480 * x.por, height = 480 * y.por)
plot.pca.genes (mi.pca, main = "PCA plot. Unpaired. CC", pch = ".", cex=1.2)
dev.off ()

### MF
mi.pca <- pcaGenes (u_mf)
names (mi.pca)
sapply (mi.pca, dim)
png (filename = "pca_mf_unpaired.png",     width = 480 * x.por, height = 480 * y.por)
plot.pca.genes (mi.pca, main = "PCA plot. Unpaired. MF", pch = ".", cex=1.2)
dev.off ()







###SAVE
save(u_bp, u_cc, u_mf, file = file.path (.job$dir$proces,"clustering_unpaired_data.RData"))

###EXIT
warnings ()
sessionInfo ()
q ("no")
