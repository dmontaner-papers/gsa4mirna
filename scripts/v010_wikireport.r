date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.0.0 RC (2013-03-26 r62418)"
library (Biobase); packageDescription ("Biobase", fields = "Version") #"2.22.0"
#library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.3"
library (TiddlyWikiR); packageDescription ("TiddlyWikiR", fields = "Version") #"1.0.4"
help (package = TiddlyWikiR)

source (".job.r"); .job

### DATA
load (file.path (.job$dir$proces, "report", "allStudyInfo.RData"))  #this object includes all necessary data to prepare tables
ls ()
datos[1:3,]
dim(datos)

t (t (colnames (datos)))


################################################################################
### WIKI REPORT: START the wiki
################################################################################
infile  <- file.path (.job$dir$docs, "index_template.html")
outfile <- file.path (.job$dir$code, "index.html")
infile
outfile

writeTiddlers (infile = infile, outfile = outfile)


################################################################################
### WIKI REPORT: REPLACE TAGS
################################################################################
mistags <- list ()

# mistags[["@@~COMPLETION_DATE@@"]] <- as.character (format (Sys.Date(), format="%d %B %Y"))
# 
# mistags[["@@~NUMBER_OF_SAMPLES@@"]] <- as.character ('17')
# 
# mistags[["@@~WORK_AUTHORS@@"]] <- twList (c("''Francisco Garcia''", "''David Montaner''"))

# my_files <- dir('/media/FGG/tiddly/gsa4mirna/supplementary_files/files')
# my_plots <- dir('/home/paco/Desktop/plots')

my_files <- dir(file.path (.job$dir$code, "supplementary_files", "files"))
my_plots <- dir(file.path (.job$dir$code, "supplementary_files", "plots"))


## 0. Table with full information for all tumors
################################################################################

tabla.cols <- c ("ID", "miRNA.Down.unpa", "miRNA.noDif.unpa", "miRNA.Up.unpa",
                 "miRNA.Down.pair", "miRNA.noDif.pair", "miRNA.Up.pair",
                 "go.Down.unpa", "go.noDif.unpa", "go.UP.unpa",
                 "go.Down.pair", "go.noDif.pair", "go.UP.pair",
                 "bp.Down.unpa",  "bp.noDif.unpa", "bp.UP.unpa",
                 "cc.Down.unpa",  "cc.noDif.unpa", "cc.UP.unpa",
                 "mf.Down.unpa",  "mf.noDif.unpa", "mf.UP.unpa",
                 "bp.Down.pair",  "bp.noDif.pair", "bp.UP.pair",
                 "cc.Down.pair",  "cc.noDif.pair", "cc.UP.pair",
                 "mf.Down.pair",  "mf.noDif.pair", "mf.UP.pair") 

tabla.nams <- c ("ID", "Down.unpaired",  "noDif.unpaired",  "Up.unpaired", 
                 "Down.paired",  "noDif.paired", "Up.paired",
                 "gos.Derg.unpaired",  "gos.noDif.unpaired",  "gos.Inh.unpaired",
                 "gos.Derg.paired",    "gos.noDif.paired",    "gos.Inh.paired",
                 "bp.Down.unpa",  "bp.noDif.unpa", "bp.UP.unpa",
                 "cc.Down.unpa",  "cc.noDif.unpa", "cc.UP.unpa",
                 "mf.Down.unpa",  "mf.noDif.unpa", "mf.UP.unpa",
                 "bp.Down.pair",  "bp.noDif.pair", "bp.UP.pair",
                 "cc.Down.pair",  "cc.noDif.pair", "cc.UP.pair",
                 "mf.Down.pair",  "mf.noDif.pair", "mf.UP.pair") 

tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams


tabla[,"resPaired"]   <- paste0 ("res_edger_paired_",   rownames (tabla), ".xlsx")
tabla[,"resUnpaired"] <- paste0 ("res_edger_unpaired_", rownames (tabla), ".xlsx")

tabla[,"resPaired_gsa_cc"]   <- paste0 ("res_gsa_paired_",     rownames (tabla), "_cc.xlsx")
tabla[,"resPaired_gsa_mf"]   <- paste0 ("res_gsa_paired_",     rownames (tabla), "_mf.xlsx")
tabla[,"resPaired_gsa_bp"]   <- paste0 ("res_gsa_paired_",     rownames (tabla), "_bp.xlsx")
tabla[,"resUnpaired_gsa_cc"] <- paste0 ("res_gsa_unpaired_",   rownames (tabla), "_cc.xlsx")
tabla[,"resUnpaired_gsa_mf"] <- paste0 ("res_gsa_unpaired_",   rownames (tabla), "_mf.xlsx")
tabla[,"resUnpaired_gsa_bp"] <- paste0 ("res_gsa_unpaired_",   rownames (tabla), "_bp.xlsx")

##
links <- tabla  
links[] <- NA  ## vaciado   #the same structure to TABLA but full of NA
links[,"resPaired"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired"])
links[,"resUnpaired"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired"])

links[,"resPaired_gsa_cc"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired_gsa_cc"])
links[,"resUnpaired_gsa_cc"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired_gsa_cc"])
links[,"resPaired_gsa_mf"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired_gsa_mf"])
links[,"resUnpaired_gsa_mf"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired_gsa_mf"])
links[,"resPaired_gsa_bp"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired_gsa_bp"])
links[,"resUnpaired_gsa_bp"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired_gsa_bp"])



links <- as.matrix (links)
#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
#For this case, the tag is: @@~TABLA_all@@
mistags[["@@~TABLA_all@@"]] <-  twTable (dat = tabla, ref = links, sortable = TRUE, align = c(rep("c",39)))
writeTags (mistags, file = outfile)


#save this general table:
tabla0 <- tabla
links0 <- links



## 1. Analyzed datasets
################################################################################

tabla.cols <- c ("ID", "N.total" , "N.cases", "N.controls", "N.paired.samples", "cancer.type")
tabla.nams <- c ("ID", "total",  "cases",  "controls",  "paired",  "description")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

#adding links to TCGA:
tcga_links <- read.csv(file.path (.job$dir$docs, "tcga_links.csv"), sep = "\t", row.names = 1)
tcga_links <- tcga_links[rownames(tabla),]
tcga_links[,"mlinks"] <- paste0("[[", tcga_links[,"description"], "|", tcga_links[,"tcga_links"], "]]")
head(tcga_links) 
tabla1 <- cbind(tabla[, c("ID", "total", "cases", "controls", "paired")], tcga_links[,"mlinks"])
colnames(tabla1) <-  c("ID", "total", "cases", "controls", "paired", "description")
head(tabla1)


#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
#For this case, the tag is: @@~TABLA1_DATASETS@@
mistags[["@@~TABLA1_DATASETS@@"]] <-  twTable (dat = tabla1, sortable = TRUE, align = c("r", rep("c",4),"l"))
writeTags (mistags, file = outfile)




## 2. MiRNA level
################################################################################

tabla.cols <- c ("ID", "miRNA.Down.unpa", "miRNA.noDif.unpa", "miRNA.Up.unpa", "miRNA.Down.pair", "miRNA.noDif.pair", "miRNA.Up.pair")
tabla.nams <- c ("ID", "Down.unpaired",  "noDif.unpaired",  "Up.unpaired",  "Down.paired",  "noDif.paired", "Up.paired")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

tabla[,"resPaired"]   <- paste0 ("res_edger_paired_",   rownames (tabla), ".xlsx")
tabla[,"resUnpaired"] <- paste0 ("res_edger_unpaired_", rownames (tabla), ".xlsx")
##
links <- tabla  
links[] <- NA  ## vaciado   #the same structure to TABLA but full of NA
links[,"resPaired"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired"])
links[,"resUnpaired"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired"])
links <- as.matrix (links)
#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
#For this case, the tag is: @@~TABLA_MIRNAS@@
mistags[["@@~TABLA2_MIRNAS@@"]] <-  twTable (dat = tabla, ref = links, sortable = TRUE, align = c(rep("c",9)))
writeTags (mistags, file = outfile)




# ############################################################################################## 
# ## MiRNA level
# 
# my_files2 <- my_files[grep("res_edger_paired", my_files)]
# my_files2
# links <- rep('NA', length(my_files2))
# data <- cbind(my_files2, links)
# class(data)
# data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
# colnames(data) <- c("files", "links")
# mistags[["@@~TABLA_MIRNA_LEVEL_PAIRED@@"]] <- twTable (dat = data, sortable = TRUE)
# 
# my_files2 <- my_files[grep("res_edger_unpaired", my_files)]
# my_files2
# links <- rep('NA', length(my_files2))
# data <- cbind(my_files2, links)
# class(data)
# data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
# colnames(data) <- c("files", "links")
# mistags[["@@~TABLA_MIRNA_LEVEL_UNPAIRED@@"]] <- twTable (dat = data, sortable = TRUE)
# 


## 3. GENE level
################################################################################

# tabla 3
######################

tabla.cols <- c ("ID", "targets.Down.unpa", "targets.Intersect.unpa", "targets.Up.unpa", "targets.Down.pair",
                 "targets.Intersect.pair", "targets.Up.pair")
tabla.nams <- c ("ID", "Down.unpaired",  "Inter.unpaired",  "Up.unpaired",  "Down.paired",  "Inter.paired", "Up.paired")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

tabla[,"transfer_index_paired"]   <- paste0 ("transfer_index_paired_",   rownames (tabla), ".xlsx")
tabla[,"transfer_index_unpaired"] <- paste0 ("transfer_index_unpaired_", rownames (tabla), ".xlsx")
##
##
links <- tabla  
links[] <- NA  ## vaciado   #the same structure to TABLA but full of NA
links[,"transfer_index_paired"]   <- paste0 ("supplementary_files/files/", tabla[,"transfer_index_paired"])
links[,"transfer_index_unpaired"] <- paste0 ("supplementary_files/files/", tabla[,"transfer_index_unpaired"])
links <- as.matrix (links)

mistags[["@@~TABLA3_GENES@@"]] <-  twTable (dat = tabla, ref = links, sortable = TRUE, align = c(rep("c",9)))
writeTags (mistags, file = outfile)




# and more information: 
######################

# my_files2 <- my_files[grep("transfer_index_paired", my_files)]
# my_files2
# links <- rep('NA', length(my_files2))
# data <- cbind(my_files2, links)
# class(data)
# data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
# colnames(data) <- c("files", "links")
# mistags[["@@~TABLA_GENE_LEVEL_TRANSFER_PAIRED@@"]] <- twTable (dat = data, sortable = TRUE)
# 
# my_files2 <- my_files[grep("transfer_index_unpaired", my_files)]
# my_files2
# links <- rep('NA', length(my_files2))
# data <- cbind(my_files2, links)
# class(data)
# data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
# colnames(data) <- c("files", "links")
# mistags[["@@~TABLA_GENE_LEVEL_TRANSFER_UNPAIRED@@"]] <- twTable (dat = data, sortable = TRUE)
# 

my_plots2 <- my_plots[grep("paired_explore_transfer", my_plots)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
# mistags[["@@~PLOTS_GENE_LEVEL@@"]] <- twTable (dat = plots, sortable = TRUE)
 mistags[["@@~PLOTS_GENE_LEVEL@@"]] <- plots[,"links"]

# TO CHECK:
# [img(80%+,)[Boxplot of the data BEFORE normalization|results/plots/raw_boxplot.png]]
# my_plots2 <- my_plots[grep("paired_explore_transfer", my_plots)]
# plots.links <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")



## 4. Gene Set 
############################################################################################## 


## TABLE 4
######################

tabla.cols <- c ("ID", "targetGOs.Down.unpa", "targetGOs.Intersect.unpa", "targetGOs.Up.unpa", "targetGOs.Down.pair",
                 "targetGOs.Intersect.pair", "targetGOs.Up.pair")
tabla.nams <- c ("ID", "Down.unpaired",  "Inter.unpaired",  "Up.unpaired",  "Down.paired",  "Inter.paired", "Up.paired")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
mistags[["@@~TABLA4_GOS@@"]] <-  twTable (dat = tabla,   sortable = TRUE, align = c(rep("c",7)))
writeTags (mistags, file = outfile)







## TABLE 5
######################

tabla <- tabla0[, c ("ID", "gos.Derg.unpaired",  "gos.noDif.unpaired",  "gos.Inh.unpaired", 
                           "gos.Derg.paired", "gos.noDif.paired", "gos.Inh.paired", 
                     "bp.Down.unpa",  "bp.noDif.unpa", "bp.UP.unpa",
                     "cc.Down.unpa",  "cc.noDif.unpa", "cc.UP.unpa",
                     "mf.Down.unpa",  "mf.noDif.unpa", "mf.UP.unpa",
                     "bp.Down.pair",  "bp.noDif.pair", "bp.UP.pair",
                     "cc.Down.pair",  "cc.noDif.pair", "cc.UP.pair",
                     "mf.Down.pair",  "mf.noDif.pair", "mf.UP.pair",
                     "resPaired_gsa_cc", "resUnpaired_gsa_cc", "resPaired_gsa_mf",
                     "resUnpaired_gsa_mf", "resPaired_gsa_bp",   "resUnpaired_gsa_bp")] 

links <- links0[, c ("ID", "gos.Derg.unpaired",  "gos.noDif.unpaired",  "gos.Inh.unpaired", 
                     "gos.Derg.paired", "gos.noDif.paired", "gos.Inh.paired", 
                     "bp.Down.unpa",  "bp.noDif.unpa", "bp.UP.unpa",
                     "cc.Down.unpa",  "cc.noDif.unpa", "cc.UP.unpa",
                     "mf.Down.unpa",  "mf.noDif.unpa", "mf.UP.unpa",
                     "bp.Down.pair",  "bp.noDif.pair", "bp.UP.pair",
                     "cc.Down.pair",  "cc.noDif.pair", "cc.UP.pair",
                     "mf.Down.pair",  "mf.noDif.pair", "mf.UP.pair",
                     "resPaired_gsa_cc", "resUnpaired_gsa_cc", "resPaired_gsa_mf",
                     "resUnpaired_gsa_mf", "resPaired_gsa_bp",   "resUnpaired_gsa_bp")]

mistags[["@@~TABLA5_GOS@@"]] <-  twTable (dat = tabla,  ref = links, sortable = TRUE, align = c(rep("c",31)))
writeTags (mistags, file = outfile)
                                       




## Gene Set level
######################

# my_files2 <- my_files[grep("res_gsa_paired", my_files)]
# my_files2
# links <- rep('NA', length(my_files2))
# data <- cbind(my_files2, links)
# class(data)
# data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
# colnames(data) <- c("files", "links")
# mistags[["@@~TABLA_GENE_SET_LEVEL_GO_PAIRED@@"]] <- twTable (dat = data, sortable = TRUE)
# 
# 
# my_files2 <- my_files[grep("res_gsa_unpaired", my_files)]
# my_files2
# links <- rep('NA', length(my_files2))
# data <- cbind(my_files2, links)
# class(data)
# data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
# colnames(data) <- c("files", "links")
# mistags[["@@~TABLA_GENE_SET_LEVEL_GO_UNPAIRED@@"]] <- twTable (dat = data, sortable = TRUE)

my_files2 <- c("common_enrichment_paired.xlsx", "common_enrichment_unpaired.xlsx" )
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLE_GENE_SET_LEVEL_ALL_STUDIES@@"]] <- twTable (dat = data, sortable = TRUE)


my_plots2 <- my_plots[grep("paired_size_effect", my_plots)]  #paired and unpaired
my_plots2 <- my_plots2[grep("un", my_plots2, invert= T)]     #only paired
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
# mistags[["@@~PLOTS_GENE_SET_LEVEL_SIZE_EFFECT_PAIRED@@"]] <- twTable (dat = plots, sortable = TRUE)
mistags[["@@~PLOTS_GENE_SET_LEVEL_SIZE_EFFECT_PAIRED@@"]] <- plots[,"links"]


my_plots2 <- my_plots[grep("unpaired_size_effect", my_plots)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
# mistags[["@@~PLOTS_GENE_SET_LEVEL_SIZE_EFFECT_UNPAIRED@@"]] <- twTable (dat = plots, sortable = TRUE)
mistags[["@@~PLOTS_GENE_SET_LEVEL_SIZE_EFFECT_UNPAIRED@@"]] <- plots[,"links"]





## Global results
######################

# Distribution of ranking index 
my_plots2 <- "paired_rindex_boxplot.png"
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
# mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_RANKING_INDEX_PAIRED@@"]] <- twTable (dat = plots, sortable = TRUE)
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_RANKING_INDEX_PAIRED@@"]] <- plots[,"links"]

my_plots2 <- "unpaired_rindex_boxplot.png"
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
# mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_RANKING_INDEX_UNPAIRED@@"]] <- twTable (dat = plots, sortable = TRUE)
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_RANKING_INDEX_UNPAIRED@@"]] <- plots[,"links"]



# Ranking index correlation
my_plots2 <- c("paired_cor_rindex0.png", "paired_cor_rindexN.png", "paired_cor_rindexT.png")
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_RANKING_INDEX_CORRELATION_PAIRED@@"]] <- plots[,"links"]

my_plots2 <- c("unpaired_cor_rindex0.png", "unpaired_cor_rindexN.png", "unpaired_cor_rindexT.png")
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_RANKING_INDEX_CORRELATION_UNPAIRED@@"]] <- plots[,"links"]


# Inhibition effect correlation distribution by ontology
my_plots2 <- my_plots[grep("_dist_of_", my_plots)]
my_plots2 <- my_plots[grep("un", my_plots2, invert = T)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_BY_ONTOLOGY_PAIRED@@"]] <- plots[,"links"]

my_plots2 <- my_plots[grep("_dist_of_", my_plots)]
my_plots2 <- my_plots2[grep("un", my_plots2)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_BY_ONTOLOGY_UNPAIRED@@"]] <-  plots[,"links"]



# Ranking index correlation vs. correlation
my_plots2 <- "paired_rindex_cor_vs_cor.png"
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_RANKING_INDEX_CORRELATION_VS_CORRELATION_PAIRED@@"]] <-  plots[,"links"]

my_plots2 <- "unpaired_rindex_cor_vs_cor.png"
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_RANKING_INDEX_CORRELATION_VS_CORRELATION_UNPAIRED@@"]] <- plots[,"links"]



############################################################################################## 
## by cancer

canceres <- c("blca", "brca", "cesc", "coad", "esca","hnsc", "kick", "kirc", "kirp", "lihc", "luad", "lusc", "paad", "pcpg", "prad",
              "read", "skcm", "stad", "thca", "ucec")


tabla <- tabla0[, c ("ID", "Down.unpaired",  "noDif.unpaired",  "Up.unpaired",
                     "Down.paired",  "noDif.paired", "Up.paired",
                     "gos.Derg.unpaired",  "gos.noDif.unpaired",  "gos.Inh.unpaired", 
                     "gos.Derg.paired", "gos.noDif.paired", "gos.Inh.paired", 
                     "bp.Down.unpa",  "bp.noDif.unpa", "bp.UP.unpa",
                     "cc.Down.unpa",  "cc.noDif.unpa", "cc.UP.unpa",
                     "mf.Down.unpa",  "mf.noDif.unpa", "mf.UP.unpa",
                     "bp.Down.pair",  "bp.noDif.pair", "bp.UP.pair",
                     "cc.Down.pair",  "cc.noDif.pair", "cc.UP.pair",
                     "mf.Down.pair",  "mf.noDif.pair", "mf.UP.pair")] 

tablai <- tabla[, c ("ID", "Down.unpaired",  "noDif.unpaired",  "Up.unpaired",
                     "gos.Derg.unpaired",  "gos.noDif.unpaired",  "gos.Inh.unpaired", 
                     "bp.Down.unpa",  "bp.noDif.unpa", "bp.UP.unpa",
                     "cc.Down.unpa",  "cc.noDif.unpa", "cc.UP.unpa",
                     "mf.Down.unpa",  "mf.noDif.unpa", "mf.UP.unpa")] 
rownames(tablai) <- paste0(rownames(tablai), "_unpaired")
colnames(tablai) <- c ("ID", "miRNA.down",  "miRNA.noDif",  "miRNA.up",
                       "gos.derg",  "gos.noDif",  "gos.inh", 
                       "bp.Down",  "bp.noDif", "bp.UP",
                       "cc.Down",  "cc.noDif", "cc.UP",
                       "mf.Down",  "mf.noDif", "mf.UP") 



tablap <- tabla0[, c ("ID", "Down.paired",  "noDif.paired", "Up.paired",
                      "gos.Derg.paired", "gos.noDif.paired", "gos.Inh.paired", 
                      "bp.Down.pair",  "bp.noDif.pair", "bp.UP.pair",
                      "cc.Down.pair",  "cc.noDif.pair", "cc.UP.pair",
                      "mf.Down.pair",  "mf.noDif.pair", "mf.UP.pair")] 
rownames(tablap) <- paste0(rownames(tablap), "_paired")
colnames(tablap) <- c ("ID", "miRNA.down",  "miRNA.noDif",  "miRNA.up",
                       "gos.derg",  "gos.noDif",  "gos.inh", 
                       "bp.Down",  "bp.noDif", "bp.UP",
                       "cc.Down",  "cc.noDif", "cc.UP",
                       "mf.Down",  "mf.noDif", "mf.UP") 

tablaf <- rbind(tablai,tablap)




  
for (i in 1:length(canceres)) {
  
  tabla <- tablaf[grep(canceres[i], rownames(tablaf)),]
  links <- as.data.frame(links0)[canceres[i],]
  mistags[[paste("@@~SUMMARY_", canceres[i], "@@", sep ="")]] <-  twTable (dat = tabla,   sortable = TRUE,
                                                                           align = c(rep("c",16)))
  
  my_plots2 <- my_plots[grep(canceres[i], my_plots)]
  links <- rep('NA', length(my_plots2))
  plots <- cbind(my_plots2, links)
  class(plots)
  #plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
  #mistags[[paste("@@~PLOTS_", canceres[i], "@@", sep ="")]] <- twTable (dat = plots, sortable = TRUE)

  plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
  mistags[[paste("@@~PLOTS_", canceres[i], "@@", sep ="")]] <- plots[,"links"]
  
  my_files2 <- my_files[grep(canceres[i], my_files)]
  my_files2
  links <- rep('NA', length(my_files2))
  data <- cbind(my_files2, links)
  class(data)
  data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
  colnames(data) <- c("files", "links")
  mistags[[paste("@@~FILES_", canceres[i], "@@", sep ="")]] <- twTable (dat = data, sortable = TRUE)
  
}             

writeTags (mistags, file = outfile)











# #TO CHECK:
# 
# 
# tabla1<-twTable (dat = tabla, ref = links, sortable = TRUE, title = "TABLE2")
# mistiddlers <- list ()
# mistiddlers[["dos"]] <- newTiddler (title = "BLCA2", content = list ("!dfañsdfja", tabla1, "sdafsdf", tabla1))
# writeTiddlers (mistiddlers, outfile = outfile)
# 





###EXIT
warnings ()
sessionInfo ()
q ("no")
