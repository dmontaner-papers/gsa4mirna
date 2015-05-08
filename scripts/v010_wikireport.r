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

my_files <- dir("/home/paco/Desktop/papers/gsa4mirna/supplementary_files/files")
my_plots <- dir("/home/paco/Desktop/papers/gsa4mirna/supplementary_files/plots")



## 1. Analyzed datasets
################################################################################

tabla.cols <- c ("ID", "N.total" , "N.cases", "N.controls", "N.paired.samples", "cancer.type")
tabla.nams <- c ("ID", "total",  "cases",  "controls",  "paired",  "description")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
#For this case, the tag is: @@~TABLA1_DATASETS@@
mistags[["@@~TABLA1_DATASETS@@"]] <-  twTable (dat = tabla, sortable = TRUE, align = c("r", rep("c",4),"l"))
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

# tabla[,"resPaired"]   <- paste0 ("res_edger_paired_",   rownames (tabla), ".xlsx")
# tabla[,"resUnpaired"] <- paste0 ("res_edger_unpaired_", rownames (tabla), ".xlsx")
# ##
# links <- tabla  
# links[] <- NA  ## vaciado   #the same structure to TABLA but full of NA
# links[,"resPaired"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired"])
# links[,"resUnpaired"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired"])
# links <- as.matrix (links)

#NOTE: this table doesn't contain xls files for target genes
#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
#For this case, the tag is: @@~TABLA_GENES@@
mistags[["@@~TABLA3_GENES@@"]] <-  twTable (dat = tabla, sortable = TRUE, align = c(rep("c",7)))
writeTags (mistags, file = outfile)




# and more information: 
######################

my_files2 <- my_files[grep("transfer_index_paired", my_files)]
my_files2
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLA_GENE_LEVEL_TRANSFER_PAIRED@@"]] <- twTable (dat = data, sortable = TRUE)

my_files2 <- my_files[grep("transfer_index_unpaired", my_files)]
my_files2
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLA_GENE_LEVEL_TRANSFER_UNPAIRED@@"]] <- twTable (dat = data, sortable = TRUE)

my_plots2 <- my_plots[grep("paired_explore_transfer", my_plots)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_LEVEL@@"]] <- twTable (dat = plots, sortable = TRUE)





## 4. Gene Set 
############################################################################################## 


## TABLE 4
######################

tabla.cols <- c ("ID", "targetsGOs.Down.unpa", "targetsGOs.Intersect.unpa", "targetsGOs.Up.unpa", "targetsGOs.Down.pair",
                 "targetsGOs.Intersect.pair", "targetsGOs.Up.pair")
tabla.nams <- c ("ID", "Down.unpaired",  "Inter.unpaired",  "Up.unpaired",  "Down.paired",  "Inter.paired", "Up.paired")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.
mistags[["@@~TABLA4_GOS@@"]] <-  twTable (dat = tabla, sortable = TRUE, align = c(rep("c",7)))
writeTags (mistags, file = outfile)



## TABLE 5
######################

tabla.cols <- c ("ID", "go.Down.unpa", "go.noDif.unpa", "go.UP.unpa", 
                       "go.Down.pair", "go.noDif.pair", "go.UP.pair")
tabla.nams <- c ("ID", "Derg.unpaired",  "noDif.unpaired",  "Inh.unpaired",  "Derg.paired",  "noDif.paired", "Inh.paired")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams

# tabla[,"resPaired"]   <- paste0 ("res_edger_paired_",   rownames (tabla), ".xlsx")
# tabla[,"resUnpaired"] <- paste0 ("res_edger_unpaired_", rownames (tabla), ".xlsx")
# ##
# links <- tabla  
# links[] <- NA  ## vaciado   #the same structure to TABLA but full of NA
# links[,"resPaired"]   <- paste0 ("supplementary_files/files/", tabla[,"resPaired"])
# links[,"resUnpaired"] <- paste0 ("supplementary_files/files/", tabla[,"resUnpaired"])
# links <- as.matrix (links)

#After generating the table, I have to wikify it. Good clue: in the template I have to indicate some tag to change this tag by the new table.

# mistags[["@@~TABLA5_GOS@@"]] <-  twTable (dat = tabla, ref = links, sortable = TRUE, align = c(rep("c",7)), 
mistags[["@@~TABLA5_GOS@@"]] <-  twTable (dat = tabla,  sortable = TRUE, align = c(rep("c",7)))
writeTags (mistags, file = outfile)
                                       




## Gene Set level
######################

my_files2 <- my_files[grep("res_gsa_paired", my_files)]
my_files2
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLA_GENE_SET_LEVEL_GO_PAIRED@@"]] <- twTable (dat = data, sortable = TRUE)


my_files2 <- my_files[grep("res_gsa_unpaired", my_files)]
my_files2
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLA_GENE_SET_LEVEL_GO_UNPAIRED@@"]] <- twTable (dat = data, sortable = TRUE)


my_plots2 <- my_plots[grep("paired_size_effect", my_plots)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_SIZE_EFFECT_PAIRED@@"]] <- twTable (dat = plots, sortable = TRUE)


my_plots2 <- my_plots[grep("unpaired_size_effect", my_plots)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_SIZE_EFFECT_UNPAIRED@@"]] <- twTable (dat = plots, sortable = TRUE)



## Gene Set level
######################

my_files2 <- c("allStudyInfo.xlsx", "common_enrichment_paired.xlsx", "common_enrichment_unpaired.xlsx" )
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLE_GENE_SET_LEVEL_ALL_STUDIES@@"]] <- twTable (dat = data, sortable = TRUE)


my_plots2 <- c("paired_rindex_boxplot.png", "unpaired_rindex_boxplot.png")
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_RANKING_INDEX@@"]] <- twTable (dat = plots, sortable = TRUE)


my_plots2 <- c("paired_rindex_boxplot.png", "unpaired_rindex_boxplot.png")
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_RANKING_INDEX@@"]] <- twTable (dat = plots, sortable = TRUE)



my_plots2 <- c("paired_cor_rindex0.png", "paired_cor_rindexN.png", "paired_cor_rindexT.png", 
  "unpaired_cor_rindex0.png", "unpaired_cor_rindexN.png", "unpaired_cor_rindexT.png")
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_RANKING_INDEX_CORRELATION@@"]] <- twTable (dat = plots, sortable = TRUE)


my_plots2 <- my_plots[grep("_dist_of_", my_plots)]
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_DISTRIBUTION_BY_ONTOLOGY@@"]] <- twTable (dat = plots, sortable = TRUE)



my_plots2 <- c("paired_rindex_cor_vs_cor.png", "unpaired_rindex_cor_vs_cor.png")
links <- rep('NA', length(my_plots2))
plots <- cbind(my_plots2, links)
class(plots)
plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
mistags[["@@~PLOTS_GENE_SET_LEVEL_RANKING_INDEX_CORRELATION_VS_CORRELATION@@"]] <- twTable (dat = plots, sortable = TRUE)




############################################################################################## 
## by cancer

canceres <- c("blca", "brca", "cesc", "coad", "esca","hnsc", "kick", "kirc", "kirp", "lihc", "luad", "lusc", "paad", "pcpg", "prad",
              "read", "skcm", "stad", "thca", "ucec")
  
for (i in 1:length(canceres)) {
  my_plots2 <- my_plots[grep(canceres[i], my_plots)]
  links <- rep('NA', length(my_plots2))
  plots <- cbind(my_plots2, links)
  class(plots)
  plots[,"links"] <- paste ("[img[","supplementary_files/plots/", my_plots2,"]]",sep="")
  mistags[[paste("@@~PLOTS_", canceres[i], "@@", sep ="")]] <- twTable (dat = plots, sortable = TRUE)

  
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











#TO CHECK:


tabla1<-twTable (dat = tabla, ref = links, sortable = TRUE, title = "TABLE2")
mistiddlers <- list ()
mistiddlers[["dos"]] <- newTiddler (title = "BLCA2", content = list ("!dfañsdfja", tabla1, "sdafsdf", tabla1))
writeTiddlers (mistiddlers, outfile = outfile)






###EXIT
warnings ()
sessionInfo ()
q ("no")
