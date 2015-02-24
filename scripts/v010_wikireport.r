date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string #"R version 3.0.0 RC (2013-03-26 r62418)"
library (Biobase); packageDescription ("Biobase", fields = "Version") #"2.22.0"
#library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"0.3.3"
library (TiddlyWikiR); packageDescription ("TiddlyWikiR", fields = "Version") #"1.0.3"
help (package = TiddlyWikiR)



################################################################################
### WIKI REPORT: START the wiki
################################################################################
infile  <- file.path ("/home","paco", "Desktop", "papers", "gsa4mirna","documents", "index_template.html")
outfile <- file.path ("/home","paco", "Desktop", "papers", "gsa4mirna",  "index.html")
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



############################################################################################## 
## MiRNA level

my_files2 <- my_files[grep("res_edger_paired", my_files)]
my_files2
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLA_MIRNA_LEVEL_PAIRED@@"]] <- twTable (dat = data, sortable = TRUE)

my_files2 <- my_files[grep("res_edger_unpaired", my_files)]
my_files2
links <- rep('NA', length(my_files2))
data <- cbind(my_files2, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files2,"supplementary_files/files",sep='|'), my_files2),"]]",sep="")
colnames(data) <- c("files", "links")
mistags[["@@~TABLA_MIRNA_LEVEL_UNPAIRED@@"]] <- twTable (dat = data, sortable = TRUE)


############################################################################################## 
## Gene level

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






############################################################################################## 
## Gene Set level

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


############################################################################################## 
## Gene Set level

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



writeTags (mistags, file = outfile)



###EXIT
warnings ()
sessionInfo ()
q ("no")
