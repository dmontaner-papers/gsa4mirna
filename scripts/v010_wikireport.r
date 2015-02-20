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
#infile  <- file.path ("/home","paco", "Desktop", "mirnas.html")
infile  <- file.path ("/home","paco", "Desktop", "index_template.html")
outfile <- file.path ("/home","paco", "Desktop", "index.html")

writeTiddlers (infile = infile, outfile = outfile)


################################################################################
### WIKI REPORT: REPLACE TAGS
################################################################################
mistags <- list ()

mistags[["@@~COMPLETION_DATE@@"]] <- as.character (format (Sys.Date(), format="%d %B %Y"))

mistags[["@@~NUMBER_OF_SAMPLES@@"]] <- as.character ('17')

mistags[["@@~WORK_AUTHORS@@"]] <- twList (c("''Francisco Garcia''", "''David Montaner''"))

my_files <- dir('/media/FGG/tiddly/gsa4mirna/supplementary_files/files')
my_plots <- dir('/home/paco/Desktop/plots')
 
##

pru <- rep('NA', length(my_files))
links <- pru
data <- cbind(my_files, links)
class(data)
data[,"links"] <- paste ("[[",file.path (paste(my_files,'/media/FGG/tiddly/gsa4mirna/supplementary_files/files',sep='|'), my_files),"]]",sep="")
mistags[["@@~FILES_TABLE@@"]] <- twTable (dat = data, sortable = TRUE)
# 
pru <- rep('NA', length(my_plots))
links <- pru
plots <- cbind(my_plots, links)
class(plots)
plots[,"links"] <- paste ("[img[","/home/paco/Desktop/plots/", my_plots,"]]",sep="")
mistags[["@@~PLOTS_TABLE@@"]] <- twTable (dat = plots, sortable = TRUE)

writeTags (mistags, file = outfile)



###EXIT
warnings ()
sessionInfo ()
q ("no")
