## .job.r
## 2010-03-09 dmontaner@cipf.es
## 2014-05-01 dmontaner@cipf.es
## script that keeps the settings for each job
## bear in mind that it is dependent on the data folder structure

## .NAME: name of the people for whom the analysis is done.
##        It will be used as the JOB_NAME (no spaces; always lower case, etc...)

## dir$data: path to the directory where all data (raw, generated and results) are stored.
##          The last directory in the path has to be JOB_NAME.

## dir$code: path to the directory where scripts, sample information and documentation are stored.
##          The last directory in the path has to be JOB_NAME.

## NOTE: datadir & codedir may be different or the SAME.

################################################################################

.NAME = "gsa4mirna" #job name. NO SPACES

################################################################################

.job <- list ()
.job$name <- .NAME
.job$dir <- list ()

### rootDir: Location in MY COMPUTER
.job$dir$data <- file.path ("~", "datos",         "2014", .job$name) #starting with ("~") if working in your home directory
.job$dir$code <- file.path ("~", "trabajos_mios", "2014", .job$name) #or ("") if working in the root directory !!!

## .job$dir$data <- file.path ("~", "Desktop", "papers", .job$name, "datos") #starting with ("~") if working in your home directory
## .job$dir$code <- file.path ("~", "Desktop", "papers", .job$name) #or ("") if working in the root directory !!!

### MORE directories
.job$dir$scripts <- file.path (.job$dir$code, "scripts")
.job$dir$docs    <- file.path (.job$dir$code, "documents")

.job$dir$rawdat      <- file.path (.job$dir$data, "data_raw")
.job$dir$annotation  <- file.path (.job$dir$data, "data_annotation")
.job$dir$proces      <- file.path (.job$dir$data, "data_processed")

.job$dir$plots  <- file.path (.job$dir$data, "results", "plots")
.job$dir$res    <- file.path (.job$dir$data, "results", "files")

### Some other parameters
.job$testmode <- FALSE ##testing mode
.job$dec <- "."
.job$idsep <- " /// "  ##separates IDs

################################################################################

rm (list = ".NAME")

##MESSAGE
cat ("\n.job.r has been sourced\n", fill = TRUE)
