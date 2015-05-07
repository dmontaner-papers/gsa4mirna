##t010_paper_tables.r
##2014-12-12 dmontaner@cipf.es
##Analysis of the miRNA data from The Cancer Genome Atlas
##This script prints some tables for the report

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.1.0 (2014-04-10)"
#library (xlsx); packageDescription ("xlsx", fields = "Version") #"0.5.5"
library (knitr); packageDescription ("knitr", fields = "Version") #"1.6"
library (xtable); packageDescription ("xtable", fields = "Version") #"1.7-3"
library (tables); packageDescription ("tables", fields = "Version") #"0.7.79"
#help (package = xtable)
#help (package = tables)

try (source (".job.r")); try (.job)

options (width = 170)

### DATA
load (file.path (.job$dir$proces, "report", "allStudyInfo.RData"))
ls ()
datos[1:3,]

t (t (colnames (datos)))

placement <- "hbt"


################################################################################
## Table 1: overwiew
################################################################################

tabla.file <-     "table1"
tabla.labe <- "tab:table1"
tabla.cols <- c ("ID", "N.total", "N.cases", "N.controls", "N.paired.samples", "cancer.type")
tabla.nams <- c ("ID", "total"  , "cases"   , "controls",   "paired",          "description")
tabla.capt <- c ("Analyzed datasets. Columns of the table display:
TCGA disease ID,
the total number of samples in the analysis,
the number of tumoral samples,
the number of control samples (solid normal tissue),
the number of paired samples available in the dataset
and the cancer type.")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams
tabla
##
writeLines (kable (tabla, caption = tabla.capt, format = "pandoc", row.names = FALSE), 
            con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".md")))
##
xt <- xtable (tabla, caption = tabla.capt, label = tabla.labe)
align (xt)[2] <- paste0 ("@{}", align (xt)[2])
print (xt,
       include.rownames = FALSE, table.placement = placement,
       size = "scriptsize",
       tabular.environment = "tabular*",
       width = "\\columnwidth",
       file = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))

## kable (tabla, format = "markdown", caption = tabla.capt)
## kable (tabla, format = "pandoc",   caption = tabla.capt)
## kable (tabla, format = "rst",      caption = tabla.capt)
## kable (tabla, format = "latex",    caption = tabla.capt)
## kable (tabla, format = "html",     caption = tabla.capt)


################################################################################
## Table 2: ressults at miRNA level
################################################################################

tabla.file <-     "table2"
tabla.labe <- "tab:table2"
tabla.cols <- c ("ID", "miRNA.Down.unpa", "miRNA.noDif.unpa", "miRNA.Up.unpa", "miRNA.Down.pair", "miRNA.noDif.pair", "miRNA.Up.pair")
tabla.nams <- c ("ID",       "Down",            "noDif",            "Up",            "Down",            "noDif",            "Up")
##tabla.nams <- tabla.cols
multicolumns <- "& \\\\multicolumn{3}{c}{Unpaired} & \\\\multicolumn{3}{c}{Paired} \\\\\\\\ \\\\cmidrule(r){2-4} \\\\cmidrule(r){5-7} \\\\\\\\"
##my.add.to.row <- list (pos = list (-1), command = multicolumns)
tabla.capt <- c ("Number of up, down and not differentially regulated miRNAS in each cancer type.") 
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams
tabla
##
writeLines (kable (tabla, caption = tabla.capt, format = "pandoc", row.names = FALSE), 
            con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".md")))
##
xt <- xtable (tabla, caption = tabla.capt, label = tabla.labe)
align (xt)[2] <- paste0 ("@{\\extracolsep{\\fill}}", align (xt)[2])
tableLines <- print (xt,
                     include.rownames = FALSE, table.placement = placement,
                     size = "footnotesize",
                     tabular.environment = "tabular*",
                     width = "\\columnwidth",
                     booktabs = TRUE)
                     ## add.to.row = my.add.to.row,
                     ## file = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))
tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
writeLines (tableLines, con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))


################################################################################
## Table 3: ressults at gene level
################################################################################

tabla.file <-     "table3"
tabla.labe <- "tab:table3"
tabla.cols <- c ("ID", "targets.Down.unpa", "targets.Intersect.unpa", "targets.Up.unpa", "targets.Down.pair", "targets.Intersect.pair", "targets.Up.pair")
tabla.nams <- c ("ID",         "Down",               "Inter",                 "Up",              "Down",               "Inter",                 "Up")
##tabla.nams <- tabla.cols
multicolumns <- "& \\\\multicolumn{3}{c}{Unpaired} & \\\\multicolumn{3}{c}{Paired} \\\\\\\\ \\\\cmidrule(r){2-4} \\\\cmidrule(r){5-7} \\\\\\\\"
##my.add.to.row <- list (pos = list (-1), command = multicolumns)
tabla.capt <- c ("Number of genes targeted by the up and down regulated miRNAS.
The intersection column (Inter) shows the number of genes which are targets of both, the up and down regulated miRNAs.")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams
tabla
##
writeLines (kable (tabla, caption = tabla.capt, format = "pandoc", row.names = FALSE), 
            con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".md")))
##
xt <- xtable (tabla, caption = tabla.capt, label = tabla.labe)
align (xt)[2] <- paste0 ("@{\\extracolsep{\\fill}}", align (xt)[2])
tableLines <- print (xt,
                     include.rownames = FALSE, table.placement = placement,
                     size = "footnotesize",
                     tabular.environment = "tabular*",
                     width = "\\columnwidth",
                     booktabs = TRUE)
                     ## add.to.row = my.add.to.row,
                     ## file = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))
tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
writeLines (tableLines, con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))


################################################################################
## Table 4: intersecction of the go terms
################################################################################

tabla.file <-     "table4"
tabla.labe <- "tab:table4"
tabla.cols <- c ("ID", "targetGOs.Down.unpa", "targetGOs.Intersect.unpa", "targetGOs.Up.unpa", "targetGOs.Down.pair", "targetGOs.Intersect.pair", "targetGOs.Up.pair")
tabla.nams <- c ("ID",           "Down",                "Inter",                    "Up",                "Down",                "Inter",                    "Up")
#tabla.nams <- tabla.cols
multicolumns <- "& \\\\multicolumn{3}{c}{Unpaired} & \\\\multicolumn{3}{c}{Paired} \\\\\\\\ \\\\cmidrule(r){2-4} \\\\cmidrule(r){5-7} \\\\\\\\"
##my.add.to.row <- list (pos = list (-1), command = multicolumns)
tabla.capt <- c ("Number of GO terms associated to the genes targeted by the up and down regulated miRNAs.
Most GO terms are targeted in cases and controls at the same time
as it can be seen in the intersection column (Inter).")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams
tabla
##
writeLines (kable (tabla, caption = tabla.capt, format = "pandoc", row.names = FALSE), 
            con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".md")))
##
xt <- xtable (tabla, caption = tabla.capt, label = tabla.labe)
align (xt)[2] <- paste0 ("@{\\extracolsep{\\fill}}", align (xt)[2])
tableLines <- print (xt,
                     include.rownames = FALSE, table.placement = placement,
                     size = "footnotesize",
                     tabular.environment = "tabular*",
                     width = "\\columnwidth",
                     booktabs = TRUE)
                     ## add.to.row = my.add.to.row,
                     ## file = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))
tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
writeLines (tableLines, con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))


################################################################################
## Table 5: ressults at gene set level
################################################################################

tabla.file <-     "table5"
tabla.labe <- "tab:table5"
tabla.cols <- c ("ID", "go.Down.unpa", "go.noDif.unpa", "go.UP.unpa", "go.Down.pair", "go.noDif.pair", "go.UP.pair")
tabla.nams <- c ("ID",    "Derg.",        "noDif",         "Inh.",       "Derg.",        "noDif",         "Inh.")
#tabla.nams <- tabla.cols
multicolumns <- "& \\\\multicolumn{3}{c}{Unpaired} & \\\\multicolumn{3}{c}{Paired} \\\\\\\\ \\\\cmidrule(r){2-4} \\\\cmidrule(r){5-7} \\\\\\\\"
tabla.capt <- c ("Number significant GO terms in the functional profiling analysis for the paired and unpaired comparisons. 
Columns \\textbf{Inh.} indicates the number of terms with a \\textbf{positive} $\\alpha$ coefficient in the logistic regression analysis.
Those are the terms inhibited or intercepted in cases.
Columns \\textbf{Derg.} indicates the number of terms with a \\textbf{negative} $\\alpha$ value.
Those are the terms inhibited in controls or \\emph{deregulated} in cases.
Columns noDif indicate the number of GOs with a not significant slope coefficient.")
tabla <- datos[,tabla.cols]
colnames (tabla) <- tabla.nams
tabla
##
writeLines (kable (tabla, caption = tabla.capt, format = "pandoc", row.names = FALSE), 
            con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".md")))
##
xt <- xtable (tabla, caption = tabla.capt, label = tabla.labe)
align (xt)[2] <- paste0 ("@{\\extracolsep{\\fill}}", align (xt)[2])
tableLines <- print (xt,
                     include.rownames = FALSE, table.placement = placement,
                     size = "footnotesize",
                     tabular.environment = "tabular*",
                     width = "\\columnwidth",
                     booktabs = TRUE)
                     ## add.to.row = my.add.to.row,
## file = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))
tableLines <- sub ("\\\\toprule\\n", paste0 ("\\\\toprule\n", multicolumns, "\n"), tableLines) ## booktabs = TRUE
tableLines <- sub ("\\\\hline\\n",   paste0 ("\\\\hline\n",   multicolumns, "\n"), tableLines) ## booktabs = FALSE
writeLines (tableLines, con = file.path (.job$dir$code, "paper", "tables", paste0 (tabla.file, ".tex")))



## datos[1:3,]
## t (t (colnames (datos)))
## colnames (datos)



###EXIT
warnings ()
sessionInfo ()
#q ("no")
