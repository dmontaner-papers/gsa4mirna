##function_signifCount_1.r
##2014-01-27 dmontaner@cipf.es
##General purpose function to count significant ressults

signifCount <- function (mylist, stat = "statistic", pval = "padj", cutoff = 0.05,
                         colNames = c("N. down-regulated", "N. not-diff.", "N. up-regulated"), verbose = FALSE) {##2014-01-27 dmontaner@cipf.es
    
    if (!is.list (mylist))        stop ("mylist should be a named list")
    if (is.null (names (mylist))) stop ("no names found for mylist")
    
    salida <- matrix (NA, nrow = length (mylist), ncol = 3)
    rownames (salida) <- names (mylist)
    colnames (salida) <- colNames
    ##
    for (el in names (mylist)) {
        li <- mylist[[el]]
        signif <- sign (li[,stat]) * (li[,pval] < cutoff)
        salida[el,] <- c (sum (signif == -1), sum (signif == 0), sum (signif == 1))
    }
    
    salida <- as.data.frame (salida, stringsAsFactors = FALSE)
    return (salida)
}

