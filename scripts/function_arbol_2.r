##function_arbol_2.r
##2010-03-10 dmontaner@cipf.es
##function to plot a clustering tree
##ADDING A COLOR AXIS

##to be called after "hclust"

## $clase ELEMENT has to be introduced in the resutl of "hclus"

###EXAMPLE
## distancia <- dist (t (datos), method = "euclidean") #dist trabaja por filas
## he <- hclust (distancia)
## he$clase <- batch  ##COLOREAMOS SEGUN LA VARIABLE batch
## arbol (cluster = he, main = "euclidean distance")

arbol <- function (cluster, ...) {#2010-03-10 dmontaner@cipf.es
  ##
  plot (cluster, hang = 0.1, axes = FALSE, ann = FALSE)
  title (...)
  ##
  ##color
  clase <- as.character (cluster$clase)
  clase.unica <- unique (clase)
  clr <- rainbow (length (clase.unica))
  names (clr) <- clase.unica
  ##
  posiciones <- 1:length (cluster$labels)
  labels.ordenados <- cluster$labels[cluster$order]
  clase.ordenada <- clase[cluster$order]
  ##
  for (cls in clase.unica) {
    touse <- clase.ordenada %in% cls
    axis (1, posiciones[touse], labels = labels.ordenados[touse],
          las = 3, col.axis = clr[cls])
  }
  ##
  #try (legend (x = 3, y = max (cluster$height), legend = names (clr), fill = clr))
  #try (legend (x = "topright", legend = names (clr), fill = clr))
}
