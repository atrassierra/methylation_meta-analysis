#########################################################################################
##                                                                                     ##
## s02.5_normalization.r                                                               ##
## Script para normalizar y realizar de nuevo un análisis exploratorio de todos las    ##
## comparaciones presentes en este estudio.                                            ##
## Autores: Antonio Manuel Trassierra Fresco, Francisco García García                  ##
##                                                                                     ##
#########################################################################################


date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string 

# Cargamos las librerías necesarias

library (limma); packageDescription ("limma", fields = "Version") #"1.2.0"
library(Biobase)
try (source (".job.r")); try (.job)

# Funciones necesarias para llevar a cabo los árboles de clustering y los pca.

source (file.path (.job$dir, "scripts", "function_arbol_2.r"))
source (file.path (.job$dir, "scripts", "function_pcaGenes_2.r"))

# Seleccionamos todas las comparaciones para analizarlas.

setwd(file.path(.job$dir, "estudios"))
estudios  <- dir()[grep("RData", dir())]
estudios
nombres <- tolower(unlist(strsplit(estudios, split = ".RData")))
nombres



for (i in 1:length(estudios)){
  # Leemos los datos
  estudios[i]
  setwd(file.path(.job$dir, "estudios"))
  load(estudios[i])
  ls()  # Dos objetos: datos y g (grupos experimentales)

  
  rownames(datos) <- datos[,1]
  datos <- datos[,-1]
  datos <- as.matrix(datos)
  
  # Normalización por cuantiles
  
  datos <- normalizeBetweenArrays(datos, method = "quantile")
  head(datos)
  summary(datos)
  cat("ID\t", file= paste(nombres[i],"_norm.txt", sep =""))
  write.table(datos, paste(nombres[i],"_norm.txt", sep = ""), quote = F, sep = "\t", append = T)
  
  # Análisis exploratorio
  
  if (file.exists(file.path (.job$dir, "resultados_normalizados", "explore"))){
    setwd(file.path(.job$dir, "resultados_normalizados", "explore"))
  } else {
    dir.create(file.path(.job$dir, "resultados_normalizados", "explore"), recursive = T)
    setwd(file.path(.job$dir, "resultados_normalizados", "explore"))    
  }
  
  # Asignamos colores a los grupos
  char.g <- as.character (unlist(g))
  tan <- unique (char.g)
  micolor <- rainbow (length (tan))
  names (micolor) <- tan
  color.g <- micolor[char.g]
  
  
  # Boxplot
  ###########
  x.por <- 2.5 
  y.por <- 1.5
  png (filename = paste(nombres[i], "_boxplot_norm.png", sep  = ""),  width = 480 * x.por, height = 480 * y.por)
  boxplot(datos, cex= 0.1, col = color.g)
  try (legend (x=1, y= 0.11, legend = levels(unlist(g)), fill = unique(color.g)))
  dev.off ()
  
  
  
  # Clustering
  #############
  x.por <- 3
  y.por <- 1
  # Clúster con distancia de correlación
  correlacion <- cor (datos)
  distancia <- as.dist ((1 - correlacion) / 2)
  hc <- hclust (distancia)
  hc$clase <- unlist(g)  #variable color
  png (filename = paste(nombres[i],"_cluster_corelation_norm.png", sep = ""),   width = 480 * x.por, height = 480 * y.por)
  arbol (cluster = hc, main = "Clustering. Correlation distance ")
  dev.off ()
  
  # Clúster con distancia euclídea 
  distancia <- dist (t (datos), method = "euclidean") 
  he <- hclust (distancia)
  he$clase <- unlist(g)  # variable color
  png (filename = paste(nombres[i], "_cluster_euclideand_norm.png", sep = ""), width = 480 * x.por, height = 480 * y.por)
  arbol (cluster = he, main = "Clustering. Euclidean  distance")
  dev.off ()
  
  
  # PCA
  ###########
  x.por <- 2
  y.por <- 1
  
  mi.pca <- pcaGenes (datos)
  names (mi.pca)
  sapply (mi.pca, dim)
  
  png (filename = paste(nombres[i], "_pca_norm.png", sep= ""),   width = 480 * x.por, height = 480 * y.por)
  plot.pca.genes (mi.pca, main = "PCA plot", col = color.g, col.class = unlist(g), pch = ".")
  dev.off ()  
}




# Salida
warnings ()
sessionInfo ()
q ("no")



