#########################################################################################
##                                                                                     ##
## s05_gsa_kegg.R                                                                      ##
## El presente script ejecuta varios GSA, uno para cada una de las comparaciones       ##
## estudiadas. Utilizamos las rutas disponibles en Kegg.                               ##
## Autores: Antonio Manuel Trassierra Fresco, Francisco García García                  ##
##                                                                                     ##
#########################################################################################


date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
rm (list = ls ())
R.version.string ##"R version 3.2.1

# Cargamos las librerías necesarias
library (mdgsa); packageDescription ("mdgsa", fields = "Version") #"1.2.0"

setwd("/home/atrassierra/meta_analysis/scripts/")
try (source (".job.r")); try (.job)



# A. Preparamos la anotación
#########################################################################################
setwd (file.path (.job$dir, "data","data_annotation"))
anot2 <- read.table("hsa_KEGG.txt", sep = "\t", header = T, as.is = T)
dim(anot2)
head(anot2)

# Algunos checkings
table(is.na(anot2[, 1]))
table(is.na(anot2[, 2]))

table(anot2[,1] == "")
table(anot2[,2] == "")

table(duplicated(anot2[, 1])) # Componentes únicos
table(duplicated(anot2[, 2])) # Keggs únicos

anot <- annotMat2list (anot2)
head (anot)
length (unique (anot2[,2]))




# B. GSEA
#########################################################################################
setwd (file.path (.job$dir, "results",  "difexp"))
ficheros <- dir()
ficheros <- ficheros[grep("RData", ficheros)]
ficheros

for (i in 1:length(ficheros)){
  setwd (file.path (.job$dir, "results",  "difexp"))
  cat ("\n ============================= \n")
  print(ficheros[i])
  load(ficheros[i])
  
  # Transformamos los genes en un ranking
  final <- final[!is.na(final$genes), ]
  rindex <- indexTransform (index = final$value, method = "normalize")
  names(rindex) <- as.vector(final$genes)
  # Transforma el ranking para que su distribucion sea adecuada como variable independiente de un modelo de regresion logistica
  summary(rindex)
  head(rindex)
  
  res <- uvGsa (index = rindex, annot = anot, fulltable = TRUE)
  print(table(res$padj < 0.05))
 
  # Realiza un analisis con un modelo de regresion logistica
  res[1:3,]
  save (list = "res", file = file.path (.job$dir, "results", "gsea", paste("kegg_", ficheros[i], sep = "")))
}
  
  
  
  
  
# Salida
warnings ()
sessionInfo ()
q ("no")



