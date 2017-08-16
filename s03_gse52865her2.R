#########################################################################################
##                                                                                     ##
## s03_gse52865her2.R                                                                  ##
## Este script es donde se buscan las DMRs (regiones diferencialmente metiladas) entre ##
## las muestras. Para ello usaremos la función del paquete minfi, bumphunter.          ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

# Cargamos las librerías necesarias.

library(doParallel)
registerDoParallel(cores = 4)
library(minfi)

setwd("/home/atrassierra/meta_analysis/estudios") # Vamos al directorio donde se encuentran los datos

data <- readGEORawFile("gse52865her2_norm.txt", sep = "\t", Uname = "unmeth", Mname = "\\.meth", row.names = 1, array = "IlluminaHumanMethylation450k")

grs <- ratioConvert(data) # Calculamos el estadístico beta con el que se computan las dmrs

g <- as.factor(c(rep("c", 17), rep("her2", 3))) # Especificamos los controles y los casos

matriz <- model.matrix(~ g) # Construimos la matriz de diseño experimental

# Bumphunter

dmrs <- bumphunter(grs, design = matriz, cutoff = 0, B= 50, type="Beta")

setwd("/home/atrassierra/GSE52865")

# Guardamos los resultados obtenidos

save(dmrs, file = "dmrs52865her2.RData")
