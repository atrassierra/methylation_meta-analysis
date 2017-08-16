#########################################################################################
##                                                                                     ##
## s03_gse59901idc.R                                                                   ##
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

data <- readGEORawFile("tableidc_norm.txt", sep = "\t", Uname = "Signal_A", Mname = "Signal_B", row.names = 1, array = "IlluminaHumanMethylation450k")

grs <- ratioConvert(data) # Calculamos el estadístico beta con el que se computan las dmrs

g <- as.factor(c(rep("c", 4), rep("idc", 10))) # Especificamos los controles y los casos

matriz <- model.matrix(~ g) # Construimos la matriz de diseño experimental

# Bumphunter

dmrs <- bumphunter(grs, design = matriz, cutoff = 0, B= 50, type="Beta")

setwd("/home/atrassierra/GSE59901")

# Guardamos los resultados obtenidos

save(dmrs, file = "dmrs59901idc.RData")
