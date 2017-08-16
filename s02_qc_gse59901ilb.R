#########################################################################################
##                                                                                     ##
## s02.5_gse59901ilb.R                                                                 ##
## Con este script comprobamos la calidad de los datos crudos y obtenemos un plot.     ##
## Recordar que los datos de metilación son buenos cuando el log de la mediana de las  ##
## intensidades de ambos canales para cada muestra es muy similar, el grado en el que  ##
## pueden divergir lo marca minfi para cada experimento dependiendo de los niveles de  ##
## intensidad.                                                                         ##
##                                                                                     ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

# Cargamos las librerías necesarias

library(minfi)

setwd("/home/atrassierra/GSE59901") # Vamos al directorio donde se encuentran los datos

data <- readGEORawFile("tableilb.txt", sep = "\t", Uname = "Signal_A", Mname = "Signal_B", row.names = 1, array = "IlluminaHumanMethylation450k")

qc <- getQC(data) # Obtenemos la calidad de las muestras

head(qc) # Podemos observarlas numéricamente

plot <- plotQC(qc) # Ploteamos las calidades

# Guardamos el plot

dev.copy(png, 'QCgse59901ilb.png')
dev.off()
