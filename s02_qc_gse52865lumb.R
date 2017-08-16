#########################################################################################
##                                                                                     ##
## s02.5_gse52865lumb.R                                                                ##
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

setwd("/home/atrassierra/GSE52865") # Vamos al directorio donde se encuentran los datos

data <- readGEORawFile("gse52865lumb.txt", sep = "\t", Uname = "unmeth", Mname = "\\.meth", row.names = 1, array = "IlluminaHumanMethylation450k")

qc <- getQC(data) # Obtenemos la calidad de las muestras

head(qc) # Podemos observarlas numéricamente

plot <- plotQC(qc) # Ploteamos las calidades

# Guardamos el plot

dev.copy(png, 'QCgse52865lumb.png')
dev.off()