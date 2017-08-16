#########################################################################################
##                                                                                     ##
## s01_gse78751.R                                                                      ##
## Este script descomprime el archivo donde se encuentra el estudio y lo procesa,      ##
## dejando las señales de metilación y no metilación para su posterior normalización.  ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

library(GEOquery)

setwd("/home/atrassierra/GSE78751") # Cambiamos el directorio de trabajo al del experimento en cuestión

gunzip('GSE78751_signal_intensities.txt.gz') # Descomprimimos el archivo para tener acceso al experimento.

# Leemos el archivo del experimento y extraemos las señales de metilación y no metilación
gse78751 = read.csv("GSE78751_signal_intensities.txt", sep = "\t")

targetID <- gse78751[, 2]

# Obtenemos las señales metilada y no metilada

metilado <- grep("\\.Meth", colnames(gse78751))

unmetil <- grep("Unmeth", colnames(gse78751))

combinado <- as.vector(rbind(metilado, unmetil))

gse78751signals <- gse78751[, combinado]

# Separamos el experimento para quedarnos con la comparación de cáncer de mama triple negativo vs control
idc <- gse78751signals[, c(1:46)]
  
control <- gse78751signals[, c(71:78)]

####################
##                ##
## idc vs control ##
##                ##
####################
gse78751idc <- cbind(targetID, control, idc) # Control vs Tumor

# Los nombres de las muestras están mal escritos y no coinciden, arreglamos esto escribiéndolos de nuevo
colnames(gse78751idc) <- c("targetID", "Control7.Meth", "Control7.Unmeth", "Control31.Meth", "Control31.Unmeth", "Control32.Meth", "Control32.Unmeth", "Control33.Meth", "Control33.Unmeth", "Paciente1.Meth", "Paciente1.Unmeth", "Paciente2.Meth", "Paciente2.Unmeth", "Paciente3.Meth", "Paciente3.Unmeth", "Paciente4.Meth", "Paciente4.Unmeth", "Paciente5.Meth", "Paciente5.Unmeth","Paciente6.Meth", "Paciente6.Unmeth","Paciente7.Meth", "Paciente7.Unmeth","Paciente8.Meth", "Paciente8.Unmeth","Paciente9.Meth", "Paciente9.Unmeth","Paciente10.Meth", "Paciente10.Unmeth","Paciente11.Meth", "Paciente11.Unmeth","Paciente12.Meth", "Paciente12.Unmeth","Paciente13.Meth", "Paciente13.Unmeth","Paciente14.Meth", "Paciente14.Unmeth","Paciente15.Meth", "Paciente15.Unmeth","Paciente16.Meth", "Paciente16.Unmeth","Paciente17.Meth", "Paciente17.Unmeth","Paciente18.Meth", "Paciente18.Unmeth","Paciente19.Meth", "Paciente19.Unmeth","Paciente20.Meth", "Paciente20.Unmeth","Paciente21.Meth", "Paciente21.Unmeth","Paciente22.Meth", "Paciente22.Unmeth","Paciente23.Meth", "Paciente23.Unmeth")

write.table(x = gse78751idc, file = "gse78751idc.txt", row.names = FALSE, sep = "\t") # Guardamos
