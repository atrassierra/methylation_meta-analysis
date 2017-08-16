#########################################################################################
##                                                                                     ##
## s01_gse59901.R                                                                      ##
## Este script descomprime el archivo donde se encuentra el estudio y lo procesa,      ##
## dejando las señales de metilación y no metilación para su posterior normalización.  ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

library(GEOquery)

setwd("/home/atrassierra/GSE59901") # Cambiamos el directorio de trabajo al del experimento en cuestión

gunzip('GSE59901_450K_MethylationData_raw.csv.gz') # Descomprimimos el archivo para tener acceso al experimento.

# Leemos el archivo del experimento y extraemos las señales de metilación y no metilación

gse59901 = read.csv("GSE59901_450K_MethylationData_raw.csv")  

targetID <- gse59901[,2]

lista <-  grep("Signal", colnames(gse59901), invert = F)

gse59901signals <- gse59901[, lista] # Signal A no metilado, signal B metilado

# Separamos el estudio para comparar los controles con los distintos tipos de cáncer de mama

c <- c(c(9:12), c(31, 32, 37, 38))

brc <- c(c(1:6), c(19, 20, 25, 26, 27, 28, 33, 34, 39, 40))

ilb <- c(c(7,8, 13:18, 21, 22, 23, 24, 29, 30, 35, 36, 41:44))

idc <- c(45:64)

####################
## control vs brc ##
####################

gse59901brc <- gse59901signals[, c(c, brc)]

gse59901brc <- cbind(targetID, gse59901brc)

write.table(x = gse59901brc, file = 'tablebrc.txt', sep = "\t", row.names = FALSE) # Guardamos

####################
## control vs ilb ##
####################

gse59901ilb <- gse59901signals[, c(c, ilb)]

gse59901ilb <- cbind(targetID, gse59901ilb)

write.table(x = gse59901ilb, file = 'tableilb.txt', sep = "\t", row.names = FALSE) # Guardamos

####################
## control vs idc ##
####################

gse59901idc <- gse59901signals[, c(c, idc)]

gse59901idc <- cbind(targetID, gse59901idc)

write.table(x = gse59901idc, file = 'tableidc.txt', sep = "\t", row.names = FALSE) # Guardamos