#########################################################################################
##                                                                                     ##
## s01_gse52865.R                                                                      ##
## Este script descomprime el archivo donde se encuentra el estudio y lo procesa,      ##
## dejando las señales de metilación y no metilación para su posterior normalización.  ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

library(GEOquery)

setwd("/home/atrassierra/GSE52865") # Cambiamos el directorio de trabajo al del experimento en cuestión

gunzip('GSE52865_complete_data.txt.gz') # Descomprimimos el archivo para tener acceso al experimento.

# Leemos el archivo del experimento y extraemos las señales de metilación y no metilación

gse52865 = read.csv("GSE52865_complete_data.txt", sep = "\t")  
metilado <- grep("\\.meth$", colnames(gse52865))
unmetil <- grep("unmeth", colnames(gse52865))
combinado <- as.vector(rbind(metilado, unmetil)) 
gse52865 <- gse52865[, c(1, combinado)] # Matriz con las intensidades de señal, el 1 es para seguir teniendo los id de las sondas.


# Separamos el estudio para comparar los controles con los distintos tipos de cáncer de mama

basal <- gse52865[,c(rep(36:41), 48, 49, 64, 65, rep(68:71), rep(80:83), 86, 87, 94, 95, rep(102:105), 112, 113)]
her2 <- gse52865[, c(42, 43, 62, 63, 66, 67)]
luma <- gse52865[, c(44, 45, rep(50:53), 58, 59, 72, 73, 84, 85, 88, 89, rep(96:101), rep(108:111))]
lumb <- gse52865[, c(60, 61, rep(74:79), rep(90:93), 114, 115)]
controles <- gse52865[, c(rep(2:35))]

######################
## control vs basal ##
######################
gse52865basal <- cbind(gse52865[, 1], controles, basal)
colnames(gse52865basal)[1] <- "TargetID"

# Guardamos el data.frame resultante en un fichero de texto
write.table(x = gse52865basal, file = "gse52865basal.txt", row.names = FALSE, sep = "\t")

#####################
## control vs her2 ##
#####################
gse52865her2 <- cbind(gse52865[, 1], controles, her2)
colnames(gse52865her2)[1] <-"TargetID"

# Guardamos el data.frame resultante en un fichero de texto
write.table(x = gse52865her2, file = "gse52865her2.txt", row.names = FALSE, sep = "\t")

#####################
## control vs luma ##
#####################
gse52865luma <- cbind(gse52865[, 1], controles, luma)
colnames(gse52865luma)[1] <- "TargetID"

# Guardamos el data.frame resultante en un fichero de texto
write.table(x = gse52865luma, file = "gse52865luma.txt", row.names = FALSE, sep = "\t")

#####################
## control vs lumb ##
#####################
gse52865lumb <- cbind(gse52865[, 1], controles, lumb)
colnames(gse52865lumb)[1] <- "TargetID"

# Guardamos el data.frame resultante en un fichero de texto
write.table(x = gse52865lumb, file = "gse52865lumb.txt", row.names = FALSE, sep = "\t")