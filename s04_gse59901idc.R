#########################################################################################
##                                                                                     ##
## s04_gse59901idc.R                                                                   ##
## Script para la anotación de las DMRs, pondremos los genes a los que apuntan por     ##
## el método del promedio y del valor absoluto. Además, prepararemos el input para el  ##
## análisis de enriquecimiento funcional.                                              ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

# Cargamos las librerías necesarias

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(bumphunter)
library(dplyr)

setwd('/home/atrassierra/GSE59901')

load('dmrs59901idc.RData')

# Preparamos la anotación y matcheamos los genes
anno <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene) # Preparamos el imput para matchGenes
genes <- matchGenes(dmrs$table, subject = anno) # Mactheamos genes con sus dmrs

data <- dmrs$table
data <- cbind(data, genes$name) # Juntamos la información de los genes en una misma tabla
names(data)[length(names(data))] <- "genes" # Cambiamos el nombre a la columna de los genes

datasorted <- data[order(-data$value),] # Ordenamos de mayor a menor nivel de metilación

# Escribimos el fichero

write.table(x = datasorted, file = "dmrgen59901idc.txt", row.names = FALSE, sep = "\t")

# Duplicados. Varias DMRs pueden apuntar a un mismo gen. Si esto ocurre, aunque no es lo más correcto,
# promediaremos la señal de metilación porque hay veces que la señal apunta a más metilaciónn en casos o 
# en controles y no sabemos como interpretarlo.

datasorted <- datasorted[, c("genes", "value")]

# Promedio

final <- datasorted %>% group_by(genes) %>% summarise(value = mean(value)) # Promediamos las DMR
final <- final[order(-final$value),] # Ordenamos de mayor a menor nivel de metilación diferencial
write.table(x = final, file = "genes59901idcprom.txt", row.names = FALSE, sep = "\t", quote = FALSE) # Input con el promedio para GSA
save(final, file = "genes59901idcprom.RData")

# Valor absoluto

final <- datasorted %>% group_by(genes) %>% summarise_each(funs(.[which.max(abs(.))])) # Nos quedamos con la DMR con mayor valor absoluto para cada gen
final <- final[order(-final$value),]
# Guardamos el resultado
write.table(x = final, file = "genes59901idc.txt", sep = "\t", row.names = FALSE, quote = FALSE) # Guardamos input para gsea con el valor absoluto
save(final, file = "genes59901idc.RData")

