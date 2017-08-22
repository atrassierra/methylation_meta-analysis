#########################################################################################
##                                                                                     ##
## s04_gse52865her2.R                                                                  ##
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

setwd('/home/atrassierra/GSE52865')

load('dmrs52865her2.RData')

# Preparamos la anotación y matcheamos los genes
anno <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene) # Preparamos el input para matchGenes
genes <- matchGenes(dmrs$table, subject = anno) # Mactheamos genes con sus dmrs

data <- dmrs$table
data <- cbind(data, genes$name) # Juntamos la información de los genes en una misma tabla
names(data)[length(names(data))] <- "genes" # Cambiamos el nombre a la columna de los genes

datasorted <- data[order(-data$value),] # Ordenamos de mayor a menor nivel de metilación

# Escribimos el fichero

write.table(x = datasorted, file = "dmrgen52865her2.txt", row.names = FALSE, sep = "\t")

# Duplicados. Varias DMRs pueden apuntar a un mismo gen. Si esto ocurre, aunque no es lo más correcto,
# promediaremos la señal de metilación porque hay veces que la señal apunta a más metilaciónn en casos o 
# en controles y no sabemos como interpretarlo.

datasorted <- datasorted[, c("genes", "value")]

# Promedio

promedio <- duplicated(data$genes)
table(promedio)
genesdup <- data$genes[promedio]
length(genesdup)

ind1 <- data$genes %in% genesdup # Creamos una variable para saber cuantos duplicados hay en el data.frame
table(ind1)

dataprom <- data[ind1,] # Creamos un nuevo data.frame con los genes duplicados
head(dataprom)
dataprom <- dataprom[order(dataprom$genes),] # Con esto podemos ver los genes duplicados
head(dataprom)
dataprom <- dataprom[, c("genes", "value")]

b <- aggregate(x = dataprom, by = list(genes = dataprom$genes), FUN = mean) # Promediamos las señales que apuntan a un mismo gen
b <- b[, -2]

genestotales <- data.frame(data$genes, data$value)
colnames(genestotales) <- c("genes", "value")

'%!in%' <- function(x, y)!('%in%'(x, y)) # Función not in creada

ind2 <- genestotales$genes %!in% genesdup # Obtenemos los genes que no están duplicados
genestotales <- genestotales[ind2,]
table(genestotales$genes %in% dataprom$genes)

final <- rbind(genestotales, b) # Juntamos todos los genes
final <- final[order(-final$value),] # Ordenamos de mayor a menor nivel de metilación diferencial
write.table(x = final, file = "genes52865her2prom.txt", row.names = FALSE, sep = "\t", quote = FALSE) # Input con el promedio para GSA

# Valor absoluto

final <- datasorted %>% group_by(genes) %>% summarise_each(funs(.[which.max(abs(.))])) # Nos quedamos con la DMR con mayor valor absoluto para cada gen
# Guardamos el resultado
write.table(x = final, file = "genes52865her2.txt", sep = "\t", row.names = FALSE, quote = FALSE) # Guardamos input para gsea con valor absoluto
save(final, file = "genes52865her2.RData")
