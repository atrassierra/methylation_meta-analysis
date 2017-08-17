#########################################################################################
##                                                                                     ##
## s04_gse52865basal.R                                                                 ##
## Script para la anotaci�n de las DMRs, pondremos los genes a los que apuntan por     ##
## el m�todo del promedio y del valor absoluto. Adem�s, prepararemos el input para el  ##
## an�lisis de enriquecimiento funcional.                                              ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

# Cargamos las librer�as necesarias

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(bumphunter)
library(dplyr)

setwd('/home/atrassierra/GSE52865')

load('dmrs52865basal.RData')

# Preparamos la anotaci�n y matcheamos los genes
anno <- annotateTranscripts(TxDb.Hsapiens.UCSC.hg19.knownGene) # Preparamos el imput para matchGenes
genes <- matchGenes(dmrs$table, subject = anno) # Mactheamos genes con sus dmrs

data <- dmrs$table
data <- cbind(data, genes$name) # Juntamos la informaci�nn de los genes en una misma tabla
names(data)[length(names(data))] <- "genes" # Cambiamos el nombre a la columna de los genes

datasorted <- data[order(-data$value),] # Ordenamos de mayor a menor nivel de metilaci�n

# Escribimos el fichero

write.table(x = datasorted, file = "dmrgen52865basal.txt", row.names = FALSE, sep = "\t")

# Duplicados. Varias DMRs pueden apuntar a un mismo gen. Si esto ocurre, aunque no es lo m�s correcto,
# promediaremos la se�al de metilaci�n porque hay veces que la se�al apunta a m�s metilaci�nn en casos o 
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

b <- aggregate(x = dataprom, by = list(genes = dataprom$genes), FUN = mean) # Promediamos las se�ales que apuntan a un mismo gen
b <- b[, -2]

genestotales <- data.frame(data$genes, data$value)
colnames(genestotales) <- c("genes", "value")

'%!in%' <- function(x, y)!('%in%'(x, y)) # Funci�n not in creada

ind2 <- genestotales$genes %!in% genesdup # Obtenemos los genes que no est�n duplicados
genestotales <- genestotales[ind2,]
table(genestotales$genes %in% dataprom$genes)

final <- rbind(genestotales, b) # Juntamos todos los genes
final <- final[order(final$value),] # Ordenamos de mayor a menor nivel de metilaci�n diferencial
write.table(x = final, file = "genes52865basalprom.txt", row.names = FALSE, sep = "\t", quote = FALSE) # Input con el promedio para GSA

# Valor absoluto

final <- datasorted %>% group_by(genes) %>% summarise_each(funs(.[which.max(abs(.))])) # Nos quedamos con la DMR con mayor valor absoluto para cada gen
final <- final[order(-final$value),]
# Guardamos el resultado
write.table(x = final, file = "genes52865basal.txt", sep = "\t", row.names = FALSE, quote = FALSE) # Guardamos input para GSA con valor absoluto
save(final, file = "genes52865basal.RData")
