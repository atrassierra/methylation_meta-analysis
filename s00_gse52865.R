#########################################################################################
##                                                                                     ##
## s00_gse52865.R                                                                      ##
## Descarga de datos para el estudio GSE52865                                          ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

library(GEOquery) # Cargamos la librería necesaria

setwd("/home/atrassierra")

getGEOSuppFiles("GSE52865") # Los datos se descargan en el directorio de trabajo actual