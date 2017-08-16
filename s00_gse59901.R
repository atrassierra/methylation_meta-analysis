#########################################################################################
##                                                                                     ##
## s00_gse59901.R                                                                      ##
## Descarga de datos para el estudio GSE59901                                          ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

library(GEOquery) # Cargamos la librería necesaria

setwd("/home/atrassierra")

getGEOSuppFiles("GSE59901") # Los datos se descargan en el directorio de trabajo actual