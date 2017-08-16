#########################################################################################
##                                                                                     ##
## s00_gse78751.R                                                                      ##
## Descarga de datos para el estudio GSE78751                                          ##
## Autor: Antonio Manuel Trassierra Fresco                                             ##
##                                                                                     ##
#########################################################################################

library(GEOquery) # Cargamos la librería necesaria

setwd("/home/atrassierra")

getGEOSuppFiles("GSE78751") # Los datos se descargan en el directorio de trabajo actual