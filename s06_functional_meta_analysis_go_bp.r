#########################################################################################
##                                                                                     ##
## s06_functional_meta_analysis_go_bp.r                                                ##
## Metaanálisis funcional para estudios de metilación                                  ##
## Autores: Antonio Manuel Trassierra Fresco, Francisco García García                  ##
## Anotación GO (bp)                                                                   ##
##                                                                                     ##
#########################################################################################

date ()
Sys.info ()[c("nodename", "user")]
commandArgs ()
R.version.string
setwd("/home/antonio/Escritorio/atrassierra/meta_analysis/scripts")
try (source (".job.r")); try (.job)


# Limpiamos el espacio de trabajo
rm (list = ls ())

# Cargamos las librerías necesarias
# Estos son los paquetes necesarios para descargar y cargar
library(Biobase); packageDescription("Biobase", fields = "Version")
library(metafor); packageDescription("metafor", fields = "Version")
library(mdgsa); packageDescription("mdgsa", fields = "Version")
library(RamiGO); packageDescription("RamiGO", fields = "Version")
library(ggplot2); packageDescription("ggplot2", fields = "Version")
library(reshape); packageDescription("reshape", fields = "Version")




# Paso 1. Preparamos el input para el metaanálisis: matrices LOR y SD
######################################################################

### Necesitamos tener previamente para cada comparación/estudio los resultados de GSA en un fichero .RData con esta estructura: 
####            N        lor        pval      padj        sd          t conv
#### GO:0000002 12  0.4373106 0.128203007 1.0000000 0.2874434  1.5213797    1
#### GO:0000018 29  0.2177107 0.238623945 1.0000000 0.1847325  1.1785188    1


### Cargamos los resultados de los GSA para cada estudio
setwd (file.path (.job$dir, "results", "gsa_cutoff0_normypromedio"))
ficheros <- dir(pattern= "go_bp_")
ficheros

# Buscamos una lista que incluya todos los GOs únicos para todos los estudios
functions<-NULL
for (fi in ficheros){
  load (fi)
  functions <- c(functions, rownames (res))
}
length (functions)
functions <- unique (functions)
length (functions)
functions <- sort (functions)


### Generamos una matriz con todos los LOR para todos los estudios
mat.lor <- matrix (NA, nrow = length (functions), ncol = length (ficheros))
rownames (mat.lor) <- functions
colnames (mat.lor) <- strsplit (ficheros, ".RData")
head (mat.lor)

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, ".RData"))
  lor <- res$lor
  names (lor) <- (rownames (res))
  mat.lor[, co] <- lor[rownames(mat.lor)] 
}

head (mat.lor)
table (is.na(mat.lor))
dim (mat.lor)
summary(mat.lor)


### Generamos una matriz con todas las SD para todos los estudios
mat.sd <- matrix (NA, nrow = length (functions), ncol = length (ficheros))
rownames (mat.sd) <- functions
colnames (mat.sd) <- strsplit (ficheros, ".RData")
head (mat.sd)

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, ".RData"))
  sd <- res$sd
  names (sd) <- (rownames (res))
  head (sd)
  mat.sd[, co] <- sd[rownames(mat.sd)] 
}

head (mat.sd)
table (is.na(mat.sd))
dim (mat.sd)
summary(mat.sd)


### Generamos una matriz con todos los p valores ajustados para todos los estudios
mat.adjp <- matrix (NA, nrow = length (functions), ncol = length (ficheros))
rownames (mat.adjp) <- functions
colnames (mat.adjp) <- strsplit (ficheros, ".RData")
head (mat.adjp)

for (fi in ficheros){
  load (fi)
  co <- as.character(strsplit (fi, ".RData"))
  adjp <- res$padj
  names (adjp) <- (rownames (res))
  head (adjp)
  mat.adjp[, co] <- adjp[rownames(mat.adjp)] 
}

head (mat.adjp)
table (is.na(mat.adjp))
dim (mat.adjp)
table(mat.adjp < 0.05)

mat.adj2 <- (mat.adjp < 0.05)
head(mat.adj2)
total <- apply(mat.adj2, 1, sum)
head(total)
summary(total)
table(total)





# Paso 2. Metaanálisis para los términos funcionales
#################################################################

# En general, suponemos que la varianza entre estudios no es cero.
# hay diferentes métodos para estimar esta varianza:
# DL (Dersimonian-Laird), REML (Restricted maximum-likelihood, default)....
# result.lor <- rma(yi = mat.lor[1, ], vi = mat.sd[1, ], method = "DL") 
# DerSimonian-Laird. Y -> log(OR) V -> var(log(OR)

# Evaluaremos también el modelo de efectos fijos (FE)


# Parámetros para seleccionar

methods = c("DL", "HE", "HS", "FE")  # Métodos para el metaanálisis

corte = 0.05  # Punto de corte para detectar resultados significativos
or    = 0.5   # Punto de corte para detectar OR 
adj.p.value = "fdr"  # Método para ajustar los p-valores

# Preparamos una matriz para los resultados
res <- matrix (NA, nrow = length (methods), ncol = 6)
colnames(res) <- c("over", "under", "sig.over", "sig.under", "sig.or.over", "sig.or.under")
rownames(res) <- methods
res

# Nota importante: tenemos que indicar que nuestra medida de variabilidad esla desviación estándar y no la varianza
# pondremos "sei" (por defecto es "vi"):

# Metaanálisis:
for (i in methods){
  meta_analisis <- lapply(1:length(rownames(mat.lor)),
                          function(x){yi = rma(mat.lor[x, ], sei =mat.sd[x, ],
                                               method = i)})
  names (meta_analisis) <- rownames(mat.lor)
  result_meta <- as.data.frame(do.call("rbind", lapply(meta_analisis,function(x)
  {c(x$ci.lb, x$b, x$ci.ub, x$pval, x$QE, x$QEp, x$se, x$tau2, x$I2, x$H2) })))
  colnames(result_meta) <- c("lower_bound", "summary_LOR", "upper_bound", "pvalue", 
                             "QE", "QEp", "SE", "tau2", "I2", "H2")
  p.adjust <- p.adjust(result_meta[, 4], method= adj.p.value)  
  result_meta <- round(cbind (result_meta, p.adjust),3)
  
  
  if (file.exists(file.path (.job$dir, "meta", "go", "bp", "files"))){
    setwd(file.path(.job$dir, "meta", "go", "bp", "files"))
  } else {
    dir.create(file.path(.job$dir, "meta", "go", "bp", "files"), recursive = T)
    setwd(file.path(.job$dir, "meta", "go", "bp", "files"))    
  }
  
  res[i, "over"]      <-  sum(result_meta[, "summary_LOR"] > 0)
  res[i, "under"]    <-  sum(result_meta[, "summary_LOR"] < 0)
  res[i,"sig.over"]   <-  sum(result_meta[, "summary_LOR"] > 0 & result_meta[,"p.adjust"] < corte)
  res[i,"sig.under"] <-  sum(result_meta[, "summary_LOR"] < 0 & result_meta[,"p.adjust"] < corte)
  res[i,"sig.or.over"]   <-  sum(result_meta[, "summary_LOR"] >  or & result_meta[,"p.adjust"] < corte)
  res[i,"sig.or.under"] <-  sum(result_meta[, "summary_LOR"] < -or & result_meta[,"p.adjust"] < corte)
  
  # Obteniendo los IDs comunes para todos los modelos 
  if (i == "DL"){
    DL <- result_meta    
    DL[,"name"] <- getGOnames(rownames(DL))
    DL[, "ID"]   <- rownames(DL)
    DL <- DL[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust", "QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(DL, "all.results.DL.txt", sep = "\t", quote = F, row.names = F)
    DLs = subset(DL, DL[,"p.adjust"] < corte)
    write.table(DLs, "sig.results.DL.txt", sep = "\t", quote = F, row.names = F)
    idDL = rownames(DLs)
  }
  if (i == "HE"){
    HE <- result_meta    
    HE[,"name"] <- getGOnames(rownames(HE))
    HE[, "ID"]   <- rownames(HE)
    HE <- HE[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(HE, "all.results.HE.txt", sep = "\t", quote = F, row.names = F)
    HEs = subset(HE, HE[,"p.adjust"] < corte)
    write.table(HEs, "sig.results.HE.txt", sep = "\t", quote = F, row.names = F)
    idHE = rownames(HEs)
  }
  if (i == "HS"){
    HS <- result_meta
    HS[,"name"] <- getGOnames(rownames(HS))
    HS[, "ID"]   <- rownames(HS)
    HS <- HS[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(HS, "all.results.HS.txt", sep = "\t", quote = F, row.names = F)
    HSs = subset(HS, HS[,"p.adjust"] < corte)
    write.table(HSs, "sig.results.HS.txt", sep = "\t", quote = F, row.names = F)
    idHS = rownames(HSs)
  }
  if (i == "FE"){
    FE <- result_meta
    FE[,"name"] <- getGOnames(rownames(FE))
    FE[, "ID"]   <- rownames(FE)
    FE <- FE[c("ID", "name", "lower_bound", "summary_LOR", "upper_bound", "pvalue",
               "p.adjust","QE", "QEp", "SE", "tau2", "I2", "H2" )]
    write.table(FE, "all.results.FE.txt", sep = "\t", quote = F, row.names = F)
    FEs = subset(FE, FE[,"p.adjust"] < corte)
    write.table(FEs, "sig.results.FE.txt", sep = "\t", quote = F, row.names = F)
    idFE = rownames(FEs)
  }
}



## Resultados globales
setwd (file.path (.job$dir, "meta", "go", "bp", "files"))


# A. Todos los resultados para los diferentes modelos
print(res)
cat ("method\t", file = "all.results.all.methods.txt")
write.table (res, file = "all.results.all.methods.txt",
             append = TRUE, quote = FALSE, sep = "\t", 
             row.names = TRUE, col.names = TRUE)



# B. Intersección de los resultados significativos
intersectmodels =  Reduce(intersect,  list(idDL, idHE, idHS, idFE))
length(intersectmodels)
write.table (intersectmodels, file = "sig_res_intersection.txt",
             quote = FALSE, sep = "\t", row.names = F, col.names = F)



# C. Unión de los resultados significativos
unionmodels = Reduce(c, list(idDL, idHE, idHS, idFE))

mat.union <- as.data.frame(table(unionmodels))
dim(mat.union)
mat.union <- mat.union[order(mat.union$Freq, decreasing = T),]
cat ("ID\tFreq\n", file = "sig_res_union.txt")
write.table (mat.union, file = "res_union.txt",append = TRUE,
             quote = FALSE, sep = "\t", row.names = F, col.names = F)






# Paso 3. Representación gráfica
#################################################################

if (file.exists(file.path (.job$dir, "meta", "go", "bp", "plots"))){
  setwd(file.path(.job$dir, "meta", "go", "bp", "plots"))
} else {
  dir.create(file.path(.job$dir, "meta", "go", "bp", "plots"), recursive = T)
  setwd(file.path(.job$dir, "meta", "go", "bp", "plots"))    
}




## 3.1. Plots para explorar los datos de entrada: matriz LOR y SD 
##################################################################
# Estos plots son los mismos para cualquier método. 
# Esto cambia cuando tenemos diferentes bases de datos biológicas 
# porque tendremos diferentes tipos de resultados del modelo logístico

gmat.sd  <- mat.sd
colnames(gmat.sd) <-  toupper(colnames(mat.sd))
g2mat.sd <- melt(gmat.sd)
colnames(g2mat.sd) <- c("function", "study", "sd")
head(g2mat.sd); dim(g2mat.sd)

gmat.lor  <- mat.lor
colnames(gmat.lor) <-  toupper(colnames(mat.lor))
g2mat.lor <- melt(gmat.lor)
colnames(g2mat.lor) <- c("function", "study", "lor")
head(g2mat.lor); dim(g2mat.lor)

x.por <- 3
y.por <- 2

jpeg (filename = "go_sd.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(g2mat.sd, aes(study, g2mat.sd[, "sd"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("Estudios") 
p <- p + ylab("Error estándar")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) # Para cambiar el tamaño del título en el eje X
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) # Para cambiar el tamaño del título en el eje Y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 20 )) # Para cambiar el tamaño de los números
p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) # Para cambiar el tamaño de los números
print(p)
dev.off ()

jpeg (filename = "go_lor.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(g2mat.lor, aes(study, g2mat.lor[, "lor"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("Estudios") 
p <- p + ylab("Log Odds Ratio")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) # Para cambiar el tamaño del título en el eje X
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) # Para cambiar el tamaño del título en el eje Y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 20 )) # Para cambiar el tamaño de los números
p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) # Para cambiar el tamaño de los números
print(p)
dev.off ()



## 3.2. Plots o gráficos para explorar los datos de salida: LOR vs FDR (volcano plot) 
######################################################################################

# Este gráfico muestra los resultados finales para cada método
# Queremos conocer la relación entre el LOR combinado y su SE

# Seleccionamos los resultados para un método específico:
res.method = DL #  results from Dersimonian-Laird method

summary(res.method[, "summary_LOR"])
summary(res.method[, "p.adjust"])
log10 <- -log(res.method[, "p.adjust"], base = 10)
summary(log10)
corte.logfdr <- 3
table(log10 > 3)
head(res.method)

# Asignamos colores. El gráfico muestra otros por defecto de manera que lo ajustamos

res.method$threshold = rep("white", nrow(res.method))
sig.up   <- (res.method[, "p.adjust"] < 0.05) & (res.method[, "summary_LOR"] > 0)
sig.down <- (res.method[, "p.adjust"] < 0.05) & (res.method[, "summary_LOR"] < 0)

res.method$threshold [sig.up]    <- "blue"
res.method$threshold [sig.down]  <- "red"
table(res.method$threshold)
res.method2 <- res.method[log10 < corte.logfdr,]
table(res.method2$threshold)

dim(res.method2)
res.method$threshold <- as.factor(res.method$threshold)

## Construimos el objeto plot
x.por <- 3
y.por <- 3
jpeg (filename = "go_volcano.jpeg",  width = 480 *x.por ,   height = 480*y.por, pointsize = 100, quality = 100) 
g <- ggplot(data=res.method2, aes(x=summary_LOR, y=-log10(p.adjust), colour= threshold)) #generate data
g <- g + geom_point(alpha=0.8, size=4) 
g <- g + xlim(c(-2, 2)) + ylim(c(0, 3.1))
g <- g + xlab("log2 Odds Ratio") 
g <- g + theme(axis.title.x= element_text(colour= "black", size = 35 )) # Para cambiar el tamaño del título en el eje X
g <- g + ylab("-log10 FDR") 
g <- g + theme(axis.title.y= element_text(colour= "black", size = 35 )) # Para cambiar el tamaño del título en el eje Y 
g <- g + labs(title = "Volcano plot")
g <- g + theme(title= element_text(colour= "black", size = 35 )) # Para cambiar el tamaño del título 
g <- g + theme(axis.text.x = element_text(colour= "black", size = 35 )) # Para cambiar el tamaño de los números
g <- g + theme(axis.text.y = element_text(colour= "black", size = 35 )) # Para cambiar el tamaño de los números
g <- g + geom_vline(xintercept= c(0), colour = "black", size = 0.5)
g <- g + geom_hline(yintercept= -log(0.05, base= 10), colour = "black", size = 0.5)
g <- g + theme(legend.position = "none") 
print(g)
dev.off ()


## 3.3. Gráfico/plot para evaluar la heterogeneidad para todos los modelos 
###########################################################################

methods<-  c(rep("DL", nrow(DL)),rep("HE", nrow(HE)),rep("HS", nrow(HS)), 
             rep("FE", nrow(FE)) )
# QEp   <- c(DL$QEp, HE$QEp, HS$QEp, SJ$QEp, ML$QEp, REML$QEp, EB$QEp, PM$QEp, FE$QEp)
QEp      <- c(DL$QEp, HE$QEp, HS$QEp,                   FE$QEp)
SE       <- c(DL$SE, HE$SE, HS$SE,                      FE$SE)
I2       <- c(DL$I2, HE$I2, HS$I2,                      FE$I2)
H2       <- c(DL$H2, HE$H2, HS$H2,                      FE$H2)
tau2     <- c(DL$tau2, HE$tau2, HS$tau2,                FE$tau2)

datos <- data.frame (cbind(as.factor(methods), QEp, SE, I2, H2, tau2))
head(datos)
dim(datos)
summary(datos)
param <- c("QEp", "SE", "I2", "H2", "tau2")
x.por <- 3
y.por <- 2


# Dibujamos de manera separada QEp porque necesitamos un abline de 0.05
jpeg (filename = "go_het_QEp.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(datos, aes(methods, datos[, "QEp"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("MÃ©todos") 
p <- p + ylab("QEp")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) # Cambiamos el tamaño del título en el eje X
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) # Cambiamos el tamaño del título en el eje Y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 30 )) # Cambiamos el tamaño de los números
p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) # Cambiamos el tamaño de los números
p <- p + geom_hline(yintercept = 0.05, col = "red")
print(p)
dev.off ()

# Y ahora el resto:

for (i in 2:length(param)){
  jpeg (filename = paste("go_het_", param[i], ".jpeg", sep = ""),  
        width = 480 *x.por ,   height = 480*y.por,
        pointsize = 150, quality = 100) 
  p <- ggplot(datos, aes(methods, datos[, param[i]]))
  p <- p + geom_boxplot(col = "blue")
  p <- p + xlab("Métodos") 
  p <- p + ylab(param[i])
  p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) # Cambiamos el tamaño del título en el eje X
  p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) # Cambiamos el tamaño del título en el eje Y
  p <- p + theme(axis.text.x = element_text(colour= "black", size = 30 )) # Cambiamos el tamaño de los números
  p <- p + theme(axis.text.y = element_text(colour= "black", size = 30 )) # Cambiamos el tamaño de los números
  print(p)
  dev.off ()
}







## 3.4. Gráfico para evaluar la heterogeneidad y la influencia de cada modelo específico y función
################################################################################################### 

# Seleccionamos los resultados para un método específico:
res.method = DLs # significant results from Dersimonian-Laird method
metodo <- "DL"   # method to estimate the variability

sig.fun <- rownames(res.method)
x.por = 2.5; y.por = 2.5


if (length(sig.fun) == 0){print ("Not significant results")} else {
  for (i in 1:length(sig.fun)){
    # Ajustamos el modelo. Lo usaremos para varios gráficos
    res <- rma(yi= mat.lor[sig.fun[i],], sei =mat.sd[sig.fun[i],], method = metodo)
    
    ## A. Gráficos de bosque (Información detallada del tamaño del efecto para cada estudio)
    png (filename = paste("go_forest_", substr(sig.fun[i],4,10),".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    par(mar=c(4,4,1,2))
    forest(res, slab = toupper(colnames(mat.lor)),
           xlab="Odds Ratio", cex=0.7,
           mlab=paste(metodo, "Model for All Studies", sep = " "), col = "red", 
           main = paste("\n", sig.fun[i], " (", res.method[i,"name"], ")",sep=""))    
    text( 9,-3, "Odds Ratio [IC 95%]", pos=2, cex = 0.7)
    dev.off()
    


    ## B. Gráficos de embudo (información detallada para la medida del efecto en cada estudio)
    # Un gráfico de embudo muestra las medidas del efecto observadas o resultantes en el eje x frente 
    # a alguna medida de precisión de la medida del efecto observado o resultante en el eje y. 
    # Basado en Sterne and Egger (2001)
    # http://www.metafor-project.org/doku.php/plots:funnel_plot_variations
    png (filename = paste("go_funnel_",substr(sig.fun[i],4,10),".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    ### Preparamos un array 2x2 para hacer el gráfico
    par(mfrow=c(2,2))
    ### Dibujamos los gráficos de embudo
    funnel(res, main="Standard Error", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    funnel(res, yaxis="vi", main="Sampling Variance", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    funnel(res, yaxis="seinv", main="Inverse Standard Error", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    funnel(res, yaxis="vinv", main="Inverse Sampling Variance", back ="darkslategray1",
           xlab = paste("LOR (", sig.fun[i], ")",sep =""))
    dev.off()
   
    ## C. Gráfico radial 
    # Los gráficos radiales fueron introducidos por Rex Galbraith (1988a, 1988b, 1994) 
    # y pueden ser útiles para en el contexto del metaanálisis para examinar los datos en busca de heterogeneidad.
    # http://www.metafor-project.org/doku.php/plots:radial_plot
    
    ### Para guardarlo como png
    png (filename = paste("GO_radial_",substr(sig.fun[i],4,10),".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    ### Ajustamos los márgenes para aprovechar mejor el espacio
    par(mar=c(5,4,0,2))
    ### Dibujamos el gráfico radial
    radial(res,back ="darkslategray1",
           main = paste("\n", sig.fun[i], " (", res.method[i,"name"], ")",sep=""))
    dev.off()
    
    
    ## D. Gráfico de influencia 
    # Muestra varias medidas de diagnóstico
    ### Para guardarlo como png
    png (filename = paste("go_influ_",substr(sig.fun[i],4,10),".png",sep =""),  
         width = 480 * x.por, height = 480 * y.por, res = 200)
    inf <- influence(res)
    plot(inf, plotfb = T)
    dev.off()
  }
}



## 3.5. Análisis de sensibilidad para un modelo específico 
########################################################### 

# Primero calculamos la desviación estándar estimada para cada método:

methods = c("DL", "HE", "HS",  "FE")  

mat.sen <- as.data.frame(matrix (NA, nrow = nrow(mat.lor), ncol = length (methods)))
colnames(mat.sen) <- methods
dim(mat.sen)
head(mat.sen)


for (i in methods){
  meta_analisis <- lapply(1:length(rownames(mat.lor)),
                          function(x){yi = rma(mat.lor[x, ], sei =mat.sd[x, ],
                                               method = i)})
  for (j in 1:length(meta_analisis)){
    mat.sen[j,i]  <- sd(leave1out(meta_analisis[[j]])$se)
    }
  print(i)
}

summary(mat.sen)
dim(mat.sen)
boxplot(mat.sen)
# Sólo modelos de efectos aleatorios:
mat.sen <- mat.sen[, c("DL", "HE", "HS")]

gmat.sen  <- mat.sen
g2mat.sen <- melt(gmat.sen)
head(g2mat.sen)
colnames(g2mat.sen) <- c("metodo", "ee_lor")
head(g2mat.sen); dim(g2mat.sen)

setwd(file.path(.job$dir, "meta", "go", "bp", "plots"))
x.por <- 3
y.por <- 2

jpeg (filename = "go_sens_allmethods.jpeg",  
      width = 480 *x.por ,   height = 480*y.por,
      pointsize = 150, quality = 100) 
p <- ggplot(g2mat.sen, aes(metodo, g2mat.sen[, "ee_lor"]))
p <- p + geom_boxplot(col = "blue")
p <- p + xlab("Método") 
p <- p + ylab("Error estándar del LOR")
p <- p + theme(axis.title.x= element_text(colour= "black", size = 40 )) # Para cambiar el tamaño del título del eje X
p <- p + theme(axis.title.y= element_text(colour= "black", size = 40 )) # Para cambiar el tamaño del título del eje Y
p <- p + theme(axis.text.x = element_text(colour= "black", size = 20 )) # Para cambiar el tamaño de los números
p <- p + theme(axis.text.y = element_text(colour= "black", size = 20 )) # Para cambiar el tamaño de los números
print(p)
dev.off ()

# Seleccionando un método:
res.method = DLs # Resultados significativos del método de Dersimonian-Laird
metodo <- "DL"   # Método para estimar la variabilidad
sig.fun <- rownames(res.method)



setwd(file.path(.job$dir, "meta", "go", "bp", "files"))   

if (length(sig.fun) == 0){print ("Not significant results")} else {
  for (i in 1:length(sig.fun)){
    # Ajustamos el modelo.Lo usarmoes para varios gráficos.
    res <- rma(yi= mat.lor[sig.fun[i],], sei =mat.sd[sig.fun[i],], method = metodo)
    res.l1out  <- leave1out(meta_analisis[[i]])
    resu.l1out <- round(print.list.rma(res.l1out),3)
    studies <- toupper(substr(colnames(mat.lor),1,4))
    resu.l1out <- cbind(studies,resu.l1out)
    write.table(resu.l1out, 
                paste("go_sensi_", substr(sig.fun[i],4,10),".txt", sep =""),
                quote = F, row.names = F)
  }
}
    
  

## 3.6. Visualización de los términos Gos significativos juntos 
################################################################


setwd(file.path(.job$dir, "meta", "go", "bp", "plots"))   

ind1 <- res.method[,"p.adjust"] < 0.05 & res.method[,"summary_LOR"] > 0  
GO_over <- rownames(res.method)[ind1]
length(GO_over)
color <- rep("green", length(GO_over))
pngRes <- getAmigoTree(goIDs=GO_over, color=color, filename= "sig_go_over.png", 
                       picType="png", saveResult=TRUE)


ind2 <- res.method[,"p.adjust"] < 0.05 & res.method[,"summary_LOR"] < 0  
table(ind2)
GO_under <- rownames(res.method)[ind2]
length(GO_under)
color <- rep("red", length(GO_under))
pngRes <- getAmigoTree(goIDs=GO_under, color=color, filename= "sig_go_under.png", 
                       picType="png", saveResult=TRUE)





# Salida
warnings ()
sessionInfo ()
q ("no")
