# WGCNA: en windows
# Es necesario actualizar R a la versión 4.0.3

args<-commandArgs(trailingOnly=TRUE)
wd<-args[2]
libdir<-args[3]


################## Actualir R en windows ################
# Instala / carga el paquete
if(!require(installr)) {
  install.packages("installr", repos='http://cran.us.r-project.org', lib = libdir)
  require(installr)
  } 


###

library("installr")

########################################################


################### Rtool para BiocManager #######################
if (!requireNamespace ("BiocManager", quietly = TRUE)) 

install.packages ("BiocManager" , repos='http://cran.us.r-project.org')

BiocManager::install("RTCGAToolbox", lib = libdir)
BiocManager::install("WGCNA", lib = libdir)
BiocManager::install("DESeq2", lib = libdir)
BiocManager::install("coexnet", lib = libdir)


install.packages ("dplyr", repos='http://cran.us.r-project.org', lib = libdir)
install.packages("factoextra", repos='http://cran.us.r-project.org', lib= libdir)
install.packages("NbClust", repos='http://cran.us.r-project.org', lib = libdir)
install.packages("DCGL", repos='http://cran.us.r-project.org', lib = libdir)

###########################################################

### cargamos las librerias:

library("WGCNA") # Análisis de redes de coexpresión de genes ponderados
library("dplyr") # Manipulación de dataframes
library("factoextra") # Extracción y visualice los resultados de análisis de datos multivariados
library("NbClust") # Determinar el mejor número de clústeres en un conjunto de datos
library("cluster") # Encontrar grupos en los datos
library("DESeq2") # Análisis de datos de RNA-seq
library("DCGL") # Análisis de coexpresión diferencial y análisis de regulación diferencial de 
                # datos de microarrays de expresión génica
library("coexnet")


