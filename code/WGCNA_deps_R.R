# WGCNA: en windows
# Es necesario actualizar R a la versión 4.0.3

# No ejecutar esto si se tiene la última versión de R

libdir<-"C:/Users/Ana/Desktop/Uni/PrimerCuatri/Biolog�aSist/TRABAJOFINAL/software/deps/"

################## Actualir R en windows ################
# Instala / carga el paquete
if(!require(installr)) {
  install.packages("installr", repos='http://cran.us.r-project.org', lib = libdir)
  require(installr)
  } 


###
# actualizar R
#install.packages("installr")

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


