# WGCNA: en windows
# Es necesario actualizar R a la versión 4.0.3

# No ejecutar esto si se tiene la última versión de R

################## Actualir R en windows ################
# Instala / carga el paquete
if(!require(installr)) {
  install.packages("installr")
  require(installr)
  } 


###
# actualizar R
#install.packages("installr")

library("installr")

########################################################


################### Rtool para BiocManager #######################
if (!requireNamespace ("BiocManager", quietly = TRUE)) 

install.packages ("BiocManager")

BiocManager :: install ("RTCGAToolbox")
BiocManager :: install("WGCNA")
BiocManager :: install ("DESeq2")

install.packages ("dplyr")
install.packages("factoextra")
install.packages("NbClust")
install.packages("DCGL")
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




