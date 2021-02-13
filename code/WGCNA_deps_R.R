# WGCNA: en windows
# Es necesario actualizar R a la versión 4.0.3

# No ejecutar esto si se tiene la última versión de R

################## Actualir R en windows ################
# Instala / carga el paquete

#if(!require(installr)) {
 # install.packages("installr")
  #require(installr)
  #} 

# setwd("D:/uni/2020-21/BiologiaSistemas/Practicas/Aurelio/trabajoFinal/Git/project_template_BS_SARS-CoV2/code")
###
# actualizar R
#install.packages("installr")

#library("installr")

########################################################


################### Rtool para BiocManager #######################
if (!requireNamespace ("BiocManager", quietly = TRUE)) 
  #if (!requireNamespace ("BiocManager", quietly = TRUE)) 
#install.packages ("BiocManager")
 install.packages("BiocManager", repos = 'http://cran.us.r-project.org')

BiocManager::install(c("RTCGAToolbox", "WGCNA", "DESeq2", "coexnet")) 

install.packages ("dplyr", repos = 'http://cran.us.r-project.org')
install.packages("factoextra", repos = 'http://cran.us.r-project.org')
install.packages("NbClust", repos = 'http://cran.us.r-project.org')
install.packages("DCGL", repos = 'http://cran.us.r-project.org')
install.packages("ggplot2", repos = 'http://cran.us.r-project.org')
###########################################################

### cargamos las librerias:

library("WGCNA") # Análisis de redes de coexpresión de genes ponderados
library("dplyr") # Manipulación de dataframes
library("ggplot2") # Gráficas
library("factoextra") # Extracción y visualice los resultados de análisis de datos multivariados
library("NbClust") # Determinar el mejor número de clústeres en un conjunto de datos
library("cluster") # Encontrar grupos en los datos
library("DESeq2") # Análisis de datos de RNA-seq
library("DCGL") # Análisis de coexpresión diferencial y análisis de regulación diferencial de 
                # datos de microarrays de expresión génica

library("coexnet") # Graficas de coexpresion de genes


