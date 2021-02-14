args = commandArgs(trailingOnly = TRUE)
wd <- args[2]
software_deps <- args[3]

# WGCNA: en windows
# Es necesario actualizar R a la versión 4.0.3

# No ejecutar esto si se tiene la última versión de R

################## Actualir R en windows ################
# Instala / carga el paquete

if(!require(installr)) {

  install.packages("installr", repos = 'http://cran.us.r-project.org')
  require(installr)
  } 

# actualizar R
library("installr")

########################################################


################### Rtool para BiocManager #######################
if (!requireNamespace ("BiocManager", quietly = TRUE)) 

install.packages("BiocManager", repos = 'http://cran.us.r-project.org')

BiocManager::install(c("RTCGAToolbox", "WGCNA", "DESeq2", "coexnet"),  lib = software_deps) 

install.packages ("dplyr", repos = 'http://cran.us.r-project.org',  lib = software_deps)
install.packages("factoextra", repos = 'http://cran.us.r-project.org',  lib = software_deps)
install.packages("DCGL", repos = 'http://cran.us.r-project.org',  lib = software_deps)
install.packages("ggplot2", repos = 'http://cran.us.r-project.org',  lib = software_deps)

###########################################################


