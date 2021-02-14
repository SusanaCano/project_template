# !/bin/bash 

##############################################
# Cargamos librerias

chmod 755 results

WORKINGDIR=$1 # Directorio de la carpeta "code"
CODEDIR=$2
DEPSDIR=$3 # Directorio de la carpeta "deps"
RESULTSDIR=$4 # Directorio de la carpeta "results"


sh "setup.sh" $WORKINGDIR $CODEDIR $DEPSDIR

#############################################


Rscript "2_AnalisisCoexWGCNA.R" $WORKINGDIR $CODEDIR $DEPSDIR $RESULTSDIR

#Rscript "5_Cytoscape.R" $WORKINGDIR $CODEDIR $DEPENDENCIASDIR $RESULTSDIR

