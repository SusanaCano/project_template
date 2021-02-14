
# !/bin/bash 

# Carga e instalaciones de librerias en R

#SOFTWARE=$
#mkdir -p $SOFTWARE
#chmod 755 project_template_BS_SARS-CoV2

#mkdir $SOFTWARE software/deps
# chmod 755 software/deps
WORKINGDIR=$1

chmod 755 $WORKINGDIR
mkdir -p $WORKINGDIR/software

CODEDIR=$2

for software in $WORKINGDIR
do
  if [ -d $WORKINGDIR/software ]
  then
    echo "La capeta $SOFTWARE ya existe."
else
	mkdir -p $WORKINGDIR/software
    if [ $WORKINGDIR -eq 0 + ]
    then
      echo "SOFTWARE se ha creaco con éxito"
    else
      echo "Ups! Algo ha fallado al crear SOFTWARE"
    fi
  fi
done

Rscript "WGCNA_deps_R.R" $WORKINGDIR $CODEDIR $WORKINGDIR/software



