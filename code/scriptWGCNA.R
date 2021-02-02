#Cargar datos
datos <- read.table(file='GSE147507_RawReadCounts_Human.tsv')

#Instalar paquetes
install.packages ("BiocManager")
BiocManager :: install ("WGCNA")
install.packages ("dplyr")
#Cargar paquetes
library(dplyr)
library(WGCNA)

#Sustituir datos perdidos por media
datosWGCNA <- read.table(file='GSE147507_RawReadCounts_Human.tsv')
transformToMean <- function(x){
  ifelse(x == 0, round(mean(x), 2), x)
}
datosNoNa=data.frame(sapply(datosWGCNA,transformToMean))
datosNoNaGenes<-cbind(Genes=rownames(datosWGCNA),datosNoNa)
datosNoNaT  = as.data.frame (t ( datosNoNaGenes [, - 1 ]))

#Filtrado 
mean.filt = 0.5
n = nrow ( datosNoNaT )
datosNoNaT[n+1,]=apply(datosNoNaT[c(1:nrow(datosNoNaT)),],2,mean)
datosNoNaT=datosNoNaT[1:n,datosNoNaT[n+1,]>mean.filt]
datosFiltradosGenes=t(datosNoNaT)
datosFiltradosGenes=data.frame(rownames(datosFiltradosGenes),datosFiltradosGenes)
datosFiltradosGenes<-select(datosFiltradosGenes,-rownames.datosFiltradosGenes.)
write.table(datosFiltradosGenes,file="data.filt.xls",row.names=rownames(datosWGCNA),col.names=T,quote=FALSE,sep="\t")

#Leer nuevos datos
data.filt <- read.table(file='data.filt.xls')
#Continuar...
