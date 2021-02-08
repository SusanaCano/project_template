setwd("C:/Users/JuanS/Documents/Proyecto BS")

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

#
library(factoextra)
library(NbClust)
library(cluster)
install.packages ("BiocManager")
BiocManager :: install ("DESeq2")
library('DESeq2')


counts_df <- read.csv("GSE147507_RawReadCounts_Human.tsv", sep="\t", row.names=1)
counts_matrix <- as.matrix(counts_df)

# YOU NEED TO CHANGE THIS TO MATCH THE REAL EXPERIMENTAL DESIGN
# Voy a hacer una clasificacion ya que lo de 39 y 39 es un ejemplo del profesor dividiendo los datos 
#por la mitad sin hacer ningun tipo de clusterinh
#DIVISION DEL PROFESOR) sample_grouping <- as.factor(c(rep("a", 39), rep("b", 39)))
#Hago clustering para ver los dos grupos
counts_dfT<-as.data.frame(t(counts_df))
kmedioids2 <- pam(counts_dfT,2)
#CLUSTERING No ejecutar tarda mucho
fviz_nbclust(normalized_counts, FUNcluster=hcut, k.max = 5, method = "wss")
#Guardo un array de los clusters
counts_clustering<-kmedioids2[["clustering"]]
#Creo un dataframe con los clusters y genes
meta <- data.frame(sample=colnames(counts_matrix), sample_type = counts_clustering)
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_matrix, colData = meta, design = ~ sample_type)
dds <- DESeq2::DESeq(dds)
normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized=TRUE)) # Getting normalized values
#Sustituyo valores perdidos por la media
transformToMean <- function(x){
  ifelse(x == 0, round(mean(x), 2), x)
}
normalized_counts_no_Na=data.frame(sapply(normalized_counts,transformToMean))
#Pongo nombres a las filas con cada gen
rownames(normalized_counts_no_Na) <- rownames(counts_df)
#Hago la transpuesta
normalized_counts_no_NaT<-as.data.frame(t(normalized_counts_no_Na))
#Uno los clusters a cada gen
normalized_counts_no_Na_cluster<-cbind(Cluster=as.factor(counts_clustering),normalized_counts_no_NaT)
# Divido el dataset dependiendo del cluster
data<-split(normalized_counts_no_Na_cluster, f = normalized_counts_no_Na_cluster$Cluster)
cluster1data<-as.data.frame(data[["1"]])
cluster2data<-as.data.frame(data[["2"]])
# Elimino la variable cluster, ya que estan los datos divididos y no sirven para nada
library(dplyr)
cluster1data <- select(cluster1data, -Cluster)
cluster2data <- select(cluster2data, -Cluster)
# Una vez tenemos el dataset dividido y normalizado ya podemos utilizar el WGCNA
library(WGCNA)





