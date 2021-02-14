
workingDir <- "C:/Users/Mariana/Desktop/PROYECTO/BS/datos nuevos"
setwd(workingDir)
getwd()
#
install.packages ("BiocManager")
BiocManager :: install ("DESeq2")
BiocManager ::install ("coexnet")

library(coexnet)
library(tidyverse)
library(ggplot2)
library(factoextra)
library(NbClust)
library(cluster)
library(dplyr)
library(DESeq2)
library(WGCNA)
library(DCGL)

library(imager)#para visualizar la parte de la neurona


#####PREPROCESAMIENTO DE DATOS Y CLUSTERING#######
##PREPROCESAMIENTO DE DATOS

#Cargamos los datos
counts_df <- read.csv("GSE147507_RawReadCounts_Human.tsv", sep="\t", row.names=1)
counts_df[counts_df==0]<-NA
counts_dfT<-as.data.frame(t(counts_df))
#Estas serían las columnas a borrar
cols_borrar <- which(colMeans(is.na(counts_dfT)) >=0.51)
#Borramos donde todas las filas sean NA
counts_dfT <- counts_dfT[,-cols_borrar]
#Transformamos de nuevo los NA en 0
funcionTo0 <- function(x){
  ifelse(is.na(x), 0, x)
}
datos<- sapply(counts_dfT,funcionTo0)
filas <- rownames(counts_dfT)

###DATOS SIN COLUMNAS QUE SON ENTERAS 0
datos<- data.frame(datos,row.names=filas)


counts_df<-t(datos)
counts_matrix <- as.matrix(counts_df)


#Hacemos clustering usando la función pam() para ver los grupos
kmedioids2 <- pam(datos,2)


#Guardamos en una nueva variable, un array de los clusters recien creados
counts_clustering<-kmedioids2[["clustering"]]

#Creo un dataframe con el conjunto de datos agrupado
meta <- data.frame(sample=colnames(counts_matrix), sample_type = counts_clustering)

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_matrix, colData = meta, design = ~ sample_type)
dds <- DESeq2::DESeq(dds)

#Normalizamos los valores 
normalized_counts <- as.data.frame(DESeq2::counts(dds, normalized=TRUE)) 

#Sustituyo valores perdidos (0) por la media de la columna a la que pertenece y creamos un nuevo dataframe
transformToMean <- function(x){
  ifelse(x == 0, round(mean(x), 2), x)
}

normalized_counts_no_Na<-data.frame(sapply(normalized_counts,transformToMean))

#Reorganizamos los nombres de las filas
rownames(normalized_counts_no_Na) <- rownames(counts_df)
normalized_counts_no_NaT<-as.data.frame(t(normalized_counts_no_Na))

##MASFACIL##
datos_normalizados<-normalized_counts_no_NaT

##CLUSTERING

#Aplicamos clustering usando la funcion hclust() con el linkage average
sampleTree <- hclust(dist(datos_normalizados), method = "average")
setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
#Creamos un pdf del arbol que muestra el clustering realizado
pdf(file = "1_SampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(3,4,2,2))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 100, col = "red")#linea roja que separa para identificar los outliers
dev.off()
#setwd(workingDir)
#getwd()



#Uno las columnas que indican el cluster al que pertenecen
datos_normalizados_cluster<-cbind(Cluster=as.factor(counts_clustering),datos_normalizados)
data<-split(datos_normalizados_cluster, f = datos_normalizados_cluster$Cluster)
cluster1data<-as.data.frame(data[["1"]])
cluster2data<-as.data.frame(data[["2"]])

# Elimino la variable cluster, ya que estan los datos divididos y no sirven para nada
cluster1data <- select(cluster1data, -Cluster)
cluster2data <- select(cluster2data, -Cluster)


#Subdividimos los dos conjuntos de datos en un mismo numero de objetos para poder aplicar el análisis de coexpresion de genes.  

#cluster1
exprs.1<-cluster1data[1:6,]
exprs.2<-cluster1data[7:12,]
exprs.3<-cluster1data[13:18,]
exprs.4<-cluster1data[19:24,]

#Cluster2
exprs.5<-cluster2data[1:6,]
exprs.6<-cluster2data[7:12,]
exprs.7<-cluster2data[13:18,]
exprs.8<-cluster2data[19:24,]
exprs.9<-cluster2data[25:30,]
exprs.10<-cluster2data[31:36,]
exprs.11<-cluster2data[37:42,]
exprs.12<-cluster2data[43:48,]
exprs.13<-cluster2data[49:54,]


expresion3<-rbind(exprs.1,exprs.2,exprs.3,exprs.4)
expresion4<-rbind(exprs.5,exprs.6,exprs.7,exprs.8)
expresion5<-rbind(exprs.6,exprs.7,exprs.8,exprs.9)
expresion6<-rbind(exprs.7,exprs.8,exprs.9,exprs.10)
expresion7<-rbind(exprs.8,exprs.9,exprs.10,exprs.11)
expresion8<-rbind(exprs.9,exprs.10,exprs.11,exprs.12)
expresion9<-rbind(exprs.10,exprs.11,exprs.12,exprs.13)


