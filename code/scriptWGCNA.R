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
library(DCGL)
#Creamos las diferentes expresiones teniendo cuidado porque debemos tener el mismo numero  en las filas asi dividimos en 6 objetos con 6 variables
#Y ademá tenemos que ir recorriendo esas variables.
#cluster1
exprs.1<-cluster1data[1:6,1:6]
exprs.2<-cluster1data[7:12,7:12]
exprs.3<-cluster1data[13:18,13:18]
exprs.4<-cluster1data[19:24,19:24]

#Cluster2
exprs.5<-cluster2data[1:6,1:6]
exprs.6<-cluster2data[7:12,7:12]
exprs.7<-cluster2data[13:18,13:18]
exprs.8<-cluster2data[19:24,19:24]
exprs.9<-cluster2data[25:30,25:30]
exprs.10<-cluster2data[31:36,31:36]
exprs.11<-cluster2data[37:42,37:42]
exprs.12<-cluster2data[43:48,43:48]
exprs.13<-cluster2data[49:54,49:54]

#UNIMOS LA MITAD DE LAS EXPRESIONES ANTERIORES Y DIVIDIMOS EN 2
expresion1<-cbind(exprs.1,exprs.2,exprs.3,exprs.4,exprs.5,exprs.6,exprs.7[13:15,])
expresion2<-cbind(exprs.7[16:18,],exprs.8,exprs.9,exprs.10,exprs.11,exprs.12,exprs.13)

#Aquí lo que hacemos es observar para seleccionar loe genes coexpresados  deferenciales según el analisis
#de red de coexpresión generica ponderada (WGCNA)
WGCNA.res1 <- WGCNA(exprs.1=expresion1, exprs.2=expresion2, power = 12, variant = "WGCNA")
WGCNA.res1[1:6]

#Vemos los enlaces y sus umbrales de correlacion para usarlo luego
Links <- qLinkfilter(expresion1, expresion2, 0.25)
names(Links)
#Links$rth contienen los dos umbrales de correlación para ambas condiciones;
#Ambos Links$cor.filtered. mantienen las matrices de correlación filtradas para las condiciones A y B.
#Vemos los umbrales de correlación:
umbral_ex1<-Links$rth.1
#tiene una correlacion de 0.803
umbral_ex2<-Links$rth.2
# tiene una correlacion de 0.527




