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
#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
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


##### ANALISIS DE COEXPRESIÓN USANDO EL PAQUETE WGCNA #####

#Usamos la funcion WGCNA para elegir los genes de mayor correlación

WGCNA.res2 <- WGCNA(exprs.1=expresion3, exprs.2=expresion4, power = 12, variant = "WGCNA")
WGCNA.res2[1:24]
Links2 <- qLinkfilter(expresion3, expresion4, 0.25)
names(Links2)
umbral_ex3<-Links2$rth.1
#vemos que tiene una correlacion de 0.994
umbral_ex4<-Links2$rth.2
#vemos que tiene una correlacion de 0.997


WGCNA.res3 <- WGCNA(exprs.1=expresion3, exprs.2=expresion5, power = 12, variant = "WGCNA")
WGCNA.res3[1:24]
Links3 <- qLinkfilter(expresion3, expresion5, 0.25)
names(Links3)
umbral_ex3<-Links3$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex5<-Links3$rth.2
#vemos que tiene una correlacion de 0.9833

WGCNA.res4 <- WGCNA(exprs.1=expresion3, exprs.2=expresion6, power = 12, variant = "WGCNA")
WGCNA.res4[1:24]
Links4 <- qLinkfilter(expresion3, expresion6, 0.25)
names(Links4)
umbral_ex3<-Links4$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex6<-Links4$rth.2
#Vamos que tiene una correlacion de 0.9291


WGCNA.res5 <- WGCNA(exprs.1=expresion3, exprs.2=expresion7, power = 12, variant = "WGCNA")
WGCNA.res5[1:24]
Links5 <- qLinkfilter(expresion3, expresion7, 0.25)
names(Links5)
umbral_ex3<-Links5$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex7<-Links5$rth.2
#Vamos que tiene una correlacion de 0.990


WGCNA.res6 <- WGCNA(exprs.1=expresion3, exprs.2=expresion8, power = 12, variant = "WGCNA")
WGCNA.res6[1:24]
Links6 <- qLinkfilter(expresion3, expresion8, 0.25)
names(Links6)
umbral_ex3<-Links6$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex8<-Links6$rth.2
#Vamos que tiene una correlacion de 0.996

WGCNA.res7 <- WGCNA(exprs.1=expresion3, exprs.2=expresion9, power = 12, variant = "WGCNA")
WGCNA.res7[1:24]
Links7 <- qLinkfilter(expresion3, expresion9, 0.25)
names(Links7)
umbral_ex3<-Links7$rth.1
#vemos que tiene una correlacion de 0.99441
umbral_ex9<-Links7$rth.2
#Vamos que tiene una correlacion de 0.0569

#no ponemos la expresion 9 porque tiene una correlacion de 0.05
datos_expresiones<-rbind(expresion3,expresion4,expresion7,expresion8)


####CONSTRUCCIÓN RED TOPOLÓGICA#####
enableWGCNAThreads()

# Escogemos un conjunto de soft-thresholding powers
powers <- c(1:30)

# Función de análisis de red topológica
sft <- pickSoftThreshold(datos_expresiones, powerVector = powers, verbose = 5)

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
#Creamos un pdf donde se vea la red topológica en libre escala
pdf(file="2_Scale independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")#esta línea corresponde a usar un corte R ^ 2 de h
#Conectividad media en función de la potencia de soft-thresholding
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#setwd(workingDir)
#getwd()


#Escoger el softPower
softPower <-sft$powerEstimate

#Hallamos la adyacencia
adjacency <- adjacency(datos_expresiones, power = softPower)

# Convertimos la adyacencia en TOM ( Matriz Topológica de Superposición)
TOM<-TOMsimilarity(adjacency)
dissTOM <- 1-TOM


# Aplicamos clustering usando la funcion hclust() con el linkage average
geneTree <- hclust(as.dist(dissTOM), method = "average")

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
#Creamos un pdf donde nos muestr ese Clusterin basado en la no similitud
pdf(file="3_Gene clustering on TOM-based dissimilarity.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
setwd(workingDir)
getwd()

#Según lo visto en los pdfs anteriores queremos que el tamaño minimo del módulo de esos clusters sean 6 (un numero relativamente alto conforme a las particiones posibles)
minModuleSize <- 5

#Identificación del módulo mediante corte de árbol dinámico:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#como vemos mayoritariamente va a haber gris porque hay 21622

# Convertir etiquetas numéricas en colores
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
#Hacemos un pdf donde se haga el dendrograma y los colores debajo
pdf(file="4_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
setwd(workingDir)
getwd()

#gen propio
MEList <- moduleEigengenes(datos_expresiones, colors = dynamicColors)
MEs <- MEList$eigengenes

# Calcular la disimilitud de los genes propios del módulo
MEDiss <- 1-cor(MEs)

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
#Creamos un pdf  que muestra el clustering realizado según el árbol de modulos de los genes
pdf(file="5_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")# h= va a ser la línea que corte el dendrograma
dev.off()
setwd(workingDir)
getwd()


#Llamada a una funcion merge que lo que nos permite es unir los datos expresiones con los colores antes hallados
merge <- mergeCloseModules(datos_expresiones, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# Colores que hay en la función merge anteruir
mergedColors <- merge$colors

# Genes propios de los nuevos módulos fusionados:
mergedMEs <-merge$newMEs

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
#Hacemos un denodrograma y color con el árbol hallado anteriormente y los colores tanto dinamicos como en función merge
pdf(file="6_Merged dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
setwd(workingDir)
getwd()


moduleColors <- mergedColors
#Construimos etiquetas numéricas que correspondan a los colores
colorOrder <- c("grey", standardColors(50))
#Ponemos el grey porque es el que mas hay
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs


nSamples<-nrow(datos_expresiones)
###p-valor COVID ######


# names (colors) of the modules
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datos_expresiones, MEs, use = "p"))
COVIDPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("COVID", modNames, sep="")
names(COVIDPvalue) <- paste("p.COVID", modNames, sep="")


#####EXPORT####

names(datos_expresiones)
probes <- names(datos_expresiones)

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)


for (mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,mod],
                         COVIDPvalue[, mod])
  names(geneInfo0) <- c(oldNames,names(geneModuleMembership)[mod],
                       names(COVIDPvalue)[mod])
}
geneOrder <-order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder, ]

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
write.table(geneInfo, file = "7_COVID.xls",sep="\t",row.names=F)
#Vemos el valor de los genes en las diferentes muestras
setwd(workingDir)
getwd()


########VISUALIZAR RED GENES MAPA DE CALOR #################
nGenes <- ncol(datos_expresiones)
nSamples <- nrow(datos_expresiones)


# Transforma dissTOM con el power para hacer que las conexiones moderadamente fuertes sean más visibles en el mapa de calor

plotTOM <- dissTOM^2
#Esto se pone para que quede un plot más estético
diag(plotTOM) <- NA
memory.limit(size = 9999999)


#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
#getwd()
#Creamos un pdf que nos muestre el mapa de calor
pdf(file="8_Network heatmap plot_all gene.pdf",width=12, height=12,compress = FALSE)
TOMplot(as.data.frame(plotTOM), geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()
setwd(workingDir)
getwd()

nSelect <- 400
set.seed(10)
select <- sample(nGenes, size = nSelect,replace = TRUE)
selectTOM <- dissTOM[select, select]

# No existe una forma sencilla de restringir un árbol de agrupamiento a un subconjunto de genes, por lo que debemos volver a agrupar.
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]


# Elevar la disimilitud a una potencia, digamos 7, hace que la trama sea más informativa al cambiar la paleta de colores; 
#establecer la diagonal en NA también mejora la claridad de la trama
plotDiss <- selectTOM^7
diag(plotDiss) <- NA

#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
pdf(file="9_Network heatmap plot_selected genes.pdf",width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()



#Creamos un pdf donde nos muestr el dendrograma de genes propios y mapa de calor de adyacencia de genes propios
pdf(file="10_Eigengene dendrogram and Eigengene adjacency heatmap.pdf", width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()



#Aunque tambien lo podemos dividir en 2 partes:
pdf(file="11_Eigengene dendrogram_2.pdf",width=6, height=6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

#Mapa de calor de adyacencia de genes propios

pdf(file="12_Eigengene adjacency heatmap_2.pdf",width=6, height=6)
# Hace la matriz del mapa de calor (nota: este gráfico sobrescribirá el gráfico de dendrograma)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
setwd(workingDir)
getwd()

#######################CYTOSCAPE##########################
# Hacemos una funcion que seleccione cada uno de los modulos de color
#memory.limit(size = 120000)
#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
for (mod in 1:nrow(table(moduleColors)))
{
  
  modules <- names(table(moduleColors))[mod]
  #Seleccionamos los modulos de las sondas(probes)
  probes <- names(datos_expresiones)
  inModule <- (moduleColors == modules)
  modProbes <- probes[inModule]
  modGenes <- modProbes
  # Seleccionamos la superposición topológica correspondiente
  modTOM <- TOM[inModule, inModule]
  
  dimnames(modTOM) <- list(modProbes, modProbes)
  #Exportamos la red a archivos de lista de nodos y de borde que Cytoscape puede leer 
  cyt <- exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", modules , ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}



# Constuimos la red

cor_pearson1 <- createNet(expData = datos_expresiones,threshold = 1,method = "correlation")
plot(cor_pearson1)
jpeg("cor_pearson.jpeg")
plot(cor_pearson1)
dev.off()
setwd(workingDir)
getwd()


#para hacer un boxplot

complete <- cofVar(datos_expresiones)
#setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
pdf(file="13_Boxplot.pdf",width=6, height=6)
boxplot(complete, main="Boxplot groups")
dev.off()
setwd(workingDir)
getwd()
vb<-boxplot(dataC1)
vb$out
#vemos así los valor atípicos que son 1000


