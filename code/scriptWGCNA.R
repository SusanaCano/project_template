workingDir <- "C:/Users/Mariana/Desktop/PROYECTO/BS"

setwd(workingDir)
getwd()
#
library(ggplot2)
library(factoextra)
library(NbClust)
library(cluster)
install.packages ("BiocManager")
BiocManager :: install ("DESeq2")
library(DESeq2)
library(WGCNA)

counts_df <- read.csv("GSE147507_RawReadCounts_Human.tsv", sep="\t", row.names=1)
counts_matrix <- as.matrix(counts_df)

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

sampleTree <- hclust(dist(normalized_counts_no_NaT), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file = "1_SampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(3,4,2,2))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#PONEMOS LA LINEA
abline(h = 790000, col = "red")##linea roja que separa para identificar los outliers
dev.off()




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
#normalizacion datos clusters

# Una vez tenemos el dataset dividido y normalizado ya podemos utilizar el WGCNA
library(WGCNA)
library(DCGL)
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


#SI SOLO TOMARAMOS EN EXPRESIONES CON UNA CONDICION O CLUSTER1 O CLUSTER 2 PARA NO JUNTAR:
expresion3<-rbind(exprs.1,exprs.2,exprs.3,exprs.4)
expresion4<-rbind(exprs.5,exprs.6,exprs.7,exprs.8)
expresion5<-rbind(exprs.6,exprs.7,exprs.8,exprs.9)
expresion6<-rbind(exprs.7,exprs.8,exprs.9,exprs.10)
expresion7<-rbind(exprs.8,exprs.9,exprs.10,exprs.11)
expresion8<-rbind(exprs.9,exprs.10,exprs.11,exprs.12)
expresion9<-rbind(exprs.10,exprs.11,exprs.12,exprs.13)


#Aquí lo que hacemos es observar para seleccionar loe genes coexpresados  deferenciales según el analisis
#de red de coexpresión generica ponderada (WGCNA)

#WGCNA.res1 <- WGCNA(exprs.1=expresion1, exprs.2=expresion2, power = 12, variant = "WGCNA")
#WGCNA.res1[1:6]
#Links <- qLinkfilter(expresion1, expresion2, 0.25)
#names(Links)
#Links$rth contienen los dos umbrales de correlación para ambas condiciones;
#Ambos Links$cor.filtered. mantienen las matrices de correlación filtradas para las condiciones A y B.
#Vemos los umbrales de correlación:
#umbral_ex1<-Links$rth.1
#vemos que tiene una correlacion de 0.803
#umbral_ex2<-Links$rth.2
# tiene una correlacion de 0.527
library(WGCNA)

WGCNA.res2 <- WGCNA(exprs.1=expresion3, exprs.2=expresion4, power = 12, variant = "WGCNA")
WGCNA.res2[1:24]
Links2 <- qLinkfilter(expresion3, expresion4, 0.25)
names(Links2)
umbral_ex3<-Links2$rth.1
#vemos que tiene una correlacion de 0.9946
umbral_ex4<-Links2$rth.2
#vemos que tiene una correlacion de 0.997


WGCNA.res3 <- WGCNA(exprs.1=expresion3, exprs.2=expresion5, power = 12, variant = "WGCNA")
WGCNA.res3[1:24]
Links3 <- qLinkfilter(expresion3, expresion5, 0.25)
names(Links3)
umbral_ex3<-Links3$rth.1
#vemos que tiene una correlacion de 0.9946
umbral_ex5<-Links3$rth.2
#vemos que tiene una correlacion de 0.9828

WGCNA.res4 <- WGCNA(exprs.1=expresion3, exprs.2=expresion6, power = 12, variant = "WGCNA")
WGCNA.res4[1:24]
Links4 <- qLinkfilter(expresion3, expresion6, 0.25)
names(Links4)
umbral_ex3<-Links4$rth.1
#vemos que tiene una correlacion de 0.9946
umbral_ex6<-Links4$rth.2
#Vamos que tiene una correlacion de 0.9303


WGCNA.res5 <- WGCNA(exprs.1=expresion3, exprs.2=expresion7, power = 12, variant = "WGCNA")
WGCNA.res5[1:24]
Links5 <- qLinkfilter(expresion3, expresion7, 0.25)
names(Links5)
umbral_ex3<-Links5$rth.1
#vemos que tiene una correlacion de 0.9946
umbral_ex7<-Links5$rth.2
#Vamos que tiene una correlacion de 0.990209


WGCNA.res6 <- WGCNA(exprs.1=expresion3, exprs.2=expresion8, power = 12, variant = "WGCNA")
WGCNA.res6[1:24]
Links6 <- qLinkfilter(expresion3, expresion8, 0.25)
names(Links6)
umbral_ex3<-Links6$rth.1
#vemos que tiene una correlacion de 0.9946
umbral_ex8<-Links6$rth.2
#Vamos que tiene una correlacion de 0.9964594

WGCNA.res7 <- WGCNA(exprs.1=expresion3, exprs.2=expresion9, power = 12, variant = "WGCNA")
WGCNA.res7[1:24]
Links7 <- qLinkfilter(expresion3, expresion9, 0.25)
names(Links7)
umbral_ex3<-Links7$rth.1
#vemos que tiene una correlacion de 0.9946
umbral_ex9<-Links7$rth.2
#Vamos que tiene una correlacion de 0.391231


####POR LO QUE PODEMOS LLEGAR A LA CONCLUSIÓN QUE LOS QUE DEL CLUSTER 1 DEBEMOS USAR LA EXPRESION 3 
####Y DEL CLUSTER 2 DEBERIAMOS USAR LA EXPRESION 4 QUE ES LA QUE MAYOR CORRELACION TIENE

datos_expresiones<-rbind(expresion3,expresion4,expresion8,expresion7,expresion5,expresion6)

####network constr#######
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers <- c(1:30)

# Call the network topology analysis function
sft <- pickSoftThreshold(datos_expresiones, powerVector = powers, verbose = 5)
# Plot the results:
#sizeGrWindow(9, 5)
pdf(file="2_Scale independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower


softPower <-sft$powerEstimate
adjacency <- adjacency(datos_expresiones, power = softPower)

memory.limit(size = 35000)#para problema de memoria
##### Turn adjacency into topological overlap
TOM <- TOMsimilarity(adjacency);
dissTOM <- 1-TOM
rm(adjacency)#eliminar adjagency.
# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "average")
#hc <- eclust(datos,FUNcluster = "hclust", k = i, hc_metric = d , hc_method = l)
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file="3_Gene clustering on TOM-based dissimilarity.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize <- 6
# Module identification using dynamic tree cut:
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
#como vemos mayoritariamente va a haber gris porque hay 21622
# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="4_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

##CALCULAMOS EL GEN PROPIO
MEList <- moduleEigengenes(datos_expresiones, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1-cor(MEs);
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file="5_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# Call an automatic merging function
merge <- mergeCloseModules(datos_expresiones, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules:
mergedMEs <-merge$newMEs
#sizeGrWindow(12, 9)
pdf(file="6_Merged dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
#Observar los colores cambiados

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs


#####
nSamples<-nrow(datos_expresiones)

### COVID 


# names (colors) of the modules
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datos_expresiones, MEs, use = "p"))
COVIDPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) <- paste("COVID", modNames, sep="")
names(COVIDPvalue) <- paste("p.COVID", modNames, sep="")


#####EXPORT

names(datos_expresiones)
probes <- names(datos_expresiones)


#################exprot COVID ############### 

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

write.table(geneInfo, file = "7_COVID.xls",sep="\t",row.names=F)
#Vemos el valor de los genes en las diferentes muestras

########VISUALIZAR RED GENES#################


nGenes <- ncol(datos_expresiones)
nSamples <- nrow(datos_expresiones)


# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#TOM = Matriz de superposición topológica (aclaracion)
plotTOM <- dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) <- NA
# Call the plot function

#sizeGrWindow(9,9)
memory.limit(size = 80000)#para problema de memoria
pdf(file="8_Network heatmap plot_all gene.pdf",width=9, height=9,compress = FALSE)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

nSelect <- 400
# For reproducibility, we set the random seed
set.seed(10)
select <- sample(nGenes, size = nSelect,replace = TRUE)
selectTOM <- dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]

# Open a graphical window
#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss <- selectTOM^7
diag(plotDiss) <- NA

pdf(file="9_Network heatmap plot_selected genes.pdf",width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()


#########Visualizing the gene network of eigengenes#############


#sizeGrWindow(5,7.5)

#Dendrograma de genes propios y mapa de calor de adyacencia de genes propios
pdf(file="10_Eigengene dendrogram and Eigengene adjacency heatmap.pdf", width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

#or devide into two parts
# Plot the dendrogram
#sizeGrWindow(6,6);
pdf(file="11_Eigengene dendrogram_2.pdf",width=6, height=6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

#Mapa de calor de adyacencia de genes propios

pdf(file="12_Eigengene adjacency heatmap_2.pdf",width=6, height=6)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()


###########################Exporting to Cytoscape all one by one ##########################
# Select each module

for (mod in 1:nrow(table(moduleColors)))
{
  
  modules <- names(table(moduleColors))[mod]
  # Select module probes
  probes <- names(datos_expresiones)
  inModule <- (moduleColors == modules)
  modProbes <- probes[inModule]
  modGenes <- modProbes
  # Select the corresponding Topological Overlap
  modTOM <- TOM[inModule, inModule]
  
  dimnames(modTOM) <- list(modProbes, modProbes)
  # Export the network into edge and node list files Cytoscape can read
  cyt <- exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", modules , ".txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", modules, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}

BiocManager ::install ("coexnet")
library(coexnet)
library(tidyverse)

# Building the network
dataC<-geneInfo0[,-c(1:2)]
memory.limit(size = 120000)#hay que aumentar el espacio de memoria

cor_pearson <- createNet(expData = dataC,threshold = 0.99,method = "correlation")
plot(cor_pearson)
jpeg("cor_pearson.jpeg")
plot(cor_pearson)
dev.off()

