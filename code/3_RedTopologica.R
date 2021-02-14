####CONSTRUCCIÓN RED TOPOLÓGICA#####
enableWGCNAThreads()

# Escogemos un conjunto de soft-thresholding powers
powers <- c(1:30)

# Función de análisis de red topológica
sft <- pickSoftThreshold(datos_expresiones, powerVector = powers, verbose = 5)

setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
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
setwd(workingDir)
getwd()


#Escoger el softPower
softPower <-sft$powerEstimate

#Hallamos la adyacencia
adjacency <- adjacency(datos_expresiones, power = softPower)

# Convertimos la adyacencia en TOM ( Matriz Topológica de Superposición)
TOM<-TOMsimilarity(adjacency)
dissTOM <- 1-TOM


# Aplicamos clustering usando la funcion hclust() con el linkage average
geneTree <- hclust(as.dist(dissTOM), method = "average")

setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
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

setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
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

setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
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

setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
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

setwd("C:/Users/Mariana/Documents/GitHub/project_template_BS_SARS-CoV2/results")
getwd()
write.table(geneInfo, file = "7_COVID.xls",sep="\t",row.names=F)
#Vemos el valor de los genes en las diferentes muestras
setwd(workingDir)
getwd()
