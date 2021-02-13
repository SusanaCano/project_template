########VISUALIZAR RED GENES MAPA DE CALOR #################


nGenes <- ncol(datos_expresiones)
nSamples <- nrow(datos_expresiones)


# Transforma dissTOM con el power para hacer que las conexiones moderadamente fuertes sean más visibles en el mapa de calor

plotTOM <- dissTOM^7
#Esto se pone para que quede un plot más estético
diag(plotTOM) <- NA
memory.limit(size = 80000)

#Creamos un pdf que nos muestre el mapa de calor
pdf(file="8_Network heatmap plot_all gene.pdf",width=9, height=9,compress = FALSE)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

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