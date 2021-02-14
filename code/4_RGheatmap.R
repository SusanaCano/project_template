
args = commandArgs(trailingOnly = TRUE)
wd <- args[1]
software_deps <- args[2]
results <- args[3]

#results <- "D:/uni/2020-21/BiologiaSistemas/Practicas/Aurelio/trabajoFinal/Git/project_template_BS_SARS-CoV2/results"
#wd <- "D:/uni/2020-21/BiologiaSistemas/Practicas/Aurelio/trabajoFinal/Git/project_template_BS_SARS-CoV2/code"
setwd(wd)

#options(stringsAsFactors = FALSE)

source(paste(wd, "3_RedTopologica.R", sep = "/"))

########VISUALIZAR RED GENES MAPA DE CALOR #################


nGenes <- ncol(datos_expresiones)
nSamples <- nrow(datos_expresiones)


# Transforma dissTOM con el power para hacer que las conexiones moderadamente fuertes sean más visibles en el mapa de calor

plotTOM <- dissTOM^7
#Esto se pone para que quede un plot más estético
diag(plotTOM) <- NA
memory.limit(size = 80000)

#Creamos un pdf que nos muestre el mapa de calor
setwd(results)
pdf(file = "8_Network_heatmap_plot_all_gene.pdf", width = 9, height = 9,compress = FALSE)
TOMplot(as.data.frame(plotTOM), geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()
setwd(wd)

nSelect <- 400
set.seed(10)
select <- sample(nGenes, size = nSelect, replace = TRUE)
selectTOM <- dissTOM[select, select]

# No existe una forma sencilla de restringir un árbol de agrupamiento a un subconjunto de genes, por lo que debemos volver a agrupar.
selectTree <- hclust(as.dist(selectTOM), method = "average")
selectColors <- moduleColors[select]


# Elevar la disimilitud a una potencia, digamos 7, hace que la trama sea más informativa al cambiar la paleta de colores; 
#establecer la diagonal en NA también mejora la claridad de la trama
plotDiss <- selectTOM^7
diag(plotDiss) <- NA

setwd(results)
pdf(file = "9_Network_heatmap_plot_selected_genes.pdf", width = 9, height = 9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()
setwd(wd)

#Creamos un pdf donde nos muestr el dendrograma de genes propios y mapa de calor de adyacencia de genes propios
setwd(results)
pdf(file = "10_Eigengene_dendrogram_and_Eigengene_adjacency_heatmap.pdf", width = 5, height = 7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0, 4, 1, 2), marHeatmap = c(3, 4, 1, 2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()
setwd(wd)

#Aunque tambien lo podemos dividir en 2 partes:
setwd(results)
pdf(file = "11_Eigengen_dendrogram_2.pdf", width = 6, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0, 4, 2, 0), plotHeatmaps = FALSE)
dev.off()
setwd(wd)

#Mapa de calor de adyacencia de genes propios
setwd(results)
pdf(file = "12_Eigengene_adjacency_heatmap_2.pdf", width = 6, height = 6)
# Hace la matriz del mapa de calor (nota: este gráfico sobrescribirá el gráfico de dendrograma)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()
setwd(wd)