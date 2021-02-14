args = commandArgs(trailingOnly = TRUE)
wd <- args[2]
software_deps <- args[3]
results <- args[4]


setwd(wd)

source(paste(wd, "4_RGheatmap.R", sep = "/"))

#######################CYTOSCAPE##########################

# Hacemos una funcion que seleccione cada uno de los modulos de color
for (mod in 1:nrow(table (moduleColors)))
{
  
  modules <- names(table (moduleColors))[mod]
  #Seleccionamos los modulos de las sondas(probes)
  probes <- names (datos_expresiones)
  inModule <- (moduleColors == modules)
  modProbes <- probes[inModule]
  modGenes <- modProbes
  # Seleccionamos la superposición topológica correspondiente
  modTOM <- TOM[inModule, inModule]
  
  dimnames(modTOM) <- list(modProbes, modProbes)
  #Exportamos la red a archivos de lista de nodos y de borde que Cytoscape puede leer 
  cyt <- exportNetworkToCytoscape(modTOM,
                                  edgeFile = paste("CytoscapeInput-edges-", modules , ".txt", sep = ""),
                                  nodeFile = paste("CytoscapeInput-nodes-", modules, ".txt", sep = ""),
                                  weighted = TRUE,
                                  threshold = 0.02,
                                  nodeNames = modProbes,
                                  altNodeNames = modGenes,
                                  nodeAttr = moduleColors[inModule])
}



# Constuimos la red
setwd(results)

memory.limit(size = 120000)

cor_pearson <- createNet(expData = datos_expresiones, threshold = 1, method = "correlation")
plot(cor_pearson)
jpeg("cor_pearson.jpeg")
plot(cor_pearson)
dev.off()
setwd(wd)