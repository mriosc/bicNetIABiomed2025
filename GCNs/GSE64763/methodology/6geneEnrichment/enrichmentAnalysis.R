# Analisis de enriquecimiento funcional
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DO.db")
BiocManager::install("enrichplot", force = TRUE)
install.packages("plyr")

# Establecer directorio de trabajo y cargar las

library(BiocManager)
library(org.Hs.eg.db)
library("hgu133plus2.db")
library("affyPLM")
library(clusterProfiler)
library(enrichplot)
library(DO.db)
library(ggplot2)
library("enrichplot")
library("plyr")
library(DOSE)

setwd("C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/methodology/7geneEnrichment")

carpeta <- "7geneEnrichment//"
archivos <- list.files(carpeta)

archivos_txt <- archivos[grep("\\.txt$", archivos)]

# Analisis para GO
for (archivo in archivos_txt) {
  genelist <- readLines(paste0(carpeta, archivo))
  entrez_ids <- mapIds(hgu133plus2.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
  go <- enrichGO(entrez_ids,  OrgDb = "hgu133plus2.db", ont = "all", pAdjustMethod = "fdr")
  if(length(go@result$ONTOLOGY) != 0){
    write.csv(go, paste0(carpeta, archivo,"_background_all_BH.csv"), row.names = FALSE)
    grafico <- dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    ggsave(filename = paste0(carpeta, archivo,"_background_all_BH.svg"), plot = grafico, device = "svg")
  }
}

# Analisis para GO
#for (archivo in archivos_txt) {
 # genelist <- readLines(paste0(carpeta,archivo))
  #entrez_ids <- mapIds(org.Hs.eg.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
  #go <- enrichGO(entrez_ids,  OrgDb = "org.Hs.eg.db", ont = "all", pAdjustMethod = "fdr")
  #if(length(go@result$ONTOLOGY) != 0){
   # write.csv(go, paste0(carpeta,archivo,"all_BH.csv"), row.names = FALSE)
    #grafico <- dotplot(go, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
    #ggsave(filename = paste0(carpeta,archivo,"all_BH.svg"), plot = grafico, device = "svg")
  #}
#}

# Analisis para KEGG
for (archivo in archivos_txt) {
  genelist <- readLines(paste0(carpeta, archivo))
  entrez_ids <- mapIds(hgu133plus2.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
  kegg <- enrichKEGG(entrez_ids, pAdjustMethod = "fdr")
  if(!is.null(kegg)){
    if(any(kegg@result$p.adjust <0.05)){
      write.csv(kegg, paste0(carpeta, archivo,"_background_KEGG_BH.csv"), row.names = FALSE)
      grafico <- dotplot(kegg) + facet_grid(scale="free")
      ggsave(filename = paste0(carpeta, archivo,"_background_KEGG_BH.svg"), plot = grafico, device = "svg")
    }
  }
}


#for (archivo in archivos_txt) {
#  genelist <- readLines(paste0(carpeta,archivo))
#  entrez_ids <- mapIds(org.Hs.eg.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
#  kegg <- enrichKEGG(entrez_ids, pAdjustMethod = "fdr")
#  if(!is.null(kegg)){
#    if(any(kegg@result$p.adjust <0.05)){
#      write.csv(kegg, paste0(carpeta,archivo,"KEGG_BH.csv"), row.names = FALSE)
#      grafico <- dotplot(kegg) + facet_grid(scale="free")
#      ggsave(filename = paste0(carpeta,archivo,"KEGG_BH.svg"), plot = grafico, device = "svg")
#    }
#  }
#}

# unidades

genelist3 <- readLines("Cluster_3.txt")
genelist7 <- readLines("Cluster_7.txt")

# PARA EL CLUSTER 3
entrez_ids3 <- mapIds(hgu133plus2.db, keys=genelist3, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids3,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  #write.csv(go, "cluster_normal_prostate//Cluster_3_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_3_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids3, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_3_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_3_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 7
entrez_ids7 <- mapIds(hgu133plus2.db, keys=genelist7, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids7,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  #write.csv(go, "cluster_normal_prostate//Cluster_2_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_7_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids5, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# cambiar graficos finales
resultados_lista <- list()

# Loop sobre cada archivo
for (archivo in archivos_txt) {
  genelist5 <- readLines(paste0(carpeta, archivo))
  entrez_ids5 <- mapIds(org.Hs.eg.db, keys=genelist, column='ENTREZID', keytype='SYMBOL')
  go <- enrichGO(entrez_ids5,  OrgDb = "org.Hs.eg.db", ont = "all", pAdjustMethod = "fdr")
  
  if(length(go@result$ONTOLOGY) != 0){
    # Agregar la información del archivo actual a la lista
    resultados_lista[[archivo]] <- go
  }
}

# Convertir todos los elementos en resultados_lista a data.frames
resultados_lista <- lapply(resultados_lista, as.data.frame)

# Combinar los resultados de todos los archivos
resultados_combinados <- rbind.fill(resultados_lista)

# Agregar la información de Shape al objeto de enriquecimiento
resultados_combinados$Shape <- factor(resultados_combinados$ONTOLOGY)

# Crear el gráfico utilizando dotplot
grafico <- dotplot(resultados_combinados, 
                   x="GeneRatio", 
                   color="Shape",
                   legend="bottom",
                   pch=16, 
                   cex=1.5)

# Guardar la gráfica en formato SVG
ggsave(filename = "resultados_combinados.svg", plot = grafico, device = "svg")




