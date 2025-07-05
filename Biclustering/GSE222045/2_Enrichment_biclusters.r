if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
BiocManager::install("org.Hs.eg.db", force=TRUE)
BiocManager::install("hgu133plus2.db", force=TRUE)
BiocManager::install("DOSE",force=TRUE)
BiocManager::install("enrichplot",force=TRUE)
BiocManager::install("clusterProfiler",force=TRUE)


install.packages(c("ggplot2","readr"))

# Cargar librerías
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(readr)

# Establecer el directorio
setwd("C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045")
setwd("/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045") # Path Aurelio

# Leer el CSV
biclusters <- read.csv("3_Filtrados_elite/normal/FC_3/FC3_normal_DEGs_symbols_deduplicated_biclusters_elite.csv", stringsAsFactors = FALSE)

# Crear carpeta de resultados
#dir.create("resultados_enrichment", showWarnings = FALSE)

# Inicializar tabla resumen
resumen <- data.frame(
  Bicluster = character(),
  Min_GO_qvalue = numeric(),
  Top_GO_term = character(),
  stringsAsFactors = FALSE
)

# Recorrer cada bicluster
for (i in 1:nrow(biclusters)) {
  id <- biclusters$Bicluster[i]
  genes_raw <- biclusters$Genes[i]
  genes <- unlist(strsplit(genes_raw, ";"))
  genes <- trimws(genes)
  
  # Mapear a ENTREZ
  entrez_ids <- mapIds(hgu133plus2.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- na.omit(entrez_ids)
  
  # Enriquecimiento GO
  go <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "ALL", pAdjustMethod = "BH")
  if (!is.null(go) && nrow(go@result) > 0) {
    min_go <- min(go@result$p.adjust, na.rm = TRUE)
    top_go <- go@result$Description[which.min(go@result$p.adjust)]
    write.csv(go@result, paste0("4_resultados_enrichment/normal/FC_3/Bicluster_", id, "_GO.csv"), row.names = FALSE)
    
    plot_go <- dotplot(go, split = "ONTOLOGY")
    
    # Verificar que 'ONTOLOGY' existe en los resultados
    if ("ONTOLOGY" %in% colnames(go@result)) {
      plot_go <- plot_go + facet_grid(ONTOLOGY ~ ., scale = "free", space = "free")
    }
    
    ggsave(
      paste0("4_resultados_enrichment/normal/FC_3/Bicluster_", id, "_GO_plot.pdf"),
      plot_go
    )
  } else {
    min_go <- NA
    top_go <- NA
  }
  
  # Añadir a la tabla resumen
  resumen <- rbind(
    resumen,
    data.frame(Bicluster = id, Min_GO_qvalue = min_go, Top_GO_term = top_go)
  )
}

# Ordenar y guardar resumen
resumen$Min_GO_qvalue <- as.numeric(resumen$Min_GO_qvalue)
resumen <- resumen[order(resumen$Min_GO_qvalue, na.last = TRUE), ]
write.csv(resumen, "4_resultados_enrichment/normal/FC_3/resumen_pvalores_biclusters.csv", row.names = FALSE)

# Mostrar el más significativo
cat("\n📌 Bicluster más significativo (GO):\n")
print(head(resumen, 1))
