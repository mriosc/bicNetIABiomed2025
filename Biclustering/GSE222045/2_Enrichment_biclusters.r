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

# Leer el CSV
biclusters <- read.csv("3_Filtrados_elite/FC_1/F1normal_DEGs_symbols_biclusters_elite.csv", stringsAsFactors = FALSE)

# Crear carpeta de resultados
#dir.create("resultados_enrichment", showWarnings = FALSE)

# Inicializar tabla resumen
resumen <- data.frame(
  Bicluster = character(),
  Min_GO_qvalue = numeric(),
  Top_GO_term = character(),
  Min_KEGG_qvalue = numeric(),
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
    write.csv(go@result, paste0("4_resultados_enrichment/FC_1/Bicluster_", id, "_GO.csv"), row.names = FALSE)
    ggsave(
      paste0("4_resultados_enrichment/FC_1/Bicluster_", id, "_GO_plot.pdf"),
      dotplot(go, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scale = "free", space = "free")
    )
  } else {
    min_go <- NA
    top_go <- NA
  }

  # Enriquecimiento KEGG
  kegg <- enrichKEGG(entrez_ids, pAdjustMethod = "BH")
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    min_kegg <- min(kegg@result$p.adjust, na.rm = TRUE)
    write.csv(kegg@result, paste0("4_resultados_enrichment/FC_1/Bicluster_", id, "_KEGG.csv"), row.names = FALSE)
    ggsave(
      paste0("4_resultados_enrichment/FC_1/Bicluster_", id, "_KEGG_plot.pdf"),
      dotplot(kegg)
    )
  } else {
    min_kegg <- NA
  }

  # Añadir a la tabla resumen
  resumen <- rbind(
    resumen,
    data.frame(Bicluster = id, Min_GO_qvalue = min_go, Top_GO_term = top_go, Min_KEGG_qvalue = min_kegg)
  )
}

# Ordenar y guardar resumen
resumen$Min_GO_qvalue <- as.numeric(resumen$Min_GO_qvalue)
resumen <- resumen[order(resumen$Min_GO_qvalue, na.last = TRUE), ]
write.csv(resumen, "4_resultados_enrichment/FC_1/resumen_pvalores_biclusters.csv", row.names = FALSE)

# Mostrar el más significativo
cat("\n📌 Bicluster más significativo (GO):\n")
print(head(resumen, 1))
