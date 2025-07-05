
# Analisis de enriquecimiento funcional
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DO.db")
BiocManager::install("enrichplot", force = TRUE)
install.packages("plyr")

# Establecer directorio de trabajo y cargar las librerias

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

#MARC GITHUB
setwd("C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/methodology/6geneEnrichment/input")
#MARC PORTATIL
#setwd("C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE17674/methodology/6geneEnrichment/input")

genelist1 <- readLines("Cluster_1.txt")
genelist2 <- readLines("Cluster_2.txt")
genelist7 <- readLines("Cluster_7.txt")
genelist20 <- readLines("Cluster_20.txt")

# PARA EL CLUSTER 1
entrez_ids1 <- mapIds(hgu133plus2.db, keys=genelist1, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids1,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../results/Cluster_1_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../results/Cluster_1_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids1, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../results/Cluster_1_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../results/Cluster_1_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 2
entrez_ids2 <- mapIds(hgu133plus2.db, keys=genelist2, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids2,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../results/Cluster_2_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../results/Cluster_2_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids2, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../results/Cluster_2_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../results/Cluster_2_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 7
entrez_ids7 <- mapIds(hgu133plus2.db, keys=genelist7, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids7,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../results/Cluster_7_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../results/Cluster_7_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids7, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../results/Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../results/Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 20
entrez_ids20 <- mapIds(hgu133plus2.db, keys=genelist20, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids20,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../results/Cluster_2_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../results/Cluster_20_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids20, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../results/Cluster_20_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../results/Cluster_20_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}