
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
setwd("C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/6geneEnrichment")
# MARC PORTATIL
#setwd("C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/6geneEnrichment")


# Unidades

genelist2 <- readLines("input/Cluster_2.txt")
genelist5 <- readLines("input/Cluster_5.txt")
genelist7 <- readLines("input/Cluster_7.txt")
genelist8 <- readLines("input/Cluster_8.txt")
genelist9 <- readLines("input/Cluster_9.txt")
genelist10 <- readLines("input/Cluster_10.txt")
genelist11 <- readLines("input/Cluster_11.txt")
genelist15 <- readLines("input/Cluster_15.txt")
genelist26 <- readLines("input/Cluster_26.txt")
genelist32 <- readLines("input/Cluster_32.txt")
genelist34 <- readLines("input/Cluster_34.txt")



# PARA EL CLUSTER 2
entrez_ids2 <- mapIds(hgu133plus2.db, keys=genelist2, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids2,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_2_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_2_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids2, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_2_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_2_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



# PARA EL CLUSTER 5
entrez_ids5 <- mapIds(hgu133plus2.db, keys=genelist5, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids5,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_5_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_5_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids5, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_5_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_5_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}

# PARA EL CLUSTER 7
entrez_ids7 <- mapIds(hgu133plus2.db, keys=genelist7, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids7,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_7_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_7_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids7, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 8
entrez_ids8 <- mapIds(hgu133plus2.db, keys=genelist8, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids8,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_8_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_8_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids8, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_8_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_8_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 9
entrez_ids9 <- mapIds(hgu133plus2.db, keys=genelist9, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids9,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_9_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_9_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids9, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_9_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_9_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 10
entrez_ids10 <- mapIds(hgu133plus2.db, keys=genelist10, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids10,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_10_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_10_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids10, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_10_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_10_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 11
entrez_ids11 <- mapIds(hgu133plus2.db, keys=genelist11, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids11,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_11_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_11_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids11, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_11_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_11_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}




# PARA EL CLUSTER 15
entrez_ids15 <- mapIds(hgu133plus2.db, keys=genelist15, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids15,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_15_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_15_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids15, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_15_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_15_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 26
entrez_ids26 <- mapIds(hgu133plus2.db, keys=genelist26, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids26,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_26_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_26_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids26, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_26_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_26_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 32
entrez_ids32 <- mapIds(hgu133plus2.db, keys=genelist32, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids32,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_32_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_32_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids32, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_32_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_32_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 34
entrez_ids34 <- mapIds(hgu133plus2.db, keys=genelist34, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids34,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_34_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_34_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids34, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_34_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_34_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}
