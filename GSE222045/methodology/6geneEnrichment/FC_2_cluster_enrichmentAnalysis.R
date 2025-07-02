
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
setwd("C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/6geneEnrichment/FC_2/input")
# MARC PORTATIL
#setwd("C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/6geneEnrichment")


# Unidades

genelist2 <- readLines("Cluster_2.txt")
genelist4 <- readLines("Cluster_4.txt")
genelist6 <- readLines("Cluster_6.txt")
genelist7 <- readLines("Cluster_7.txt")
genelist9 <- readLines("Cluster_9.txt")
genelist12 <- readLines("Cluster_12.txt")
genelist13 <- readLines("Cluster_13.txt")
genelist15 <- readLines("Cluster_15.txt")
genelist16 <- readLines("Cluster_16.txt")
genelist22 <- readLines("Cluster_22.txt")
genelist36 <- readLines("Cluster_36.txt")



# PARA EL CLUSTER 2
entrez_ids2 <- mapIds(hgu133plus2.db, keys=genelist2, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids2,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_2_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_2_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids2, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_2_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_2_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



# PARA EL CLUSTER 4
entrez_ids4 <- mapIds(hgu133plus2.db, keys=genelist4, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids4,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_4_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_4_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids4, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_4_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_4_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



# PARA EL CLUSTER 6
entrez_ids6 <- mapIds(hgu133plus2.db, keys=genelist6, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids6,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_6_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_6_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids6, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_6_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_6_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



# PARA EL CLUSTER 7
entrez_ids7 <- mapIds(hgu133plus2.db, keys=genelist7, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids7,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_7_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_7_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids7, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_7_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



# PARA EL CLUSTER 9
entrez_ids9 <- mapIds(hgu133plus2.db, keys=genelist9, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids9,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_9_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_9_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids9, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_9_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_9_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 12
entrez_ids12 <- mapIds(hgu133plus2.db, keys=genelist12, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids12,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_12_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_12_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids12, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_12_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_12_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 13
entrez_ids13 <- mapIds(hgu133plus2.db, keys=genelist13, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids13,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_13_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_13_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids13, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_13_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_13_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



# PARA EL CLUSTER 15
entrez_ids15 <- mapIds(hgu133plus2.db, keys=genelist15, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids15,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_15_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_15_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids15, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_15_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_15_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 16
entrez_ids16 <- mapIds(hgu133plus2.db, keys=genelist16, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids16,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_16_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_16_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids16, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_16_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_16_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 22
entrez_ids22 <- mapIds(hgu133plus2.db, keys=genelist22, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids22,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_22_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_22_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids22, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_22_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_22_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}


# PARA EL CLUSTER 36
entrez_ids36 <- mapIds(hgu133plus2.db, keys=genelist36, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids36,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_36_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_36_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids36, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_36_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_36_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}
