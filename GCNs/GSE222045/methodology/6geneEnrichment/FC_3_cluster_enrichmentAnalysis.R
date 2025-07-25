
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
setwd("C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/6geneEnrichment/FC_3/input/")
# MARC PORTATIL
#setwd("C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/6geneEnrichment")


# Unidades

genelist1 <- readLines("Cluster_1.txt")
genelist2 <- readLines("Cluster_2.txt")
genelist3 <- readLines("Cluster_3.txt")
genelist5 <- readLines("Cluster_5.txt")
genelist6 <- readLines("Cluster_6.txt")
genelist22 <- readLines("Cluster_22.txt")




# PARA EL CLUSTER 1
entrez_ids1 <- mapIds(hgu133plus2.db, keys=genelist1, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids1,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "Cluster_1_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "Cluster_1_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids1, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "Cluster_1_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "Cluster_1_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}



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





# PARA EL CLUSTER 3
entrez_ids3 <- mapIds(hgu133plus2.db, keys=genelist3, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids3,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_3_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_3_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids3, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_3_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_3_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
    #ggsave(filename = "importanteBreast_KEGG_BH_finla.svg", plot = grafico, device = "svg")
  }
}




# PARA EL CLUSTER 5
entrez_ids5 <- mapIds(hgu133plus2.db, keys=genelist5, column='ENTREZID', keytype='SYMBOL')

go <- enrichGO(entrez_ids5,  OrgDb = hgu133plus2.db, ont = "all", pAdjustMethod = "BH")
if(length(go@result$ONTOLOGY) != 0){
  write.csv(go, "../Cluster_5_all_BH.csv", row.names = FALSE)
  grafico <- dotplot(go, split="ONTOLOGY", font.size = 10) + 
    facet_grid(ONTOLOGY~., scale ='free', space = 'free')
  ggsave(filename = "../Cluster_5_all_BH_final.pdf", plot = grafico, device = "pdf", height = 9)
  #ggsave(filename = "importanteBreast_all_BH_final.svg", plot = grafico, device = "svg", height = 11)
  
}

kegg <- enrichKEGG(entrez_ids5, pAdjustMethod = "BH")
if(!is.null(kegg)){
  if(any(kegg@result$p.adjust <0.05)){
    write.csv(kegg, "../Cluster_5_normalSarcoma_MajorRevision_resolution_KEGG_BH.csv", row.names = FALSE)
    grafico <- dotplot(kegg) + facet_grid(scale="free")
    ggsave(filename = "../Cluster_5_normalSarcoma_MajorRevision_resolution_KEGG_BH_final.pdf", plot = grafico, device = "pdf")
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


