## Procesamiento datos Microarray
install.packages("gplots")
install.packages("RColorBrewer")

# instalacion de Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install() # gestiona paquetes de bioconductor

# instalacion de paquetes a partir de BiocManager
BiocManager::install("preprocessCore")
BiocManager::install("affy", force = TRUE)
BiocManager::install("limma", force = TRUE)
BiocManager::install("annaffy", force = TRUE)
BiocManager::install("affyPLM", force = TRUE)
BiocManager::install("arrayQualityMetrics", force = TRUE)
BiocManager::install("hgu133plus2.db", force = TRUE)

# Cambiar ruta directorio de trabajo
#setwd("C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE17674/methodology/1getExpressionData") # Marc prueba
setwd("C:/Users/marcr/OneDrive/Escritorio/IABioMed/GSE64763/methodology/1getExpressionData") # Marc
#setwd("/Users/aurelio/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/1getExpressionData") # Aurelio MacOSX 
#setwd("/home/principalpc/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/1getExpressionData") # Aurelio Ubuntu 

# iniciacion de los paquetes
library("BiocManager")
library("affy")
library("limma")
library("dplyr")
library("readr")
library("tibble")
library("annaffy")
library("affyPLM")
library("arrayQualityMetrics")
library("hgu133plus2.db")
library("gplots")


# 1) Construccion del microarray a partir de datos en formato CEL
# ------------------------------------------------------------
#SDRF <- read.delim("C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE17674/dataset/data_info.csv", sep = ",") # Marc Prueba
SDRF <- read.delim("C:/Users/marcr/OneDrive/Escritorio/IABioMed/GSE64763/dataset/datainfo.csv", sep = ",") # Marc
#SDRF <- read.delim("/Users/aurelio/git-repositories/sarcoma-iwbbio2025/GSE17674/dataset/data_info.csv", sep = ",") # Aurelio MacOSX
#SDRF <- read.delim("/home/principalpc/git-repositories/sarcoma-iwbbio2025/GSE17674/dataset/data_info.csv", sep = ",") # Aurelio Ubuntu

rownames(SDRF) <- SDRF$Ids # Pone el nombre de los samples como identificador de las filas del archivo dataInfo.
SDRF <- AnnotatedDataFrame(SDRF) # Convierte la matrix en un objeto AnnotatedDataFrame para que sea compatible con la funcion ReadAffy()
microarray.raw.data <- ReadAffy(filenames = SDRF$Array.Data.File, verbose=TRUE, phenoData = SDRF) 
microarray.raw.data


### 2) QUALITY CONTROL (SOLO PARA MICROARRAYS)
# -----------------------------------------
# imagenes de escaneres (daños fisicos)  
dev.off() # Reiniciar el sistema plot de RStudio
image(microarray.raw.data[,1], col=rainbow(100))
image(microarray.raw.data[,2], col=rainbow(100))
image(microarray.raw.data[,3], col=rainbow(100))
image(microarray.raw.data[,4], col=rainbow(100))
image(microarray.raw.data[,5], col=rainbow(100))
image(microarray.raw.data[,6], col=rainbow(100))
image(microarray.raw.data[,7], col=rainbow(100))
image(microarray.raw.data[,8], col=rainbow(100))
image(microarray.raw.data[,9], col=rainbow(100))
image(microarray.raw.data[,10], col=rainbow(100))
image(microarray.raw.data[,11], col=rainbow(100))
image(microarray.raw.data[,12], col=rainbow(100))

# metricas de calidad
arrayQualityMetrics(expressionset = microarray.raw.data, outdir = "results/qualityMetrics", force = TRUE, do.logtransform = TRUE, intgroup = "type")

### 3) DATA PREPROCESSING
# ---------------------------------------------
# Comprobar que los valores del microarray no están normalizados
boxplot(microarray.raw.data, col = rainbow(12), las=2, ylab="Fluorescence")

# NORMALIZAR: Robust multiarray average (rma) realiza la correccion de la fluorescencia de fondo, 
# normaliza y calcula los niveles de expresion con transformacion de los datos en log2
microarray.processed.data <-rma(microarray.raw.data)

# Comprobar de nuevo que los valores del microarray ya se han normalizado
boxplot(microarray.processed.data, col = rainbow(12), las=2, ylab="Fluorescence A.U.")

# Guardar matriz de expresion normalizada
write.table(exprs(microarray.processed.data), file="results/rawAffy_expressionLevel_rmaData.csv", sep = ",", 
            quote = FALSE, row.names = TRUE, col.names = TRUE)

# Anadir la palabra "affy_ids" al CSV
modArchivo <- readLines("results/rawAffy_expressionLevel_rmaData.csv")
modArchivo[1] <- paste("affy_ids", modArchivo[1], sep = ",")
writeLines(modArchivo, "results/rawAffy_expressionLevel_rmaData.csv")

# Creacion de la matriz de expresion normalizada
expression.level <- exprs(microarray.processed.data)
head(expression.level)
dim(expression.level)

# Guardar toda la info por cada uno de los genes de la matriz de expresion
totalGenes <- aafTableAnn(rownames(expression.level), "hgu133plus2.db", aaf.handler())
saveText(totalGenes, file="results/aafTableAnn_majorRevision.txt")

# modificar nombre de las columnas
sampleID <- c("Fibroid_1","Fibroid_2","Fibroid_3","Fibroid_4",
              "Fibroid_5","Fibroid_6","Fibroid_7","Fibroid_8",
              "Fibroid_9","Fibroid_10","Fibroid_11","Fibroid_12",
              "Fibroid_13", "Fibroid_14", "Fibroid_15", "Fibroid_16",
              "Fibroid_17", "Fibroid_18", "Fibroid_19", "Fibroid_20",
              "Fibroid_21", "Fibroid_22", "Fibroid_23", "Fibroid_24",
              "Fibroid_25", "Tumor_1", "Tumor_2", "Tumor_3",
              "Tumor_4", "Tumor_5", "Tumor_6", "Tumor_7",
              "Tumor_8", "Tumor_9", "Tumor_10", "Tumor_11",
              "Tumor_12", "Tumor_13", "Tumor_14", "Tumor_15",
              "Tumor_16", "Tumor_17", "Tumor_18", "Tumor_19",
              "Tumor_20", "Tumor_21", "Tumor_22", 
              "Tumor_23", "Tumor_24", "Tumor_25",
              "NormalMyometrium_1", "NormalMyometrium_2", "NormalMyometrium_3", 
              "NormalMyometrium_4", "NormalMyometrium_5", "NormalMyometrium_6",
              "NormalMyometrium_7", "NormalMyometrium_8", "NormalMyometrium_9",
              "NormalMyometrium_10", "NormalMyometrium_11", "NormalMyometrium_12", 
              "NormalMyometrium_13", "NormalMyometrium_14", "NormalMyometrium_15", 
              "NormalMyometrium_16", "NormalMyometrium_17", "NormalMyometrium_18",
              "NormalMyometrium_19", "NormalMyometrium_20", "NormalMyometrium_21",
              "NormalMyometrium_22", "NormalMyometrium_23", "NormalMyometrium_24",
              "NormalMyometrium_25", "NormalMyometrium_26", "NormalMyometrium_27",
              "NormalMyometrium_28", "NormalMyometrium_29")

colnames(expression.level) <- sampleID
head(expression.level)

### 4) ANALISIS DEGs (USANDO Limma)
# ---------------------------------------------
# Separar tumor y control
controlExpression <- expression.level[, c(1:24, 50:79)]  # Fibroid (1:24) + NormalMyometrium (50:79)
tumorExpression <- expression.level[, 25:49]             # Tumor (25:49)

# Media de cada condicion
tumorMean = apply(tumorExpression, 1, mean)
controlMean = apply(controlExpression, 1, mean)
head(tumorMean)
head(controlMean)

# Maximum of all the means
limit = max(tumorMean, controlMean)

# Scatter plot
plot(controlMean ~ tumorMean, xlab = "TUMOR", ylab = "CONTROL",
     main = "Dataset", xlim = c(0, limit), ylim = c(0, limit))
abline(0, 1, col = "red")

# Capturar fold-change (significancia biologica - diferencia de las medias entre condiciones)
fold = tumorMean - controlMean

# Histograma de diferencias de fold-change
hist(fold, col = "gray")

# Paquete Limma
experimental.design <- model.matrix(~ -1 + factor(c(
  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,   # Tumor_1 to Tumor_25
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,   # Fibroid_1 to Fibroid_25
  1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1  # NormalMyometrium_1 to _29
)))
 # 1 -> 43 columnas tumores // 2 -> 18 columnas normales
colnames(experimental.design) <- c("tumor", "normal")
head(experimental.design)
linear.fit <- lmFit(expression.level, experimental.design)
contrast.matrix <-makeContrasts(tumor - normal, levels = c("tumor","normal"))  # especificacion de contrastes
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix) # calculo del Fold-Change
contrast.results <- eBayes(contrast.linear.fit) # calculo p-valor/q-valor
results <- topTable(contrast.results, number = 22277, coef = 1, sort.by = "logFC") # Coeficiente 1 ya que el interés es el tumoral. 54675 se refiere al numero de genes.
head(results)
fold.change.results <- results$logFC
adj.Pval.results <- results$adj.P.Val
genes.ids.results <- rownames(results)

# metodo combinado para los genes activados y reprimidos
activatedGenes <- genes.ids.results[fold.change.results > 1 & adj.Pval.results < 0.05]
repressedGenes <- genes.ids.results[fold.change.results < -1 & adj.Pval.results < 0.05]
length(activatedGenes)
length(repressedGenes)

# visualizacion de los datos volcano plot
names(fold.change.results) <- genes.ids.results
log.padj.results <- -log10(adj.Pval.results)
names(log.padj.results) <- genes.ids.results
plot(fold.change.results,log.padj.results, ylab="-log10(p value)", xlab= "log2 fold change", pch=19, cex=0.5,
     col = "grey", xlim=c(-6,6))
points(fold.change.results[activatedGenes],log.padj.results[activatedGenes], 
       pch = 19, cex = 0.5, col = "red")
points(fold.change.results[repressedGenes],log.padj.results[repressedGenes], 
       pch = 19, cex = 0.5, col = "blue")


# DEGs
allDEGs <- genes.ids.results[abs(fold.change.results) > 1 & adj.Pval.results < 0.05]
head(allDEGs)

# Normal DEGs
normalDEGs <- expression.level[, c("Fibroid_1","Fibroid_2","Fibroid_3","Fibroid_4",
              "Fibroid_5","Fibroid_6","Fibroid_7","Fibroid_8",
              "Fibroid_9","Fibroid_10","Fibroid_11","Fibroid_12",
              "Fibroid_13", "Fibroid_14", "Fibroid_15", "Fibroid_16",
              "Fibroid_17", "Fibroid_18", "Fibroid_19", "Fibroid_20",
              "Fibroid_21", "Fibroid_22", "Fibroid_23", "Fibroid_24",
              "Fibroid_25", "NormalMyometrium_1", "NormalMyometrium_2", 
              "NormalMyometrium_3", "NormalMyometrium_4", "NormalMyometrium_5", "NormalMyometrium_6",
              "NormalMyometrium_7", "NormalMyometrium_8", "NormalMyometrium_9",
              "NormalMyometrium_10", "NormalMyometrium_11", "NormalMyometrium_12", 
              "NormalMyometrium_13", "NormalMyometrium_14", "NormalMyometrium_15", 
              "NormalMyometrium_16", "NormalMyometrium_17", "NormalMyometrium_18",
              "NormalMyometrium_19", "NormalMyometrium_20", "NormalMyometrium_21",
              "NormalMyometrium_22", "NormalMyometrium_23", "NormalMyometrium_24",
              "NormalMyometrium_25", "NormalMyometrium_26", "NormalMyometrium_27",
              "NormalMyometrium_28", "NormalMyometrium_29")]

normalDEGs <- normalDEGs[rownames(normalDEGs) %in% allDEGs,]
normalDEGs <- cbind(attr_name = rownames(normalDEGs), normalDEGs)

write.table(normalDEGs, file="results/3normal_DEGs.csv", sep="\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

tumorDEGs <- expression.level[, c("Tumor_1","Tumor_2","Tumor_3","Tumor_4",
                                                "Tumor_5","Tumor_6","Tumor_7","Tumor_8",
                                                "Tumor_9","Tumor_10","Tumor_11","Tumor_12",
                                                "Tumor_13", "Tumor_14", "Tumor_15", "Tumor_16",
                                                "Tumor_17", "Tumor_18", "Tumor_19", "Tumor_20",
                                                "Tumor_21", "Tumor_22", "Tumor_23", "Tumor_24",
                                                "Tumor_25")]
tumorDEGs <- tumorDEGs[rownames(tumorDEGs) %in% allDEGs,]
tumorDEGs <- cbind(attr_name = rownames(tumorDEGs), tumorDEGs)

write.table(tumorDEGs, file="results/3tumor_DEGs.csv", sep="\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


################################################################################

# OBTENCIÓN DE DEGs SOLO CON TUMOR Y CONTROLES SANOS

################################################################################


# 1. Selección de muestras (solo Tumor vs Normal)
tumor_samples <- 26:50
normal_samples <- 51:79

# 2. Subconjunto de matriz de expresión
subset_expr <- expression.level[, c(tumor_samples, normal_samples)]

# 3. Diseño experimental para Tumor vs Normal
grupo <- factor(c(rep("tumor", length(tumor_samples)), 
                  rep("normal", length(normal_samples))))
design <- model.matrix(~0 + grupo)
colnames(design) <- levels(grupo)

# 4. Análisis con limma
fit <- lmFit(subset_expr, design)
contrast.matrix <- makeContrasts(tumor - normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 5. Extracción de resultados ordenados por logFC
results_sub <- topTable(fit2, coef = 1, number = Inf, sort.by = "logFC")

# 6. Obtener DEG usando criterios estándar
fold.change.sub <- results_sub$logFC
adj.pval.sub <- results_sub$adj.P.Val
genes.sub <- rownames(results_sub)

activatedGenes_sub <- genes.sub[fold.change.sub > 2 & adj.pval.sub < 0.05]
repressedGenes_sub <- genes.sub[fold.change.sub < -2 & adj.pval.sub < 0.05]

cat("Genes activados:", length(activatedGenes_sub), "\n")
cat("Genes reprimidos:", length(repressedGenes_sub), "\n")

# 7. Volcano plot
log.padj.sub <- -log10(adj.pval.sub)
plot(fold.change.sub, log.padj.sub, pch = 19, cex = 0.5,
     col = "grey", xlim = c(-6, 6), 
     xlab = "log2 Fold Change", ylab = "-log10 adj.Pval",
     main = "Tumor vs Normal Myometrium")
points(fold.change.sub[activatedGenes_sub], log.padj.sub[activatedGenes_sub], 
       col = "red", pch = 19, cex = 0.5)
points(fold.change.sub[repressedGenes_sub], log.padj.sub[repressedGenes_sub], 
       col = "blue", pch = 19, cex = 0.5)
abline(h = -log10(0.05), col = "blue", lty = 2)
abline(v = c(-1, 1), col = "red", lty = 2)

# 1. Obtener DEGs con los nuevos resultados
allDEGs_sub <- genes.sub[abs(fold.change.sub) > 2 & adj.pval.sub < 0.05]

# 2. Extraer expresión de muestras normales (solo Normal Myometrium)
normal_samples_names <- colnames(expression.level)[51:79]
normalDEGs <- expression.level[, normal_samples_names]
normalDEGs <- normalDEGs[rownames(normalDEGs) %in% allDEGs_sub, ]
normalDEGs <- cbind(attr_name = rownames(normalDEGs), normalDEGs)

write.table(normalDEGs, file = "results/DEGs_NormalMyometrium.csv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# 3. Extraer expresión de muestras tumorales
tumor_samples_names <- colnames(expression.level)[26:50]
tumorDEGs <- expression.level[, tumor_samples_names]
tumorDEGs <- tumorDEGs[rownames(tumorDEGs) %in% allDEGs_sub, ]
tumorDEGs <- cbind(attr_name = rownames(tumorDEGs), tumorDEGs)

write.table(tumorDEGs, file = "results/DEGs_Tumor.csv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)












#### NO EJECUTAR A PARTIR DE AQUI: ESTA ES OTRA OPCION PARA CALCULAR LOS DEGS
### 4) ANALISIS DEGs (USANDO T-TEST)
# ---------------------------------------------
# Separar tumor y control
tumorExpression <- expression.level[,1:44]
controlExpression <- expression.level[,45:ncol(expression.level)]

# Media de cada condicion
tumorMean = apply(tumorExpression, 1, mean)
controlMean = apply(controlExpression, 1, mean)
head(tumorMean)
head(controlMean)

# Maximum of all the means
limit = max(tumorMean, controlMean)

# Scatter plot
plot(controlMean ~ tumorMean, xlab = "TUMOR", ylab = "CONTROL",
     main = "Dataset", xlim = c(0, limit), ylim = c(0, limit))
abline(0, 1, col = "red")

# Capturar fold-change (significancia biologica - diferencia de las medias entre condiciones)
fold = tumorMean - controlMean

# Histograma de diferencias de fold-change
hist(fold, col = "gray")

# Ejecutar significancia estadistica (usando t-test)
pvalue = NULL
tstat = NULL

for(i in 1 : nrow(tumorExpression)) { # Por cada gen : 
  x = tumorExpression[i,] # WT of gene number i
  y = controlExpression[i,] # KO of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

# Histograma de p-values (-log10)
hist(-log10(pvalue), col = "gray")

# Volcano: poner la significancia biologica (fold-change) y significancia estadistica (p-value) en un plot
plot(fold, -log10(pvalue), main = "Dataset - Volcano")

fold_cutoff = 3
pvalue_cutoff = 0.01
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(expression.level[filter_by_fold, ])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(expression.level[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue
filtered = expression.level[filter_combined,]
dim(filtered)
head(filtered)

# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
plot(fold, -log10(pvalue), main = "Dataset - Volcano #2")
points (fold[filter_combined], -log10(pvalue[filter_combined]),
        pch = 16, col = "red")

# Highlighting up-regulated in red and down-regulated in blue
plot(fold, -log10(pvalue), main = "Dataset - Volcano #3")
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "blue")

# Cluster the rows (genes) & columns (samples) by correlation
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))

# Generate a heatmap
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7)

# Enhanced heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")

# Save the heatmap to a PDF file
pdf ("results/datasetHeatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
dev.off()

# Save the DE genes to a text file
write.table (filtered, "dataset.txt", sep = "\t",
             quote = FALSE)