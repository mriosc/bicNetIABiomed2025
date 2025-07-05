import os
import pandas as pd
import numpy as np

"""
A partir del mejor bicluster normal segun su p-value:
    
1) Calcular el Pearson correlación (Baak et al., 2020) entre los genes presentes en el bicluster importante; es decir, 
si el bicluster consta de 50 genes con 6 muestras entonces creará un 50 × 50 matriz de correlación dependiendo de los valores de expresión de 6 muestras.
    
2) Convertir el matriz de correlación en una matriz binaria considerando los valores entre −0,5 to 0,5 como 0 ′ s y los restos como 1 ′ s.
    
3) Calcule el grado de cada gen para identificar los genes sospechosos que tienen el mismo grado que los genes de élite y guárdelos para una mayor investigación 
porque es probable que los genes con el mismo grado que los genes de élite se conviertan en un biomarcador potencial.
    
"""


# === 0. Establecer directorio de trabajo ===
#os.chdir(r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045")
os.chdir("/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE64763") # Ruta Aurelio


# === 1. Cargar archivos ===
resumen = pd.read_csv("4_resultados_enrichment/normal/FC_1/resumen_pvalores_biclusters.csv")
biclusters = pd.read_csv("3_Filtrados_elite/normal/FC_1/FC1_normal_DEGs_symbols_biclusters_elite.csv")

with open("3_Filtrados_elite/genes_elite.txt", "r") as f:
    elite_genes = [g.strip() for g in f.read().replace("\n", "").split(",")]

# Leer la matriz de expresión separada por tabuladores
expr = pd.read_csv("1_DEGs/FC_1/FC1_normal_DEGs_symbols.csv", sep="\t")
expr = expr.set_index(expr.columns[0])

# === 2. Seleccionar el bicluster con mejor p-valor ===
best_bicluster_id = resumen.loc[resumen["Min_GO_qvalue"].idxmin(), "Bicluster"]
genes_raw = biclusters[biclusters["Bicluster"] == best_bicluster_id]["Genes"].values[0]
genes = [g.strip() for g in genes_raw.split(";")]

# === 3. Filtrar expresión para genes del bicluster ===
expr_bicluster = expr.loc[expr.index.intersection(genes)]

# === 4. Calcular matriz de correlación Pearson (genes x genes) ===
corr_matrix = expr_bicluster.T.corr()

# === 5. Convertir a matriz binaria (1 si corr > 0.5 o < -0.5, 0 en otro caso) ===
binary_matrix = ((corr_matrix < -0.5) | (corr_matrix > 0.5)).astype(int)
binary_matrix.to_csv("5_degree/normal/FC_1/matriz_binaria_bicluster.csv")

# === 6. Calcular grados de los genes ===
degrees = binary_matrix.sum(axis=1)

# === 7. Guardar grados de todos los genes de los biclusters ===
degrees_df = pd.DataFrame({"Gene": degrees.index, "Degree": degrees.values})
degrees_df["Elite"] = degrees_df["Gene"].isin(elite_genes)
degrees_df.to_csv("5_degree/normal/FC_1/grados_todos_genes_bicluster.csv", index=False)

# === 8. Identificar genes élite y comparar grados con tolerancia ===
tolerancia = 0
elite_in_bicluster = set(genes).intersection(elite_genes)
elite_degrees = degrees.loc[list(elite_in_bicluster)]

matched_genes_list = []
for elite_gene, elite_deg in elite_degrees.items():
    match = degrees[(degrees >= elite_deg - tolerancia) & (degrees <= elite_deg + tolerancia)]
    for gene, deg in match.items():
        matched_genes_list.append({"Gene": gene, "Degree": deg, "Elite_Degree": elite_deg, "Elite_Gene": elite_gene})

matched_genes_df = pd.DataFrame(matched_genes_list)
matched_genes_df.to_csv("5_degree/normal/FC_1/genes_mismo_grado_elite.csv", index=False)

print("\n✅ Script completado.")
print(f"Bicluster más significativo: {best_bicluster_id}")
print(f"Genes con el mismo grado que genes élite: {matched_genes_df.shape[0]} → guardados en 'genes_mismo_grado_elite.csv'")
