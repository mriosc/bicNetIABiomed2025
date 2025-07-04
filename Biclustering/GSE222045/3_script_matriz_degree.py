import os
import pandas as pd
import numpy as np

# === 0. Establecer directorio de trabajo ===
os.chdir(r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045")

# === 1. Cargar archivos ===
resumen = pd.read_csv("4_resultados_enrichment/FC_1/resumen_pvalores_biclusters.csv")
biclusters = pd.read_csv("3_Filtrados_elite/FC_1/F1normal_DEGs_symbols_biclusters_elite.csv")

with open("3_Filtrados_elite/genes_elite.txt", "r") as f:
    elite_genes = [g.strip() for g in f.read().replace("\n", "").split(",")]

expr = pd.read_csv("1_DEGs/FC_1/normal_DEGs_symbols.csv", sep="\t")
expr = expr.set_index(expr.columns[0])

# === 2. Seleccionar el bicluster con mejor p-valor ===
best_bicluster_id = resumen.loc[resumen["Min_GO_qvalue"].idxmin(), "Bicluster"]
genes_raw = biclusters[biclusters["Bicluster"] == best_bicluster_id]["Genes"].values[0]
genes = [g.strip() for g in genes_raw.split(";")]

# === 3. Filtrar expresión para genes del bicluster ===
expr_bicluster = expr.loc[expr.index.intersection(genes)]

# === 4. Correlación y matriz binaria ===
corr_matrix = expr_bicluster.T.corr()
binary_matrix = ((corr_matrix < -0.5) | (corr_matrix > 0.5)).astype(int)
binary_matrix.to_csv("matriz_binaria_bicluster.csv")

# === 5. Calcular grados ===
degrees = binary_matrix.sum(axis=1)

# === 6. Asociar todos los genes con los genes élite con los que comparten grado ===
elite_in_bicluster = set(genes).intersection(elite_genes)
elite_degrees = degrees.loc[list(elite_in_bicluster)]

# Creamos una lista para almacenar los resultados
resultados = []

for gene, deg in degrees.items():
    matched_elites = elite_degrees[elite_degrees == deg].index.tolist()
    if matched_elites:
        resultados.append({
            "Gene": gene,
            "Degree": deg,
            "MatchedEliteGenes": ",".join(matched_elites),
            "IsElite": gene in elite_genes
        })

# === 7. Guardar resultados ===
df_resultados = pd.DataFrame(resultados)
df_resultados.to_csv("genes_mismo_grado_elite.csv", index=False)

# También exportamos todos los grados
degrees_df = pd.DataFrame({
    "Gene": degrees.index,
    "Degree": degrees.values,
    "Elite": degrees.index.isin(elite_genes)
})
degrees_df.to_csv("grados_todos_genes_bicluster.csv", index=False)

print("\n✅ Script corregido completado.")
print(f"Bicluster más significativo: {best_bicluster_id}")
print(f"Genes con grado igual a élite: {len(df_resultados)} → guardados en 'genes_mismo_grado_elite.csv'")