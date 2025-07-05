import pandas as pd
import os

# === Rutas ===
# Marc
#ruta_expr = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/1_DEGs/FC_1/tumor_DEGs_symbols.csv"
#ruta_bicluster = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/5_degree/tumor/FC_1/mapping_bicluster_tumor.csv"
#output_dir = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/5_degree/tumor/FC_1"

# Aurelio
ruta_expr = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/1_DEGs/FC_1/tumor_DEGs_symbols.csv"
ruta_bicluster = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/5_degree/tumor/FC_1/mapping_bicluster_tumor.csv"
output_dir = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/5_degree/tumor/FC_1"

# === Crear carpeta si no existe ===
os.makedirs(output_dir, exist_ok=True)

# === 1. Leer archivos ===
expr = pd.read_csv(ruta_expr, sep="\t")
expr.set_index("attr_name", inplace=True)

bicluster = pd.read_csv(ruta_bicluster)
genes_bicluster = bicluster.loc[0, "Genes"].split(";")
genes_bicluster = [g.strip() for g in genes_bicluster if g.strip() in expr.index]

# === 2. Submatriz de expresión ===
expr_bicluster = expr.loc[genes_bicluster]

# === 3. Matriz de correlación Pearson ===
corr_matrix = expr_bicluster.T.corr()

# === 4. Matriz binaria ===
binary_matrix = ((corr_matrix > 0.5) | (corr_matrix < -0.5)).astype(int)

# === 5. Guardar resultados ===
corr_matrix.to_csv(os.path.join(output_dir, "matriz_correlacion_mapping_bicluster_tumor.csv"))
binary_matrix.to_csv(os.path.join(output_dir, "matriz_binaria_mapping_bicluster_tumor.csv"))

print("✅ Matrices generadas y guardadas correctamente en:")
print(output_dir)


# === 6. Calcular grados ===
grados = binary_matrix.sum(axis=1)
grados_df = pd.DataFrame({"Gene": grados.index, "Degree": grados.values})

# === 7. Guardar todo ===
corr_matrix.to_csv(os.path.join(output_dir, "matriz_correlacion_mapping_bicluster_tumor.csv"))
binary_matrix.to_csv(os.path.join(output_dir, "matriz_binaria_mapping_bicluster_tumor.csv"))
grados_df.to_csv(os.path.join(output_dir, "grados_mapping_bicluster_tumor.csv"), index=False)

print("✅ Todo generado correctamente:")
print(f"- Matriz de correlación: {corr_matrix.shape}")
print(f"- Matriz binaria: {binary_matrix.shape}")
print(f"- Grados guardados: {grados_df.shape[0]} genes")
