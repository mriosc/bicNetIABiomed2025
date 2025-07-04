import pandas as pd

# === Parámetros ===
ruta_expr = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/1_DEGs/FC_1/FC1_tumor_DEGs_symbols_deduplicated.csv"
ruta_bicluster = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/5_degree/tumor/FC_1/mapping_bicluster_tumor.csv"
output_dir = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/5_degree/tumor/FC_1/"

# === 1. Cargar datos ===
expr = pd.read_csv(ruta_expr)
expr = expr.set_index(expr.columns[0])

bicluster = pd.read_csv(ruta_bicluster)
genes = bicluster.loc[0, "Genes"].split(";")
genes = [g.strip() for g in genes if g.strip() in expr.index]

# === 2. Submatriz de expresión del bicluster ===
expr_bicluster = expr.loc[genes]

# === 3. Matriz de correlación de Pearson ===
corr_matrix = expr_bicluster.T.corr()
corr_matrix.to_csv(output_dir + "matriz_correlacion_mapping_bicluster_tumor.csv")

# === 4. Matriz binaria ===
binary_matrix = ((corr_matrix > 0.5) | (corr_matrix < -0.5)).astype(int)
binary_matrix.to_csv(output_dir + "matriz_binaria_mapping_bicluster_tumor.csv")



print("✅ Matrices guardadas: correlación y binaria.")
