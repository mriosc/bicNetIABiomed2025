import pandas as pd
import os

"""
A partir de los biclusters de mapeo:
    
1) Calcular el Pearson correlación (Baak et al., 2020) entre los genes presentes en el bicluster importante; es decir, 
si el bicluster consta de 50 genes con 6 muestras entonces creará un 50 × 50 matriz de correlación dependiendo de los valores de expresión de 6 muestras.
    
2) Convertir el matriz de correlación en una matriz binaria considerando los valores entre −0,5 to 0,5 como 0 ′ s y los restos como 1 ′ s.
    
3) Calcule el grado de genes en el bicluster de mapeo de enfermedades recuperado en el paso 10 para encontrar genes sospechosos que 
tengan el mismo grado que los genes de élite y guárdelos para un examen más detallado. Debido a la variabilidad de los valores de expresión, 
los genes que tienen el mismo grado que un gen de élite en una condición normal pueden comportarse de manera diferente en el estado de 
enfermedad.
    
"""

# === Rutas ===
# Marc
#ruta_expr = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/1_DEGs/FC_1/FC1_tumor_DEGs_symbols.csv"
#ruta_bicluster = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/5_degree/tumor/FC_1/mapping_bicluster_tumor.csv"
#output_dir = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/5_degree/tumor/FC_1"

# Aurelio
ruta_expr = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE64763/1_DEGs/FC_1/FC1_tumor_DEGs_symbols.csv"
ruta_bicluster = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE64763/5_degree/tumor/FC_1/mapping_bicluster_tumor.csv"
output_dir = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE64763/5_degree/tumor/FC_1"

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

print("Matrices generadas y guardadas correctamente en:")
print(output_dir)


# === 6. Calcular grados ===
grados = binary_matrix.sum(axis=1)
grados_df = pd.DataFrame({"Gene": grados.index, "Degree": grados.values})

# === 7. Guardar todo ===
corr_matrix.to_csv(os.path.join(output_dir, "matriz_correlacion_mapping_bicluster_tumor.csv"))
binary_matrix.to_csv(os.path.join(output_dir, "matriz_binaria_mapping_bicluster_tumor.csv"))
grados_df.to_csv(os.path.join(output_dir, "grados_mapping_bicluster_tumor.csv"), index=False)

print("Todo generado correctamente:")
print(f"- Matriz de correlación: {corr_matrix.shape}")
print(f"- Matriz binaria: {binary_matrix.shape}")
print(f"- Grados guardados: {grados_df.shape[0]} genes")
