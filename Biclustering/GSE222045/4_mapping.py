import pandas as pd
import os

"""
Extraer biclusters de mapeo del conjunto de bicluster tumorales que corresponde al más bajo p-valor bicluster obtenido del conjunto de bicluster
normal. Los biclusters que tienen los mismos genes de élite, que están presentes tanto en condiciones normales como patológicas, se denominan bicluster de mapeo.
"""

#ruta_entrada = "C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/" # Ruta Marc
#ruta_salida = "C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/5_degree/tumor/FC_1" # Ruta Marc

ruta_entrada = "/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/" # Ruta Aurelio
ruta_salida = "/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/5_degree/tumor/FC_1" # Ruta Aurelio

# Archivo de genes élite del bicluster normal
archivo_genes_elite = ruta_entrada+"5_degree/normal/FC_1/genes_mismo_grado_elite.csv"

# Archivo de biclusters tumorales
archivo_biclusters_tumor = ruta_entrada+"2_biclusters/FC_1/tumor_DEGs_symbols_biclusters.csv"

# Ruta de salida
os.makedirs(ruta_salida, exist_ok=True)  # Crea los directorios si no existen

# Leer lista de genes élite
elite_df = pd.read_csv(archivo_genes_elite)
elite_genes = set(elite_df["Gene"].dropna().str.strip().str.upper())

# Leer biclusters tumorales
df = pd.read_csv(archivo_biclusters_tumor)
df["Genes"] = df["Genes"].astype(str)

# Buscar el bicluster con mayor número de genes élite en común
mejor_bicluster = None
max_matches = 0
genes_encontrados = set()

for idx, row in df.iterrows():
    genes = set(g.strip().upper() for g in row["Genes"].split(";"))
    interseccion = genes.intersection(elite_genes)
    if len(interseccion) > max_matches:
        max_matches = len(interseccion)
        mejor_bicluster = row
        genes_encontrados = interseccion

# Guardar el mapping bicluster
if mejor_bicluster is not None:
    print(f"Mapping bicluster encontrado con {max_matches} genes élite en común.")
    print("Genes en común:", ", ".join(genes_encontrados))

    output = pd.DataFrame([mejor_bicluster])
    output_path = os.path.join(ruta_salida, "mapping_bicluster_tumor.csv")
    output.to_csv(output_path, index=False)
    print(f"Guardado como: {output_path}")
else:
    print("No se encontró ningún bicluster con genes élite en común.")
