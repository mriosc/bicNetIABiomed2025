import pandas as pd
import os

# Ruta ABSOLUTA al archivo de genes élite del bicluster normal
archivo_genes_elite = "C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/5_degree/normal/FC_1/genes_mismo_grado_elite.csv"

# Ruta ABSOLUTA al archivo de biclusters tumorales
archivo_biclusters_tumor = "C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/2_biclusters/FC_1/tumor_DEGs_symbols_biclusters.csv"

# Ruta de salida
ruta_salida = "C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/5_degree/tumor/FC_1"
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
    print(f"✅ Mapping bicluster encontrado con {max_matches} genes élite en común.")
    print("Genes en común:", ", ".join(genes_encontrados))

    output = pd.DataFrame([mejor_bicluster])
    output_path = os.path.join(ruta_salida, "mapping_bicluster_tumor.csv")
    output.to_csv(output_path, index=False)
    print(f"🔽 Guardado como: {output_path}")
else:
    print("⚠️ No se encontró ningún bicluster con genes élite en común.")
