import pandas as pd

# === Parámetros ===
tolerancia = 5  # Por ejemplo, +/-100 en grado
ruta_grados = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/5_degree/normal/FC_1/grados_todos_genes_bicluster.csv"
ruta_elite = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE222045/3_Filtrados_elite/genes_elite.txt"

# === Leer grados ya calculados ===
grados_df = pd.read_csv(ruta_grados)
grados_df.columns = grados_df.columns.str.strip()  # limpieza

# Leer lista de genes élite
with open(ruta_elite, "r") as f:
    elite_genes = [g.strip() for g in f.read().replace("\n", "").split(",")]

# Obtener grados de genes élite presentes
elite_df = grados_df[grados_df["Gene"].isin(elite_genes)]

# Buscar genes con grado dentro de la tolerancia
matched = []

for _, elite_row in elite_df.iterrows():
    elite_gene = elite_row["Gene"]
    elite_degree = elite_row["Degree"]
    lower = elite_degree - tolerancia
    upper = elite_degree + tolerancia

    genes_similares = grados_df[
        (grados_df["Degree"] >= lower) &
        (grados_df["Degree"] <= upper) &
        (grados_df["Gene"] != elite_gene)
    ].copy()

    genes_similares["EliteGeneReference"] = elite_gene
    genes_similares["EliteDegree"] = elite_degree
    matched.append(genes_similares)

# Guardar resultados
if matched:
    matched_df = pd.concat(matched, ignore_index=True)
    matched_df.to_csv("genes_mismo_grado_elite_mapping_bicluster_tumor.csv", index=False)
    print(f"\n✅ Script completado. Genes encontrados: {matched_df.shape[0]}")
else:
    print("\n⚠️ No se encontraron genes con grados coincidentes dentro de la tolerancia.")
