import pandas as pd

# === Rutas de entrada ===
ruta_grados_normal = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/5_degree/normal/FC_1/grados_todos_genes_bicluster.csv"
ruta_grados_tumor = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/5_degree/tumor/FC_1/grados_mapping_bicluster_tumor.csv"
ruta_elite = r"C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763/3_Filtrados_elite/genes_elite.txt"

# === Leer archivos ===
grados_normal = pd.read_csv(ruta_grados_normal)
grados_tumor = pd.read_csv(ruta_grados_tumor)

# Asegurar nombres de columnas
grados_normal.columns = grados_normal.columns.str.strip()
grados_tumor.columns = grados_tumor.columns.str.strip()

# Leer genes élite
with open(ruta_elite, "r") as f:
    genes_elite = [g.strip() for g in f.read().replace("\n", "").split(",")]

# Intersección de genes entre ambas condiciones
genes_comunes = set(grados_normal["Gene"]).intersection(set(grados_tumor["Gene"]))

# Filtrar solo genes comunes
grados_normal = grados_normal[grados_normal["Gene"].isin(genes_comunes)]
grados_tumor = grados_tumor[grados_tumor["Gene"].isin(genes_comunes)]

# === Resultado final ===
resultados = []

for gen_elite in genes_elite:
    # Obtener grados del gen élite
    grado_normal = grados_normal[grados_normal["Gene"] == gen_elite]["Degree"]
    grado_tumor = grados_tumor[grados_tumor["Gene"] == gen_elite]["Degree"]

    if grado_normal.empty or grado_tumor.empty:
        continue  # Si no hay datos en alguna condición, saltar

    sentido = "aumenta" if grado_tumor.values[0] > grado_normal.values[0] else "disminuye"

    # Fusionar grados por gene
    grados_merge = pd.merge(grados_normal, grados_tumor, on="Gene", suffixes=("_normal", "_tumor"))

    # Eliminar el gen élite del análisis
    grados_merge = grados_merge[grados_merge["Gene"] != gen_elite]

    # Calcular sentido de cambio para cada gen
    grados_merge["Sentido"] = grados_merge.apply(
        lambda row: "aumenta" if row["Degree_tumor"] > row["Degree_normal"]
        else "disminuye" if row["Degree_tumor"] < row["Degree_normal"]
        else "igual",
        axis=1
    )

    # Filtrar genes con el mismo sentido que el gen élite
    biomarcadores = grados_merge[grados_merge["Sentido"] == sentido].copy()
    biomarcadores["Elite"] = gen_elite
    biomarcadores["SentidoElite"] = sentido

    resultados.append(biomarcadores[["Gene", "Degree_normal", "Degree_tumor", "Sentido", "Elite", "SentidoElite"]])

# Guardar resultado
if resultados:
    df_final = pd.concat(resultados, ignore_index=True)
    df_final.to_csv("potenciales_biomarcadores_por_elite.csv", index=False)
    print("✅ CSV generado con biomarcadores potenciales:", df_final.shape[0])
else:
    print("⚠️ No se encontraron biomarcadores con cambio coherente.")
