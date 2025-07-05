import pandas as pd
import matplotlib.pyplot as plt
import os

# Cargar datos
ruta_normal = "/mnt/data/FC1_normal_DEGs_symbols_deduplicated.csv"
ruta_tumor = "/mnt/data/FC1_tumor_DEGs_symbols_deduplicated.csv"
ruta_biomarcadores = "/mnt/data/potenciales_biomarcadores_por_elite.csv"

df_normal = pd.read_csv(ruta_normal, sep="\t", index_col=0)
df_tumor = pd.read_csv(ruta_tumor, sep="\t", index_col=0)
df_biomarcadores = pd.read_csv(ruta_biomarcadores)

# Crear carpeta de salida
output_dir = "/mnt/data/expresion_genes_png"
os.makedirs(output_dir, exist_ok=True)

# Obtener lista única de genes biomarcadores
genes_interes = df_biomarcadores["Gene"].unique()

# Crear gráfico para cada gen
for gene in genes_interes:
    if gene in df_normal.index and gene in df_tumor.index:
        plt.figure(figsize=(6, 4))
        plt.plot(range(df_normal.shape[1]), df_normal.loc[gene].values, marker='o', label='Normal')
        plt.plot(range(df_tumor.shape[1]), df_tumor.loc[gene].values, marker='o', label='Tumor')
        plt.title(f"Gene: {gene}")
        plt.xlabel("Sample ID")
        plt.ylabel("Gene Expression Value")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{gene}.png")
        plt.close()

import zipfile

# Comprimir las imágenes en un ZIP
zip_path = "/mnt/data/expresion_genes_png.zip"
with zipfile.ZipFile(zip_path, 'w') as zipf:
    for root, _, files in os.walk(output_dir):
        for file in files:
            file_path = os.path.join(root, file)
            zipf.write(file_path, arcname=file)

zip_path