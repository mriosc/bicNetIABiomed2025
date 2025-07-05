import pandas as pd
import matplotlib.pyplot as plt
import os

"""
Los biólogos ya han establecido que si el valor de expresión genética de un gen particular aumenta (o disminuye) e
n todas las muestras mientras pasa del estado normal al estado patológico, entonces este gen específico tiene un extremo probabilidad 
convertirse en un biomarcador potencial de esa enfermedad (Mandal y otros, 2020; Xu y otros, 2006; Kakati y otros, 2018). 
Ahora aplicamos esta regla de la siguiente manera: consideramos los genes obtenidos de los pasos 12 y 14 (anteriormente considerados como 
biomarcadores potenciales sospechosos) e inspeccionamos si los valores de expresión de estos genes varían significativamente (aumentan o disminuyen) 
para todo tipo de muestras.

Los genes que muestran la respuesta positiva en este ultimo paso pueden considerarse biomarcadores potenciales.
"""

# Cargar datos
#Marc
#ruta_normal = "/mnt/data/FC1_normal_DEGs_symbols.csv"
#ruta_tumor = "/mnt/data/FC1_tumor_DEGs_symbols.csv"
#ruta_biomarcadores = "/mnt/data/potenciales_biomarcadores_por_elite.csv"

# Aurelio
ruta_normal = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/1_DEGs/FC_3/normal_DEGs_symbols.csv"
ruta_tumor = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/1_DEGs/FC_3/tumor_DEGs_symbols.csv"
ruta_biomarcadores = r"/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/6_final_biomarkers/FC_3/potenciales_biomarcadores_por_elite.csv"


df_normal = pd.read_csv(ruta_normal, sep="\t", index_col=0)
df_tumor = pd.read_csv(ruta_tumor, sep="\t", index_col=0)
df_biomarcadores = pd.read_csv(ruta_biomarcadores)

# Crear carpeta de salida
output_dir = "/home/principalpc/git-repositories/bicNetIABiomed2025/Biclustering/GSE222045/7_expression_biomarkers/FC_3"
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