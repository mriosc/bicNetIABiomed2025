import pandas as pd
import os

# Load the CSV file Aurelio
# df = pd.read_csv('/Users/aurelio/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/4clustermaker/clusters.csv') # Aurelio MacOS
# df = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/4clustermaker/clusters.csv') # Aurelio Windows
# output_path = 'C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/5cytohubba/results/filtered_clusters.csv' # Aurelio Windows
# output_txt_dir = 'C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/6geneEnrichment/input' # Aurelio Windows

# Load the CSV file Marc
df = pd.read_csv('C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/4clustermaker/FC_3/clusters.csv') # Marc Git
output_path = 'C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/5cytohubba/FC_3/filtered_clusters.csv' # Marc Git
output_txt_dir = 'C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/6geneEnrichment/FC_3/input' # Marc Git
#df = pd.read_csv('C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/5cytohubba/clusters.csv') # Marc prueba
#output_path = 'C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/5cytohubba/results/filtered_clusters.csv' # Marc prueba
#output_txt_dir = 'C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/6geneEnrichment/input' # Marc prueba

# Display the first few rows of the dataframe to understand its structure
df.head()

# Agrupar por la columna '__glayCluster' y eliminar grupos con menos de 10 filas
filtered_df = df.groupby('__glayCluster').filter(lambda x: len(x) >= 10)

# Guardar el dataframe filtrado en un nuevo archivo CSV

filtered_df.to_csv(output_path, index=False)

print(f"Archivo guardado en: {output_path}")

# Guardar un archivo .txt por cada cluster, con los valores de la columna 'name'
for cluster_value, group in filtered_df.groupby('__glayCluster'):
    # Obtener los valores de la columna 'name'
    cluster_names = group['name']
    
    # Crear el nombre del archivo .txt
    txt_file_path = os.path.join(output_txt_dir, f'Cluster_{cluster_value}.txt')
    
    # Guardar los nombres en el archivo .txt sin cabecera
    cluster_names.to_csv(txt_file_path, index=False, header=False)

print(f"Archivos guardados en: {output_txt_dir}")