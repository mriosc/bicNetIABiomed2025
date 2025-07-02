import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os



txt_clusters_dir = 'C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/6geneEnrichment/FC_3/input'
output_path = 'C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/5cytohubba/FC_3/infoHubs.txt'

#txt_clusters_dir = 'C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/6geneEnrichment/input'
#output_path = 'C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/5cytohubba/results/infoHubs.txt'

# txt_clusters_dir = 'C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/6geneEnrichment/input'  # Ruta donde se generaron los archivos .txt de los clusters
# output_path = 'C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/5cytohubba/results/infoHubs.txt'

df = pd.read_csv('C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/3venny/FC_3/tumorExclussive.csv', delimiter= ',') # Marc Git
# df = pd.read_csv('C:/Users/marcr/OneDrive/Escritorio/iwbbio/GSE99671/methodology/3venny/tumorExclussive.csv', delimiter= ',') # Marc prueba
# df = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE17674/methodology/3venny/tumorExclussive.csv', delimiter=',') # Aurelio Windows
G = nx.Graph()

for index, row in df.iterrows():
    G.add_edge(row['GEN1'], row['GEN2'])
    
degrees = [d for n, d in G.degree()]


# calculo de hub (FER. art.)
Q1 = np.percentile(degrees, 25)
Q3 = np.percentile(degrees, 75)
IQR = Q3 - Q1

# Calcular el umbral para detectar hubs
threshold = Q3 + 1.5 * IQR

# Identificar los nodos hub
hubs = [node for node, deg in G.degree() if deg > threshold]

print(f"Umbral de grado para hubs: {threshold}")
print(f"Nodos hub: {hubs}")
print(f"Total hubs: {len(hubs)}")

# Buscar en qué cluster están los hubs
hub_clusters = {}

# Iterar sobre cada archivo de cluster y buscar si los hubs están en ese cluster
for txt_file in os.listdir(txt_clusters_dir):
    if txt_file.endswith('.txt'):
        cluster_value = txt_file.split('_')[1].split('.')[0]  # Obtener el número de cluster del nombre del archivo
        with open(os.path.join(txt_clusters_dir, txt_file), 'r') as f:
            cluster_genes = f.read().splitlines()
            for hub in hubs:
                if hub in cluster_genes:
                    if hub not in hub_clusters:
                        hub_clusters[hub] = []
                    hub_clusters[hub].append(cluster_value)

# Guardar la información en un archivo .txt
with open(output_path, 'w') as f:
    f.write(f"Nodos hub: {hubs}\n")
    f.write(f"Total hubs: {len(hubs)}\n\n")
    f.write("Hub - Clusters:\n")
    for hub, clusters in hub_clusters.items():
        f.write(f"{hub}: {' '.join(clusters)}\n")

print(f"Información de los hubs guardada en: {output_path}")

# Visualización opcional de la red y los hubs
pos = nx.spring_layout(G)
plt.figure(figsize=(12, 12))
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_nodes(G, pos, nodelist=G.nodes(), node_size=20, node_color='blue')
nx.draw_networkx_nodes(G, pos, nodelist=hubs, node_size=50, node_color='red')
nx.draw_networkx_labels(G, pos, labels={node: node for node in hubs}, font_size=12)
plt.title("Red con nodos hub destacados")
plt.show()
