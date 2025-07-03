import pandas as pd

# Reemplazar Affymetrix ID por Symbol en la red
df_csv = pd.read_csv('c:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/results/2networks/original_networks/FC_3/normal/normal_network_affy_095.csv')  # Cargar el archivo csv_5filas.csv
df_txt = pd.read_csv('c:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/results/1getExpressionData/aafTableAnn_majorRevision_no_tabs.txt', sep='\t')  # Cargar el archivo aafTableAnn_majorRevision.txt con tabuladores

probe_to_symbol = dict(zip(df_txt['Probe'], df_txt['Symbol']))
columns_to_rename = ['interaction', 'name', 'shared interaction', 'shared name']
for col in columns_to_rename:
    df_csv[col] = df_csv[col].map(probe_to_symbol).fillna(df_csv[col])

# Eliminar aquellas interacciones en el que se relacionen consigo mismo, es decir, en el que la columna interaction y name tengan el mismo Symbol.
df_csv = df_csv[df_csv['interaction'] != df_csv['name']]

# Eliminar interacciones duplicadas (interaction y name) y nos quedamos con aquellas interacciones con mayor weight.
df_csv['interaction_pair'] = df_csv.apply(lambda row: tuple(sorted([row['interaction'], row['name']])), axis=1)

# Ahora agrupamos por la nueva columna 'interaction_pair' y nos quedamos con la fila que tenga el valor máximo de 'weight'
df_csv_deduplicated = df_csv.loc[df_csv.groupby('interaction_pair')['Weight'].idxmax()]

# Eliminar la columna auxiliar 'interaction_pair'
df_csv_deduplicated = df_csv_deduplicated.drop(columns=['interaction_pair'])

# Guardar fichero.
df_csv_deduplicated.to_csv('c:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/results/2networks/results_networks/FC_3/normal/normal_network_affy_095.csv', index=False)

print("Renombrado completado y guardado en 'csv_5filas_renombrado.csv'.")