import pandas as pd
import mygene

def convert_to_symbol(column_name, df, mg_instance):
    df[column_name] = df[column_name].astype('str')
    attr_names = df[column_name].tolist()
    query_result = mg_instance.querymany(attr_names, scopes='entrezgene', fields='symbol', species='human')
    gene_symbols = pd.DataFrame(query_result)
    gene_symbols['query'] = gene_symbols['query'].astype(str)
    gene_symbols['symbol'] = gene_symbols['symbol'].astype(str)
    df_with_symbols = df.merge(gene_symbols[['query', 'symbol']], left_on=column_name, right_on='query', how='left')
    df_with_symbols[column_name] = df_with_symbols['symbol']
    
    if column_name == 'interaction':
        df_with_symbols['shared interaction'] = df_with_symbols['symbol']
    elif column_name == 'name':
        df_with_symbols['shared name'] = df_with_symbols['symbol']    
    
    df_with_symbols.drop(columns=['query', 'symbol'], inplace=True)
    return df_with_symbols

#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_07.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_075.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_08.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_085.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_09.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_095.csv',sep=',')


#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_07.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_075.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_08.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_085.csv',sep=',')
#data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_09.csv',sep=',')
data = pd.read_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_095.csv',sep=',')

# 1) Convertir genes a Symbol (interaction column)
mg = mygene.MyGeneInfo()
columns_to_convert = ['interaction', 'name']
for col in columns_to_convert:
    data = convert_to_symbol(col, data, mg)

# 2) Eliminar aquellas filas que no se han encontrado Symbol.
data = data.dropna(subset=['interaction', 'name'])
data = data[data['interaction'] != 'nan']
data = data[data['name'] != 'nan']

# 3) Identificar Symbol duplicados y quedarnos con el maximo valor de expresion
data = data.sort_values(by='Weight', ascending=False).drop_duplicates(subset=['interaction', 'name'], keep='first')



#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_07_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_075_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_08_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_085_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_09_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/tumor/tumor_network_affy_095_Symbol.csv', sep=',', index=False)

#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_07_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_075_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_08_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_085_Symbol.csv', sep=',', index=False)
#data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_09_Symbol.csv', sep=',', index=False)
data.to_csv('C:/Users/aurel/git-repositories/sarcoma-iwbbio2025/GSE222045/results/2networks/original_networks/FC_3/normal/normal_network_affy_095_Symbol.csv', sep=',', index=False)


