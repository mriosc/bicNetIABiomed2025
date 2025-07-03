import pandas as pd
### EJECUTAR ESTE TROZO DE CÓDIGO TAN SÓLO UNA VEZ
# Eliminar múltiples tabuladores juntos
with open('c:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/methodology/1getExpressionData/results/aafTableAnn_majorRevision.txt', 'r') as file:
    content = file.read()

content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')
content = content.replace('\t\t', '\t')

with open('c:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE17674/methodology/1getExpressionData/results/aafTableAnn_majorRevision_no_tabs.txt', 'w') as file:
    file.write(content)