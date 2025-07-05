import networkx as nx
import csv

# Carga de datos (Marc)
normalSarcoma = nx.read_edgelist('C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/3venny/FC_3/normal_network_07_Symbol_NameInteraction.csv', delimiter=(';'), nodetype=str)
tumorSarcoma = nx.read_edgelist('C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/3venny/FC_3/tumor_network_07_Symbol_NameInteraction.csv', delimiter=(';'), nodetype=str)

# Carga de datos (aurelio)
#normalSarcoma = nx.read_edgelist('/Users/aurelio/git-repositories/sarcoma-iwbbio2025/GSE99671/methodology/3venny/network_normal_07_NameInteraction.csv', delimiter=(';'), nodetype=str)
#tumorSarcoma = nx.read_edgelist('/Users/aurelio/git-repositories/sarcoma-iwbbio2025/GSE99671/methodology/3venny/network_tumor_07_NameInteraction.csv', delimiter=(';'), nodetype=str)

# Invertir nodos del normal
invertedNormal = [(target, source) for source, target in normalSarcoma.edges()]
invertedNormal = nx.Graph(invertedNormal)

# Restar tumor menos el normal
interactionTumorSarcomaExclussive = set(tumorSarcoma.edges()) - set(normalSarcoma.edges()) # interacciones exclusivas tumorales
tumorSarcomaExclussive = nx.Graph(list(interactionTumorSarcomaExclussive))

# Restar resultado anterior menos el invertido
interactionTumorSarcomaExclussive = set(tumorSarcomaExclussive.edges()) - set(invertedNormal.edges()) # interacciones exclusivas tumorales
tumorSarcomaExclussive = nx.Graph(list(interactionTumorSarcomaExclussive))

# Extraer las aristas de tumorSarcomaExclussive
edges_tumor_exclusive = list(tumorSarcomaExclussive.edges())

# Guardar las aristas en un fichero CSV (Aurelio)
"""with open('/Users/aurelio/git-repositories/sarcoma-iwbbio2025/GSE99671/methodology/3venny/tumorExclussive.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["GEN1", "GEN2"])
    writer.writerows(edges_tumor_exclusive)
"""

# Guardar las aristas en un fichero CSV (Marc)
with open('C:/Users/marcr/OneDrive/Documentos/github_repositorios/sarcoma-iwbbio2025/GSE222045/results/3venny/FC_3/tumorExclussive.csv', mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["GEN1", "GEN2"])
    writer.writerows(edges_tumor_exclusive)
