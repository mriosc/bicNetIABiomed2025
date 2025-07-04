import os
import pandas as pd

# Cambiar al directorio de trabajo
os.chdir("C:/Users/marcr/OneDrive/Escritorio/IABioMed/Biclusters UniBic/GSE64763")

# Leer lista de genes élite
with open("3_Filtrados_elite/genes_elite.txt", "r") as f:
    elite_set = set(line.strip().upper() for line in open("3_Filtrados_elite/genes_elite.txt") if line.strip())

# Función modificada: conserva biclusters completos si al menos uno de sus genes está en la lista de élite
def conservar_biclusters_si_alguno_elite(input_path, elite_set, output_suffix="_elite"):
    df = pd.read_csv(input_path)

    def contiene_gen_elite(genes_str):
        genes = [g.strip().upper() for g in genes_str.split(";")]
        return any(g in elite_set for g in genes)

    df_filtrado = df[df["Genes"].map(contiene_gen_elite)].reset_index(drop=True)

    output_path = input_path.replace(".csv", f"{output_suffix}.csv")
    df_filtrado.to_csv(output_path, index=False)
    return output_path


# Aplicar esta versión modificada al archivo único disponible
bicluster_file = "2_biclusters/FC_1/FC1_tumor_DEGs_symbols_deduplicated_biclusters.csv"
output_path = conservar_biclusters_si_alguno_elite(bicluster_file, elite_set)
print("Archivo filtrado guardado en:", output_path)
