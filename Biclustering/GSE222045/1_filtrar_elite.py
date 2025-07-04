# Función modificada: conserva biclusters completos si al menos uno de sus genes está en la lista de élite
def conservar_biclusters_si_alguno_elite(input_path, elite_set, output_suffix="_elite_any"):
    df = pd.read_csv(input_path)
    
    def contiene_gen_elite(genes_str):
        genes = genes_str.split(";")
        return any(g in elite_set for g in genes)

    # Filtrar biclusters donde al menos uno de los genes es de la lista élite
    df_filtrado = df[df["Genes"].map(contiene_gen_elite)].reset_index(drop=True)

    output_path = input_path.replace(".csv", f"{output_suffix}.csv")
    df_filtrado.to_csv(output_path, index=False)
    return output_path

# Aplicar esta versión modificada al archivo único disponible
bicluster_file = "/mnt/data/normal_DEGs_symbols_biclusters.csv"
output_path = conservar_biclusters_si_alguno_elite(bicluster_file, elite_set)
output_path