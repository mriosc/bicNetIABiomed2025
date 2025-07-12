import pandas as pd

# Cargar lista de genes élite desde un archivo .txt (separados por coma)
def cargar_genes_elite(path_txt):
    with open(path_txt, "r") as f:
        genes = [g.strip() for g in f.read().split(",")]
    return set(genes)

# Función principal: conservar bicluster si al menos un gen es élite
def conservar_biclusters_si_alguno_elite(input_csv, elite_genes, output_suffix="_elite_any"):
    df = pd.read_csv(input_csv)

    def contiene_gen_elite(genes_str):
        genes = genes_str.split(";")
        return any(g in elite_genes for g in genes)

    df_filtrado = df[df["Genes"].map(contiene_gen_elite)].reset_index(drop=True)
    output_path = input_csv.replace(".csv", f"{output_suffix}.csv")
    df_filtrado.to_csv(output_path, index=False)
    print(f"✅ Guardado: {output_path}")

# USO EJEMPLO
if __name__ == "__main__":
    # Cambia las rutas según tu entorno local
    path_genes_elite = "C:/Users/ivan0/OneDrive - Universidad Pablo de Olavide de Sevilla/Universidad/Investigación/Congreso/Datasets/genes_elite.txt"
    archivos_biclusters = [
        "C:/Users/ivan0/OneDrive - Universidad Pablo de Olavide de Sevilla/Universidad/Investigación/Congreso/Datasets/Biclusters Dataset Nuevo/FC1_normal_DEGs_symbols_deduplicated_biclusters.csv",
        # añade aquí más archivos si lo deseas
    ]

    elite_genes = cargar_genes_elite(path_genes_elite)

    for archivo in archivos_biclusters:
        conservar_biclusters_si_alguno_elite(archivo, elite_genes)
