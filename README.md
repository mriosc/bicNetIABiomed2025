# IAbiomed

### GSE64763
Tipo de estudio: Microarray de expresión génica (Affymetrix U133A 2.0)
Muestras:
25 Tumores uterinos (ULMS)
24 Fibroides (benignos)
30 Miometrio normal

### GSE222045
Tipo de estudio: RNA-seq
Muestras:
45 Tumores de osteosarcoma
45 Tejidos sanos (de los mismos pacientes)


## Metodología

### PASO 6,7,8)
(3_script_matriz_degree.py)
En el archivo genes_mismo_grado_elite.csv faltaría por hacer lo siguiente:

- Meter el parámetro de % coincidencia para los grados.
SHMT1: 53 (Faltaría guardar en este CSV con que gen elite está relacionado).

### PASO 9) 
Ejecutar algoritmo de biclustering para dataset tumoral.** (Ejecutado al principio)

### PASO 10) 
(script: 4_mapping.py)
Coge el mejor bicluster del normal y busca si TODOS sus genes élites están presentes en algún bicluster tumoral. Si es así, el bicluster sería un bicluster de mapeo. Objetivo: Reducir el numero de bicluster del tumor. 

### PASO 11 + PASO13) 
(5_tumor_matriz_degree.py)
Para el bicluster tumoral, ejecutar el mismo script del paso 6,7,8.


GENE2: 60 <-> GE1 - ASPM.

¿Con cual te quedas? Con el que tenga el mismo grado, es decir, solo el Gene6.

¿Por qué? Debido a la variabilidad de los valores de expresión, los genes que tienen el mismo grado que un gen de élite en una condición normal pueden comportarse de manera diferente en el estado de enfermedad.


### PASO 12 (REPETIR POR CADA UNO DE LOS GENES ÉLITES)

GE1 - ASPM
################

1) Comparar el orden del gen elite G1.-ASPM de normal hacia tumoral:
GE1 - ASPM (normal): 77
GE1 - ASPM (tumoral): 91
Sentido: Aumentando.

2) Quedarnos con aquellos genes que estén presentes tanto en bicluster normal como en tumoral. Es decir, si hay genes que están SOLO en el bicluster normal o en el tumoral, se eliminan.

3) A partir de aquellos genes que están presentes en los dos biclusters (normal y tumoral), hacer la siguiente operación:

Comprobar el orden de cada gen si hay cambio:

BICLUSTER NORMAL / BICLUSTER TUMORAL:
GENE1: 68 / 70 (aum)
GENE2: 54 / 60 (aum)
GENE3: 80 / 70 (dis)
GENE4: 26 / 40 (aum)
GENE5: 67 / 60 (dis)
GENE6: 72 / 85 (aum)

Potenciales biomarcadores: Gene1, Gene2, Gene4, Gene6 porque tienen el mismo sentido (aumento) que el gen elite G1-ASPM.

**Hacer el mismo proceso con cada uno de los genes élites.**

---------------------------------------------------------

### PASO 15)

Potencial biomarcador del paso 12: Gene1, Gene2, Gene4, Gene6
Potencial biomarcador del paso 14: Gene6.

Potenciales biomarcadores totales: Gene1, Gene2, Gene4, Gene6

#######################################################

Pintar la gráfica 3 del paper con los potenciales biomarcadores totales (nivel de expresión génica de los DEGs).

PASO 16) Ir a GO/KEGG, bibliografía, para comprobar si realmente tienen sentido como potencial biomarcador.


