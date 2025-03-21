#Script para agregar informacion en archivos PROPD de instancias de motivos anotadas en ELM databse en 3 columnas, segun si el peptido fue identificado con el motivo LxCxE o E2F_mimic o LxSxE.
#Los peptidos que figuran en libreria HD2 "heredan" esta informacion cuando se cruzan datos de los archivos para completas Propd. Sin embargo, aquellos peptidos de libreria 'HD' o 'NA' no contaran con esta informacion. Con este script se completa la informacion de esos peptidos que de otra manera, no serian analizados. 
#El umbral considerado sera de 9 residuos de overlap (56% de identidad de secuencia) dado que encontramos true positives con un overlap de 9 residuos. 
#FEBRERO 2023: CAMBIO EL UMBRAL A 8 RESIDUOS, PORQUE E2F3 (QUE ES UN TP) COINCIDE CON LA SECUENCIA DE ELM EN 8 RESIDUOS.                 
# propdViejos= "C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\Ilva vs Davey\\ProPD_files\\ProP-PD_DAVEY HD2\\Columnas_'counts'_sumadas\\20220805_164005_p107_hd2.csv"
# with open(propdViejos) as oldFiles:
#     oldies_df= pd.read_csv(oldFiles, sep=';')
#     print(oldies_df.columns)


import os  # Proporciona funciones para interactuar con el sistema operativo, como leer archivos en directorios.
import sys  # Se usa para acceder a funciones y variables del sistema (aunque en este código no se usa explícitamente).
import pandas as pd  # Librería para manipulación y análisis de datos, especialmente en estructuras tipo DataFrame.
import numpy as np  # Biblioteca para trabajar con arreglos numéricos y operaciones matemáticas (no se usa en este código).
import difflib  # Permite comparar secuencias y encontrar similitudes entre cadenas de texto.
from datetime import datetime  # Manejo de fechas y tiempos en Python.


# Obtiene la fecha y hora actual en el formato 'YYYYMMDD_HHMMSS'
hoy = datetime.today()
hoy = hoy.strftime('%Y%m%d_%H%M%S')

# Función para calcular el solapamiento entre dos secuencias
# Devuelve True si la proporción de solapamiento de s2 en s1 es mayor o igual al umbral
def get_overlap(s1, s2, umbral):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    # Verifica si el segmento solapado está en un extremo de alguna de las secuencias
    if (pos_b == 0 and (pos_a + size == len(s1))) or ((pos_b + size == len(s2)) and pos_a == 0):
        overlap = ((size/len(s2)) >= umbral)
    else:
        overlap = False
    return overlap

# Función para marcar la presencia de motivos específicos en un DataFrame
def completeHD2(df1, i, df2, j):
    if df1.at[i, 'ELMIdentifier'] == 'LIG_Rb_LxCxE_1':
        df2.at[j, 'LxCxE'] = True
    elif df1.at[i, 'ELMIdentifier'] == 'LIG_Rb_pABgroove_1':
        df2.at[j, 'E2F_mimic'] = True
    elif df1.at[i, 'ELMIdentifier'] == 'LxSxE':
        df2.at[j, 'LxSxE'] = True
    return df2

# Definición de rutas de archivos de entrada y salida
ELMdir = 'C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\files_ELM\\elm_instances+motivos.tsv'
propdDir = 'C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\IlvaVSDavey\\ProPD_files\\Propd_ListaDefinitivaDePeptidos'
path_out = 'c:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\IlvaVSDavey\\ProPD_files\\Propd_ListaDefinitivaDePeptidos\\ProPD_conELM'

# Carga el archivo de instancias ELM
with open(ELMdir) as ELMinstances:
    df_ELM = pd.read_csv(ELMinstances, skiprows=5, sep='\t')

# Itera sobre los archivos de ProPD en el directorio especificado
for file in os.listdir(propdDir):
    if file.endswith('.csv'):
        propd_file = propdDir + '\\' + file
        
        # Carga el archivo CSV en un DataFrame
        propd_df = pd.read_csv(propd_file, sep=';')
        
        # Inicializa las columnas de los motivos como False
        propd_df['LxCxE'] = False
        propd_df['E2F_mimic'] = False
        propd_df['LxSxE'] = False

        # Itera sobre las entradas de ELM para identificar coincidencias
        for k, acc_elm in enumerate(df_ELM['Primary_Acc']):
            if df_ELM.at[k, 'InstanceLogic'] == 'true positive':
                seq_ELM = df_ELM.at[k, 'Sequence']
                
                # Si el identificador de acceso está en el DataFrame de ProPD
                if acc_elm in propd_df["Accession"].values:
                    indices_propd = propd_df.loc[propd_df["Accession"] == acc_elm].index
                    print(indices_propd, acc_elm, file)
                    
                    # Itera sobre los índices coincidentes
                    for x in indices_propd:
                        seq_propd = propd_df.at[x, 'Sec_C']
                        
                        # Verifica solapamiento con umbral del 50% (8 residuos)
                        if get_overlap(seq_ELM, seq_propd, 0.50):
                            propd_df = completeHD2(df_ELM, k, propd_df, x)

        # Guarda el archivo con la información agregada
        propd_File = file.replace('10323_ProPdcompl_PI_PF_', 'ELM_')
        fileout = path_out + '\\' + hoy + '_' + propd_File
        propd_df.to_csv(fileout, header=True, index=None, sep=';', mode='w')

