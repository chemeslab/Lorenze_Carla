
import pandas as pd  # Librería para manipulación y análisis de datos en estructuras tipo DataFrame.
import numpy as np  # Biblioteca para operaciones numéricas y manejo de arreglos.
from datetime import datetime  # Se usa para trabajar con fechas y horas en Python.

# Obtiene la fecha y hora actual en formato 'YYYYMMDD_HHMMSS', útil para generar nombres de archivos con timestamp.
hoy = datetime.today()
hoy = hoy.strftime('%Y%m%d_%H%M%S')

# Definición de rutas de archivos
fileHD2 = 'c:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\files_all\\HD2_lines.txt'
hd2_DupSeq_dir = 'C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\files_Proteomic\\'
hd2_DupSeq_file = '20230120_132823_HD2_Proteomic_DupSeq.csv'
HD2_Pathout = 'c:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\files_all\\'

# Inicialización de una variable para controlar la primera iteración del bucle
count = 0

# Apertura del archivo fileHD2 que contiene nombres de archivos a procesar
with open(fileHD2) as archivos:
    for archivo in archivos:  # Itera sobre cada línea del archivo (cada línea representa un archivo CSV)
        archivo = archivo.rstrip()  # Elimina espacios en blanco y saltos de línea al final de la cadena
        archivo = HD2_Pathout + archivo + '.csv'  # Construye la ruta completa del archivo CSV a leer
        
        # Carga el archivo CSV en un DataFrame de pandas
        df_small = pd.read_csv(archivo, sep=';')
        
        if count != 0:
            # Si no es la primera iteración, concatena el nuevo DataFrame con el DataFrame acumulador
            df_DB = pd.concat([df_DB, df_small], ignore_index=True)
        else:
            # Si es la primera iteración, inicializa df_DB con el primer DataFrame leído
            df_DB = df_small
            count = 1  # Cambia count para indicar que ya se ha leído el primer archivo

    # Crea una nueva columna 'UniqueID' combinando Accession_n, Posicion_inicio y Sec_pep (secuencia mutada, con A)
    df_DB['UniqueID'] = df_DB['Accession_n'].values + '_' + \
        [str(i) for i in df_DB['Posicion_inicio']] + \
        '_' + df_DB['Sec_pep'].values

# Lectura del archivo de secuencias duplicadas y procesamiento de datos
with open(hd2_DupSeq_dir + hd2_DupSeq_file) as HD2_DupSeq:
    df_DupSeq = pd.read_csv(HD2_DupSeq, sep=';')  # Carga el archivo en un DataFrame

    # Agrupa por la columna 'Sec_db' (secuencia sin mutar, con C) y concatena valores únicos separados por comas
    df_grouped = df_DupSeq.groupby(['Sec_db'], as_index=False, sort=False).agg(
        lambda x: ','.join(np.unique([str(i) for i in x])))

    # Cuenta la cantidad de veces que cada 'Sec_db' aparece en 'Accession_n'
    df_counts = df_DupSeq.groupby(['Sec_db'], as_index=False, sort=False)['Accession_n'].count()
    
    # Renombra la columna 'Accession_n' a 'Rep_secDB' para indicar la cantidad de repeticiones
    df_counts.rename(columns={"Accession_n": "Rep_secDB"}, inplace=True)

    # Restablece los índices del DataFrame agrupado (aunque no se guarda en la variable)
    df_grouped.reset_index(drop=True)

    # Une la información de conteo con el DataFrame agrupado
    df_grouped = df_grouped.join(df_counts['Rep_secDB'], how='outer')

# Termina la ejecución del script en este punto 
exit()

# Construcción del nombre del archivo de salida reemplazando el nombre base con '_HD2final.csv'
final_file = hd2_DupSeq_file.replace(
    '20230120_132823_HD2_Proteomic_DupSeq.csv', '_HD2final.csv')

# Ruta completa para guardar el archivo de salida
fileout = HD2_Pathout + '\\' + hoy + final_file

# Guarda el DataFrame final en un archivo CSV con el formato especificado
df_grouped.to_csv(fileout, header=True, index=None, sep=';', mode='w')