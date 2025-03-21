# IUPred promedio y desvio estandar para cada peptido de ProPD
#Toma archivos propd y archivos de scores IUPred para toda la sec de la proteina. 
#se crea un archivo csv con las columnas 'Accession', 'SeqStart', 'SeqStop', 'Mean', "Std" y  "Disorder_%" para cada acc conteniendo al peptido que fue hit en propd. 

import os as os
import pandas as pd
import numpy as np
from datetime import datetime

# Obtener la fecha y hora actual para nombrar el archivo de salida
hoy = datetime.today()
hoy = hoy.strftime('%Y%m%d_%H%M%S')

# Definir las rutas de los directorios donde se encuentran los archivos
iupred_dir = "C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\IUPred\\MARZO_2023\\Accession.iupred_scores\\"
propd_dir = "C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\IlvaVSDavey\\ProPD_files\\Propd_ListaDefinitivaDePeptidos\\ProPD_PFAM\\"
path_out = "C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\IlvaVSDavey\\ProPD_files\\Propd_ListaDefinitivaDePeptidos\\IUPred scores_Propd\\"

# Iterar a través de los archivos en el directorio de ProPD
for archivo in os.listdir(propd_dir):
    # Filtrar archivos que comienzan con una fecha específica
    if archivo.startswith('20230330_150516'):
        # Construir la ruta completa del archivo de ProPD
        propd_files = propd_dir + archivo
        # Extraer el nombre del pocket del archivo
        pocket_name = archivo.replace('20230330_150516_PFam_', '')
        # Leer el archivo CSV de ProPD en un DataFrame
        df_propd = pd.read_csv(propd_files, sep=";")
        # Obtener la lista de accesiones de ProPD
        propd_acc = df_propd['Accession'].values
        
        # Crear un DataFrame vacío para almacenar los resultados
        df_pocket = pd.DataFrame(columns=[
            "Accession_number", "seq_Start", "seq_Stop", "Mean", "Std", "Disorder_%"])

        # Iterar sobre cada accesión
        for i, acc in enumerate(propd_acc):
            # Construir la ruta del archivo de scores IUPred para la accesión
            score_file = iupred_dir + acc + ".iupred"
            # Leer el archivo de scores IUPred en un DataFrame
            df_iupred = pd.read_csv(score_file, sep="\t", skiprows=6)
            
            # Obtener los índices de inicio y fin de la secuencia desde el DataFrame de ProPD
            seqStart = df_propd['Posicion_i'][i]
            seqStop = df_propd['Posicion_f'][i]
            
            # Almacenar la información de la accesión en el DataFrame de resultados
            df_pocket.loc[i, "Accession_number"] = acc
            df_pocket.loc[i, "seq_Start"] = seqStart
            df_pocket.loc[i, "seq_Stop"] = seqStop
            
            # Obtener los valores de IUPred para la región de interés
            score_iupred = df_iupred.loc[seqStart - 1: seqStop - 1, 'IUPRED2'].values
            
            # Calcular la media y la desviación estándar de los scores IUPred
            score_mean = score_iupred.mean()
            score_std = score_iupred.std()
            
            # Calcular el porcentaje de desorden en la secuencia
            porcentaje = sum(1 for x in score_iupred if x >= 0.4) * 100 / len(score_iupred)
            
            # Almacenar los resultados en el DataFrame
            df_pocket.loc[i, "Mean"] = score_mean
            df_pocket.loc[i, "Std"] = score_std
            df_pocket.loc[i, "Disorder_%"] = porcentaje
        
        # Construir el nombre del archivo de salida
        fileout = path_out + hoy + '_IUPred_scores_' + pocket_name
        # Guardar el DataFrame de resultados en un archivo CSV
        df_pocket.to_csv(fileout, header=True, index=None, sep=';', mode='w')