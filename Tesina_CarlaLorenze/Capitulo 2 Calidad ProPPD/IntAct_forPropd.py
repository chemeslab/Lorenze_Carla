#Script para agregar columnas de informacion IntAct a los archivos Propd que contienen datos de las dos librerias completas (HD y HD2)


#Creo una columna "HD2_protein" de TRUE/FALSE en cada uno de los archivos intact si el interactor de la lista intact esta en HD2, chequeo por accession number 
#df_intact['HD2_protein'] = newIntact.isin(df_hd2['Accession_n'].values)

from hashlib import new
import os
import sys
import pandas as pd
import numpy as np
from datetime import datetime

hoy = datetime.today()
hoy = hoy.strftime('%Y%m%d_%H%M%S')

#archivos con umbral 56% y con columnas ELM
propd_dir = "C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\IlvaVSDavey\\ProPD_files\\Propd_ListaDefinitivaDePeptidos\\ProPD_conELM\\"
intact_raw = 'c:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\RawData\\IntAct\\'
intact_lines = 'IntAct_files.txt'
pathOUT = "C:\\Users\\Carla\\Desktop\\Tesina\\GitHub\\phageD\\Data\\IlvaVSDavey\\ProPD_files\\Propd_ListaDefinitivaDePeptidos\\ProPD_conIntAct"

#Abro los archivos propd
for file in os.listdir(propd_dir):
    if file.startswith("20230303_130308_"):
        archivo= propd_dir + "\\" + file
        #print(archivo)
        propd_df = pd.read_csv(archivo, sep=';')
        #print(propd_df.head())

    #     #abro los archivos de datos descargados de intact para cada pocket, int directa o indirecta.
        with open(intact_raw + intact_lines) as intact:
            for line in intact:
                #Cada 'line' es el nombre con el que identificamos los archivos intact (Rb directo/indirecto, p107 directo/indirecto, p130 directo/indirecto)
                #rstrip elimina los espacios finales del line
                line = line.rstrip()
                id_target = ''
                # Identifico en 'id_target 'los acc number de las pocket proteins
                if line == 'Rb_Indirect' or line == 'Rb_Direct':
                    id_target = 'uniprotkb:P06400'
                elif line == 'p107_Indirect' or line == 'p107_Direct':
                    id_target = 'uniprotkb:P28749'
                elif line == 'p130_Indirect' or line == 'p130_Direct':
                    id_target = 'uniprotkb:Q08999'
                #Abro los archivos intact y los paso a DF
                archivos_intact = intact_raw + line + '.txt'
                df_intact = pd.read_csv(archivos_intact, sep='\t')
                newIntact = []
                for i in range(0, len(df_intact['# ID(s) interactor A'])):
                    # Este for recorre los interactores para ver si el interactor A o el B es la bait (pocket protein) o no porque en intact A y B están desordenados y no siempre A es la pocket o viceversa.
                    if df_intact['# ID(s) interactor A'][i].startswith('intact:') or df_intact['ID(s) interactor B'][i].startswith('intact:'):
                        newIntact.append("-")
                    else:
                        if df_intact['# ID(s) interactor A'][i] != id_target:
                            newIntact.append(df_intact['# ID(s) interactor A'][i])
                        else:
                            newIntact.append(df_intact['ID(s) interactor B'][i])
                newIntact = pd.Series(newIntact).str.replace('uniprotkb:',"",regex=False).str.replace(r'-[0-9]',"",regex=True)
                #uso los accession numbers de los files propd con la funcion explode para que separe los que estan en un mismo peptido separado con coma 
                #explode_df =propd_df['Accession'].str.split(",").explode()
                propd_df[line]= False
                for intact_acc in newIntact.values:
                    if intact_acc in propd_df['Accession'].values:
                        indices_propd = propd_df.loc[propd_df['Accession'] == intact_acc].index
                        for k in indices_propd:
                            #Completo columna en propd por cada archivo intact (que lleva el nombre de "pocket_interaccion") de TRUE/FALSE si el acc de propd esta en la lista de interactores que descargamos de intact
                            propd_df.at[k, line] = True
    #print(propd_df['p107_Indirect'].any())

        file = file.replace("20230303_130308_ELM", "_IntAct")
        fileout = pathOUT + "\\" + hoy +  file 
        propd_df.to_csv(fileout, header=True, index=None, sep=';', mode='w')