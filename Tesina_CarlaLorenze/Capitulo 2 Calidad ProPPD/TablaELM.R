# Cargar la librería 'dplyr' que proporciona funciones para manipulación de datos
library(dplyr)

# Leer el archivo CSV con la lista de instancias ELM (Lista de base de datos ELM + instancias curadas manualmente)
elm <- read.csv(file= 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/RepositorioGitHub/ELM/ELM_Motivos.csv', 
                sep= ';', 
                header = T) 

# Leer el archivo TSV (tab-separated values) que contiene informacion completa de las instancias ELM, descargado de PDB utilizando los accession numbers. 
elmMapping <-read.csv(file= 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/RepositorioGitHub/ELM/ELMmapping.tsv', 
                      sep= '\t', 
                      header = T) 

# Filtrar los datos para mantener solo las filas donde la columna 'InstanceLogic' tenga el valor 'true positive'
elmHomo <- elm[elm$InstanceLogic == 'true positive',]

# Unir los datos de 'elmHomo' con 'elmMapping' usando la columna 'Primary_Acc' de 'elmHomo' y 'Entry' de 'elmMapping'
# Se seleccionan solo algunas columnas del archivo 'elmMapping' para agregar a 'elmHomo'
elmHomo <- elmHomo %>%
  left_join(elmMapping %>% select(Entry,
                                  Entry.Name, 
                                  Gene.Names, 
                                  Protein.names ), 
            by = c("Primary_Acc" = 'Entry')) # Establecer la correspondencia entre las columnas 'Primary_Acc' y 'Entry'

# Seleccionar las columnas de interés y renombrarlas para tener nombres más claros
elmHomo <- elmHomo %>%
  select(
    UniprotID = Entry.Name,    
    AccessionN = Primary_Acc,  
    NombreProt = Protein.names,
    NombreGen = Gene.Names,    
    FuncionSLiM = ELMType,     
    SLiMID = ELMIdentifier,    
    # PosInicial = Start, 
    # PosFinal = End, 
    Secuencia = Sequence)

# Escribir los datos filtrados y modificados en un nuevo archivo CSV
write.table(x = elmHomo, 
            file = paste('C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/RepositorioGitHub/ELM/', 'TablaELM.csv', sep=''), 
            sep = ";", 
            append = FALSE, 
            quote = F, 
            na = "NA", 
            row.names = F, 
            col.names = T) 