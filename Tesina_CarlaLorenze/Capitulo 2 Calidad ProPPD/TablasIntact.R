library(dplyr)  # Carga la librería dplyr para manipulación de datos

# Directorios donde se encuentran los archivos
pocketDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/'
intactDir <- 'C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/files_intact/'
mappingDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/RepositorioGitHub/InteractoresIntAct/'
proteomic <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/InteractoresProteomica.csv'

# Carga de archivo de proteómica (Sanidas)
proteomicDF <- read.csv(file = proteomic, 
                        sep = ';', 
                        header = TRUE)

# Definición de archivos de interacción para cada proteína pocket
# Archivos para p107
directp107 <- '20230119_153404_p107_Direct.csv'
indirectp107 <- '20230119_153404_p107_Indirect.csv'
p107mappingIDs <- 'intact_p107.tsv'
p107ProPPD <- 'ProP-PD_p107_HD2.csv'

# Archivos para Rb
directRb <- '20230119_153404_Rb_Direct.csv'
indirectRb <- '20230119_153404_Rb_Indirect.csv'
RbmappingIDs <- 'intact_Rb.tsv'
RbProPPD <- 'ProP-PD_Rb_HD2.csv'

# Archivos para p130
directp130 <- '20230119_153404_p130_Direct.csv'
indirectp130 <- '20230119_153404_p130_Indirect.csv'
p130mappingIDs <- 'intact_p130.tsv'
p130ProPPD <- 'ProP-PD_p130_HD2.csv'

# Función para obtener numeros de acceso de interactores
getAccession <- function(directEvidence, indirectEvidence, IDinteractor){
  # Carga de datos de interacciones directas e indirectas
  pocket_D <- read.csv(file = paste(intactDir, directEvidence, sep=''), sep = ';', header = TRUE)
  pocket_I <- read.csv(file = paste(intactDir, indirectEvidence, sep=''), sep = ';', header = TRUE)
  
  # Identifica el interactor correcto dependiendo de la columna donde se encuentra
  pocket_D <- pocket_D %>%
    mutate(Interactor = if_else(X..ID.s..interactor.A != IDinteractor, X..ID.s..interactor.A, ID.s..interactor.B))
  pocket_I <- pocket_I %>%
    mutate(Interactor = if_else(X..ID.s..interactor.A != IDinteractor, X..ID.s..interactor.A, ID.s..interactor.B))
  
  # Limpieza de identificadores de UniProt y eliminaciones de sufijos no deseados
  pocket_D$Interactor <- gsub('uniprotkb:', '', pocket_D$Interactor)
  pocket_I$Interactor <- gsub('uniprotkb:', '', pocket_I$Interactor)
  pocket_I$Interactor <- gsub('-[12]', '', pocket_I$Interactor)
  
  # Obtención de identificadores únicos
  accession <- unique(c(pocket_D$Interactor, pocket_I$Interactor))
  tablaIntAct <- data.frame("AccessionN" = accession)
  
  # Eliminación de identificador específico no deseado
  tablaIntAct$AccessionN <- gsub('intact:EBI-971875', '', tablaIntAct$AccessionN)
  
  # Clasificación de interacciones según tipo de evidencia
  tablaIntAct <- tablaIntAct %>%
    mutate(TipoInteraccion = case_when(
      AccessionN %in% pocket_D$Interactor & AccessionN %in% pocket_I$Interactor ~ "Directa-Indirecta",
      AccessionN %in% pocket_D$Interactor ~ "Directa",
      AccessionN %in% pocket_I$Interactor ~ "Indirecta",
      TRUE ~ NA_character_
    ))
  
  cat(accession)  # Muestra los identificadores obtenidos
  return(tablaIntAct)
}

# Función para crear una tabla combinando datos de IntAct con información de mapeo y ProP-PD
creoTabla <- function(mappingFile, proppdFile, fileName, intAct){
  
  mappingDF <- read.csv(file = paste(mappingDir, 
                                     mappingFile, 
                                     sep=''), 
                        sep = '\t', 
                        header = TRUE)
  
  proPPD_DF <- read.csv(file = paste(pocketDir, 
                                     proppdFile, 
                                     sep=''), 
                        sep = ';', 
                        header = TRUE)
  
  # Agrega información de UniProt
  intAct <- intAct %>%
    left_join(mappingDF %>% select(Entry, 
                                   Entry.Name, 
                                   Gene.Names, 
                                   Protein.names), 
              by = c("AccessionN" = 'Entry'))
  
  # Identifica si los numeros de acceso están presentes en ProP-PD
  intAct <- intAct %>%
    mutate(HitrProPPD = AccessionN %in% proPPD_DF$Accession)
  
  # Identifica si los numeros de acceso están en datos de proteómica
  intAct <- intAct %>%
    mutate(Proteomica = AccessionN %in% proteomicDF$Accession)
  
  # Renombra columnas para mayor claridad
  intAct <- intAct %>%
    rename(UniprotID = Entry.Name, 
           NombreGen = Gene.Names, 
           NombreProt = Protein.names)
  
  # Limpia nombres de proteínas
  intAct$NombreProt <- gsub("\\s*\\([^\\)]+\\)", "", intAct$NombreProt)
  
  return(intAct)
}

# Función para guardar la tabla final en un archivo
guardoTabla <- function(intAct, fileName){
  
  pathOut <- paste(pocketDir, 
                   'RepositorioGitHub/InteractoresIntAct/', 
                   fileName, 
                   sep='')
  
  write.table(x = intAct, 
              file = pathOut, 
              sep = ";", 
              append = FALSE, 
              quote = FALSE, 
              na = "NA", 
              row.names = FALSE, 
              col.names = TRUE)
}

# Obtención de numeros de acceso para cada proteína pocket
accessions_p107 <- getAccession(directp107, indirectp107, "uniprotkb:P28749")
accessions_Rb <- getAccession(directRb, indirectRb, "uniprotkb:P06400")
accessions_p130 <- getAccession(directp130, indirectp130, "uniprotkb:Q08999")

# Creación de tablas finales combinando IntAct, mapeo y ProP-PD
intactRb <- creoTabla(RbmappingIDs, RbProPPD, 'Rb_IntAct.csv', accessions_Rb)
intactRb <- intactRb[-14,]  # Se elimina esta fila específica porque es un interactor no identificado ("[int]")
intactp107 <- creoTabla(p107mappingIDs, p107ProPPD, 'p107_IntAct.csv', accessions_p107)
intactp130 <- creoTabla(p130mappingIDs, p130ProPPD, 'p130_IntAct.csv', accessions_p130)

# Guarda las tablas finales en archivos CSV
guardoRb <- guardoTabla(intAct = intactRb, 
                        fileName = 'Rb_IntAct.csv')
guardop107 <- guardoTabla(intAct = intactp107, 
                        fileName = 'p107_IntAct.csv')
guardop130 <- guardoTabla(intAct = intactp130, 
                          fileName = 'p130_IntAct.csv')
