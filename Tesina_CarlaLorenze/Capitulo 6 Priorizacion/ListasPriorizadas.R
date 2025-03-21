# Carga de bibliotecas necesarias para la manipulación de datos
library(readxl)      # Para leer archivos de Excel
library(dplyr)       # Para manipulación de datos

# Definición de las rutas de los directorios y archivos de datos
mainDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/'
regexDir <- 'Regex/NuevoRegex/Predictions.tsv/'
foldXDir <- 'FoldX/'
propdDir <- 'RepositorioGitHub/ProPPDHits/' 
intactDir <- 'RepositorioGitHub/InteractoresIntAct/'
benchmarkingDir <- 'RepositorioGitHub/FoldX/'
proteomicFile <- 'InteractoresProteomica.csv'

# Carga de datos de interacciones proteómicas
proteomicRb <- read.csv(file = paste(mainDir, proteomicFile, sep = ''), sep = ';', header = TRUE)

# Archivos de datos adicionales
hd2File <- 'LibraryHD2.csv'
rbPropdFile <- 'RbHitsProPPD.csv'
rbValidationFile <- 'RB_Benchmarking.xlsx'
rb1GUX8 <- 'Retinoblastoma/HD2/FoldXScan/LxCxE/1GUX_8'
rb1N4M5 <- 'Retinoblastoma/HD2/FoldXScan/E2F/1N4M_5'
regexRBfile <- 'Rb_Predictions.tsv'

p107PropdFile <- 'p107HitsProPPD.csv'
p107ValidationFile <- 'p107_Benchmarking.xlsx'
p1071GUX8 <- 'p107/LxCxE/FoldXScan/1GUX_8'
p1071N4M5 <- 'p107/E2F/FoldXScan/1N4M_5'
regexp107file <- 'p107_Predictions.tsv'

intactp107 <- 'p107_IntAct.csv'
intactRb <- 'Rb_IntAct.csv'

# Carga de datos de la biblioteca HD2
hd2Path <- paste(mainDir, hd2File, sep = '')
hd2DF <- read.csv(file = hd2Path, sep = ';', header = TRUE)

# Función para crear una tabla de resultados para peptidos que contienen coincidencias de regex
tablaLxCxE <- function(regexFile, propdFile, valFile, matrixScanFile) {
  
  # Inicialización de un dataframe vacío para almacenar datos de FoldX
  FoldXDF <- data.frame("SeqID" = NA, 
                        "Sequence" = NA, 
                        "Matrix" = NA, 
                        "Start" = NA, 
                        "End" = NA, 
                        "FoldX" = NA, 
                        "Min" = NA, 
                        "Sub.Sequence" = NA, 
                        "FoldXperRes" = NA, 
                        "RelativeMin" = NA, 
                        "MoreThan5" = NA, 
                        "RegexMatch" = NA, 
                        "RegexPattern" = NA)
  
  # Carga de datos de FoldX
  filePath <- paste(mainDir, foldXDir, matrixScanFile, sep = "")
  files <- list.files(path = filePath)
  
  # Carga de datos de regex
  regexPath <- paste(mainDir, regexDir, regexFile, sep = '')
  regexDF <- read.csv(file = regexPath, sep = '\t', header = TRUE)
  
  # Carga de datos de ProPPD
  propdPath <- paste(mainDir, propdDir, propdFile, sep = '')
  propdDF <- read.csv(file = propdPath, sep = ';', header = TRUE)
  
  # Carga de datos de benchmarking
  validationPath <- paste(mainDir, valFile, sep = '')
  benchDF <- read_xlsx(path = validationPath)
  
  # Procesamiento de los archivos de FoldX
  for (f in files) {
    archivoFoldX <- paste(filePath, f, sep = '/')
    singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
    
    # Filtrado de datos que cumplen las condiciones
    df <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$RegexMatch == 'True'), ]
    
    # Selección del mínimo valor de FoldX
    minimo <- min(df$FoldX)
    df <- df[which(df$FoldX == minimo), ]
    
    # Agregación de los datos seleccionados al dataframe
    FoldXDF <- rbind(FoldXDF, df)
  }
  
  # Eliminación de la primera fila vacía
  FoldXDF <- FoldXDF[-1,]
  
  # Definición de columnas lógicas para regex
  logicColumns <- c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 
                    'LxCxE_Acidic_Hidrofobic', 'LxCxE_WC1', 'LxCxE_WC2')
  
  # Conversión de columnas a valores lógicos
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Filtrado de hits de LxCxE
  lxcxeHits <- regexDF[which(regexDF$LxCxE), c("ID", "SEQ", "LxCxE", "LxCxE_Seq", 
                                               "LxCxE_Acidic", "LxCxE_Acidic_Seq", 
                                               "LxCxE_Hidrofobic", "LxCxE_Hidrofobic_Seq", 
                                               "LxCxE_Acidic_Hidrofobic", "LxCxE_Acidic_Hidrofobic_Seq", 
                                               "LxCxE_WC1", "LxCxE_WC1_Seq", 
                                               "LxCxE_WC2", "LxCxE_WC2_Seq")]
  
  # Asignación de categorías a los hits
  lxcxeHits$WC <- FALSE
  lxcxeHits$WC[which(lxcxeHits$LxCxE_WC1 | lxcxeHits$LxCxE_WC2)] <- TRUE
  lxcxeHits$Categoria <- NA
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE)] <- 'Core'
  lxcxeHits$Categoria[which(lxcxeHits$WC)] <- 'PV'
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE_Acidic)] <- 'Acid'
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE_Acidic & lxcxeHits$WC)] <- 'Acid_PV'
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE_Hidrofobic)] <- 'Hidrofobic'
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE_Acidic_Hidrofobic)] <- 'Acid_Hidrofobic'
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE_Hidrofobic & lxcxeHits$WC)] <- 'Hidrofobic_PV'
  lxcxeHits$Categoria[which(lxcxeHits$LxCxE_Acidic_Hidrofobic & lxcxeHits$WC)] <- 'Acid_PV_Hidrofobic'
  
  # Ordenación de categorías
  orderedCat <- c('Acid_PV_Hidrofobic', 'Hidrofobic_PV', 'Acid_Hidrofobic', 'Hidrofobic', 'Acid_PV', 'Acid', 'PV', 'Core')
  lxcxeHits$Categoria <- factor(x = lxcxeHits$Categoria, levels = orderedCat)
  
  # Unión de datos de FoldX
  lxcxeHits <- lxcxeHits %>%
    left_join(FoldXDF %>% select(SEQ = Sequence, FoldX), by = 'SEQ')
  
  # Asignación de fortalezas en datos de benchmarking
  benchDF$Strength <- NA
  benchDF$Strength[which(benchDF$PD_COUP_SEC_interaction == "N" | benchDF$AS_interaction == "N")] <- "N"
  benchDF$Strength[which(benchDF$PD_COUP_SEC_interaction == "W" | benchDF$AS_interaction == "W1" | benchDF$AS_interaction == "W2")] <- "W"
  benchDF$Strength[which(benchDF$PD_COUP_SEC_interaction == "SP" | benchDF$AS_interaction == "SP" | benchDF$AS_interaction == "P")] <- "SP"
  
  # Unión de datos de ProPPD y benchmarking
  lxcxeHits <- lxcxeHits %>%
    left_join(propdDF %>% select(SecC, AccessionN, RSA, IUPred, NResPfam, UniprotID, NombreProt, NombreGen, LocCelular, Proteomica, IntAct), 
              by = c("SEQ" = "SecC")) %>%
    left_join(benchDF %>% select(Peptide, Strength), by = c("SEQ" = "Peptide"))
  
  # Selección final de columnas
  lxcxeHits <- lxcxeHits %>%
    select(UniprotID, AccessionN, NombreProt, NombreGen, Secuencia = SEQ, Categoria, 
           RSA, IUPred, NResPfam, ValorFoldX = FoldX, ValidacionExp = Strength, 
           Proteomica, IntAct)
  
  return(lxcxeHits)  # Devuelve el dataframe procesado
}

# Función para crear una tabla de resultados para peptidos que contienen coincidencias de E2F
tablaE2F <- function(regexFile, propdFile, valFile, matrixScanFile) {
  
  # Inicialización de un dataframe vacío para almacenar datos de FoldX
  FoldXDF <- data.frame("SeqID" = NA, 
                        "Sequence" = NA, 
                        "Matrix" = NA, 
                        "Start" = NA, 
                        "End" = NA, 
                        "FoldX" = NA, 
                        "Min" = NA, 
                        "Sub.Sequence" = NA, 
                        "FoldXperRes" = NA, 
                        "RelativeMin" = NA, 
                        "MoreThan5" = NA, 
                        "RegexMatch" = NA, 
                        "RegexPattern" = NA)
  
  # Carga de datos de FoldX
  filePath <- paste(mainDir, foldXDir, matrixScanFile, sep = "")
  files <- list.files(path = filePath)
  
  # Carga de datos de regex
  regexPath <- paste(mainDir, regexDir, regexFile, sep = '')
  regexDF <- read.csv(file = regexPath, sep = '\t', header = TRUE)
  
  # Carga de datos de ProPPD
  propdPath <- paste(mainDir, propdDir, propdFile, sep = '')
  propdDF <- read.csv(file = propdPath, sep = ';', header = TRUE)
  
  # Carga de datos de benchmarking
  validationPath <- paste(mainDir, valFile, sep = '')
  benchDF <- read_xlsx(path = validationPath)
  
  # Procesamiento de los archivos de FoldX
  for (f in files) {
    archivoFoldX <- paste(filePath, f, sep = '/')
    singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
    
    # Filtrado de datos que cumplen las condiciones
    df <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$RegexMatch == 'True'), ]
    
    # Selección del mínimo valor de FoldX
    minimo <- min(df$FoldX)
    df <- df[which(df$FoldX == minimo), ]
    
    # Agregación de los datos seleccionados al dataframe
    FoldXDF <- rbind(FoldXDF, df)
  }
  
  # Eliminación de la primera fila vacía
  FoldXDF <- FoldXDF[-1,]
  
  # Definición de columnas lógicas para regex
  logicColumns <- c('E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Filtrado de hits de E2F
  e2fHits <- regexDF[which(regexDF$E2F), c("ID", "SEQ", "E2F", "E2F_Seq", 
                                           "E2F_Acidic", "E2F_Acidic_Seq", 
                                           "E2F_Acidic_Aromatic", "E2F_Acidic_Aromatic_Seq", 
                                           "E2F_Aromatic", "E2F_Aromatic_Seq")]
  
  # Asignación de categorías a los hits
  e2fHits$Categoria <- NA
  e2fHits$Categoria[which(e2fHits$E2F)] <- 'Core'
  e2fHits$Categoria[which(e2fHits$E2F_Acidic)] <- 'Acid'
  e2fHits$Categoria[which(e2fHits$E2F_Aromatic)] <- 'Aromatic'
  e2fHits$Categoria[which(e2fHits$E2F_Acidic_Aromatic)] <- 'Acid_Aromatic'
  
  # Ordenación de categorías
  orderedCat <- c('Acid_Aromatic', 'Aromatic', 'Acid', 'Core')
  e2fHits$Categoria <- factor(x = e2fHits$Categoria, levels = orderedCat)
  
  # Unión de datos de FoldX
  e2fHits <- e2fHits %>%
    left_join(FoldXDF %>% select(SEQ = Sequence, FoldX), by = 'SEQ')
  
  # Asignación de fortalezas en datos de benchmarking
  benchDF$Strength <- NA
  benchDF$Strength[which(benchDF$PD_COUP_SEC_interaction == "N" | benchDF$AS_interaction == "N")] <- "N"
  benchDF$Strength[which(benchDF$PD_COUP_SEC_interaction == "W" | benchDF$AS_interaction == "W1" | benchDF$AS_interaction == "W2")] <- "W"
  benchDF$Strength[which(benchDF$PD_COUP_SEC_interaction == "SP" | benchDF$AS_interaction == "SP" | benchDF$AS_interaction == "P")] <- "SP"
  
  # Unión de datos de ProPPD y benchmarking
  e2fHits <- e2fHits %>%
    left_join(propdDF %>% select(SEQ = SecC, AccessionN, NombreProt, 
                                 NombreGen, LocCelular, RSA, IUPred, 
                                 NResPfam, DomPfam, UniprotID, 
                                 Proteomica, IntAct), by = 'SEQ') %>%
    left_join(benchDF %>% select(SEQ = Peptide, Strength), by = 'SEQ')
  
  # Selección final de columnas
  e2fHits <- e2fHits %>%
    select(UniprotID, AccessionN, NombreProt, NombreGen, Secuencia = SEQ, 
           Categoria, RSA, IUPred, NResPfam, ValorFoldX = FoldX, 
           ValidacionExp = Strength, Proteomica, IntAct)
  
  return(e2fHits)  # Devuelve el dataframe procesado
}

# Función para encontrar péptidos que coinciden en ambas tablas
getBothMatch <- function(dfLxcxe, dfE2F, fileLxcxe, fileE2F) {
  
  dfLxcxe$BothMatch <- FALSE  # Inicialización de columna para coincidencias
  dfE2F$BothMatch <- FALSE     # Inicialización de columna para coincidencias
  
  # Comprobación de coincidencias en la tabla LxCxE
  for (i in 1:nrow(dfLxcxe)) {
    sequence <- dfLxcxe$Secuencia[i]
    dfE2F$BothMatch[which(dfE2F$Secuencia == sequence)] <- TRUE
  }
  
  # Comprobación de coincidencias en la tabla E2F
  for (i in 1:nrow(dfE2F)) {
    sequence <- dfE2F$Secuencia[i]
    dfLxcxe$BothMatch[which(dfLxcxe$Secuencia == sequence)] <- TRUE
  }
  
  # Guardado de resultados en archivos
  write.table(x = dfLxcxe, file = paste("C:/Users/Carla/Desktop/", fileLxcxe, sep = ''), 
              sep = ";", append = FALSE, quote = FALSE, na = "NA", 
              row.names = FALSE, col.names = TRUE)
  
  write.table(x = dfE2F, file = paste("C:/Users/Carla/Desktop/", fileE2F, sep = ''), 
              sep = ";", append = FALSE, quote = FALSE, na = "NA", 
              row.names = FALSE, col.names = TRUE)
}

tabla_rbLxcxe <- tablaLxCxE(regexFile = regexRBfile,
                      propdFile = rbPropdFile, 
                      valFile = rbValidationFile,
                      matrixScanFile = rb1GUX8)

tabla_p107Lxcxe <- tablaLxCxE(regexFile = regexp107file,
                        propdFile = p107PropdFile, 
                        valFile = p107ValidationFile,
                        matrixScanFile = p1071GUX8)

tabla_rbE2F <- tablaE2F(regexFile = regexRBfile,
                  propdFile = rbPropdFile, 
                  valFile = rbValidationFile,
                  matrixScanFile = rb1N4M5)

tabla_p107E2F <- tablaE2F(regexFile = regexp107file,
                    propdFile = p107PropdFile, 
                    valFile = p107ValidationFile,
                    matrixScanFile = p1071N4M5)


rb <- getBothMatch(dfLxcxe = tabla_rbLxcxe,
                   dfE2F = tabla_rbE2F, 
                   fileLxcxe = 'Rb_LxCxE.csv',
                   fileE2F = 'Rb_E2F.csv' )

p107 <- getBothMatch(dfLxcxe = tabla_p107Lxcxe, 
                     dfE2F = tabla_p107E2F, 
                     fileLxcxe = 'p107_LxCxE.csv',
                     fileE2F ='p107_E2F.csv')

