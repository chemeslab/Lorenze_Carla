# Carga de bibliotecas necesarias para la visualización y manipulación de datos
library(ggplot2)      # Para crear gráficos
library(svglite)      # Para guardar gráficos en formato SVG
library(readxl)       # Para leer archivos de Excel

# Definición de las rutas de los directorios y archivos de datos
foldXDir <- 'C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/experimentalTP/'  
rbLxcxePath <- 'NuevoScript/Rb/LxCxEhits/foldxFiles/'
rbE2FPath <- 'NuevoScript/Rb/E2Fhits/foldxFiles/'
rbBenchFile <- 'RB_Benchmarking.xlsx'
p107LxcxePath <- 'NuevoScript/p107/LxCxEhits/FoldXscan'
p107E2FPath <- 'NuevoScript/p107/E2Fhits/FoldXScan'
p107BenchFile <- 'p107_Benchmarking.xlsx'

# Configuración del tema para los gráficos
tema_plots <- theme_bw() + theme(
  axis.text.x = element_text(angle = 0, 
                             hjust = 0.5, 
                             vjust = 0.5, 
                             size = 22, 
                             colour = 'black'),
  axis.text.y = element_text(angle = 0, 
                             hjust = 0.5, 
                             vjust = 0.5, 
                             size = 22, colour = 'black'),
  axis.title = element_text(size = 30),
  axis.ticks = element_line(linewidth = 1, 
                            colour = 'black'),
  panel.border = element_rect(linewidth = 2, 
                              colour = 'black'),
  panel.grid = element_blank(),
  plot.margin = margin(t = 20, 
                       r = 10, 
                       b = 10, 
                       l = 10, 
                       unit = "pt")
)

# Función para procesar datos de FoldX y benchmarking
bindersData <- function(matrix, BenchFile) {
  
  # Inicialización de un dataframe vacío para almacenar datos de FoldX
  FoldXDF <- data.frame("SeqID" = NA, "Sequence" = NA, "Matrix" = NA, "Start" = NA, 
                        "End" = NA, "FoldX" = NA, "Min" = NA, "Sub.Sequence" = NA, 
                        "FoldXperRes" = NA, "RelativeMin" = NA, "MoreThan5" = NA, 
                        "RegexMatch" = NA, "RegexPattern" = NA)
  
  # Carga del archivo de benchmarking
  benchPath <- paste(foldXDir, BenchFile, sep = '/')
  benchamarkDF <- read_xlsx(path = benchPath)
  
  # Carga de archivos FoldX
  filePath <- paste(foldXDir, matrix, sep = "/")
  files <- list.files(path = filePath)
  
  # Procesamiento de archivos de FoldX
  for (f in files) {
    if (startsWith(f, '2024')) {  # Procesa solo archivos de 2024
      archivoFoldX <- paste(filePath, f, sep = '/')
      singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
      
      # Filtrado de datos que cumplen ciertas condiciones
      df <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$RegexMatch == 'True'), ]
      
      # Selección del mínimo valor de FoldX
      minimo <- min(df$FoldX)
      df <- df[which(df$FoldX == minimo), ]
      
      # Agregación de los datos seleccionados al dataframe
      FoldXDF <- rbind(FoldXDF, df)
    }
  }
  
  # Eliminación de la primera fila vacía
  FoldXDF <- FoldXDF[-1,]
  
  # Impresión del rango de valores de FoldX
  cat('Rango de valores de matriz ')
  cat(matrix)
  cat('\t')
  cat(range(FoldXDF$FoldX))
  
  # Asignación de fortalezas a las interacciones
  benchamarkDF$Strength <- NA
  benchamarkDF$Strength[which(benchamarkDF$PD_COUP_SEC_interaction == "N" | benchamarkDF$AS_interaction == "N")] <- "N"
  benchamarkDF$Strength[which(benchamarkDF$PD_COUP_SEC_interaction == "W" | benchamarkDF$AS_interaction == "W1" | benchamarkDF$AS_interaction == "W2")] <- "W"
  benchamarkDF$Strength[which(benchamarkDF$PD_COUP_SEC_interaction == "SP" | benchamarkDF$AS_interaction == "SP" | benchamarkDF$AS_interaction == "P")] <- "SP"
  
  # Asignación de fortalezas en FoldXDF
  FoldXDF$Strength <- NA
  for (i in 1:nrow(FoldXDF)) {
    protName <- FoldXDF$SeqID[i]
    FoldXDF$Strength[i] <- benchamarkDF$Strength[which(benchamarkDF$Name == protName)]
  }
  
  return(FoldXDF)  # Devuelve el dataframe de FoldX
}

# Función para calcular métricas de recall y especificidad
recall_specific <- function(umbral, df) {
  umbralArray <- c(1:10)  # Definición de umbrales para evaluación
  
  # Inicialización de un dataframe para almacenar resultados
  recall_esp <- data.frame('Umbral' = umbralArray, 'Recall' = 0, 
                           'Especificidad' = 0, 'SPW_menorAumbral' = 0, 
                           'N_mayorAumbral' = 0)
  
  # Cálculo de totales de SPW y N
  recall_esp$SPW_total <- length(df$SeqID[which(df$Strength == "SP" | df$Strength == "W")]) 
  recall_esp$N_total <- length(df$SeqID[which(df$Strength == "N")])
  
  # Cálculo de recall y especificidad para cada umbral
  for (i in 1:nrow(recall_esp)) {
    umbral <- recall_esp$Umbral[i]
    
    recall_esp$SPW_menorAumbral[which(recall_esp$Umbral == umbral)] <- length(df$SeqID[which((df$Strength == "SP" | df$Strength == "W") & df$FoldX <= umbral)])
    recall_esp$Recall[which(recall_esp$Umbral == umbral)] <- (recall_esp$SPW_menorAumbral[i] / recall_esp$SPW_total[i]) * 100
    recall_esp$N_mayorAumbral[which(recall_esp$Umbral == umbral)] <- length(df$SeqID[which(df$Strength == "N" & df$FoldX >= umbral)])
    recall_esp$Especificidad[which(recall_esp$Umbral == umbral)] <- recall_esp$N_mayorAumbral[i] / recall_esp$N_total[i]
  }
  
  return(recall_esp)  # Devuelve el dataframe de métricas
}

# Función para trazar el Recall y la Especificidad
plotRecall <- function(df, xintercept, fileOutName) {
  
  # Creación del gráfico de Recall
  plotRecall <- ggplot(df, 
                       aes(x = Umbral, 
                           y = Recall)) +
    geom_vline(xintercept = xintercept, 
               linetype = 'dashed', 
               color = '#7393B3', 
               linewidth = 1) +
    geom_line(aes(y = Recall), 
              color = 'lightblue') +
    geom_point(aes(y = Recall), 
               color = "black", 
               pch = 21, 
               size = 3, 
               fill = 'darkblue') +
    geom_line(aes(y = Especificidad * 100), 
              color = 'pink') +
    geom_point(aes(y = Especificidad * 100), 
               color = "black", 
               pch = 24, 
               size = 3, 
               fill = 'darkred') +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, 110), 
                       breaks = seq(0, 100,
                                    by = 20), 
                       name = "Recall (%)", 
                       sec.axis = sec_axis(trans = ~ . / 100, name = "Especificidad", 
                                           breaks = seq(0, 1, by = 0.2))) +
    scale_x_continuous(limits = c(0, 10), 
                       breaks = seq(0, 10, by = 2), 
                       expand = c(0, 0)) + 
    tema_plots
  
  # Guardado del gráfico
  ggsave(filename = fileOutName, plot = plotRecall, device = "png", dpi = 600, width = 16, height = 12, units = "cm")
}

# Función para crear un histograma de FoldX
histograma <- function(df, FileOutName, yLims, yBy, labelPosition) {
  
  binder <- c('SP', 'W', 'N')  # Clases para la variable Strength
  df$Strength <- factor(df$Strength, levels = binder, labels = binder, ordered = TRUE)
  
  # Definición de colores para las instancias
  colores <- c("SP" = "darkgreen", "W" = 'lightgreen', "N" = 'pink')
  
  max_plot <- yLims[2] - 0.1 * yLims[2]  # Ajuste para el gráfico
  
  # Creación del gráfico de histograma
  foldxPlot <- ggplot(df, 
                      aes(x = FoldX, 
                          fill = Strength)) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = yLims, 
                       breaks = seq(yLims[1], 
                                    yLims[2], 
                                    by = yBy), 
                       sec.axis = sec_axis(trans = ~./max_plot, 
                                           name = "Fr. acumulada", 
                                           breaks = seq(0, 1, by = 0.2), 
                                           labels = scales::percent)) +
    scale_x_continuous(limits = c(-3, 12), 
                       breaks = seq(-2, 12, by = 2), 
                       expand = c(0, 0)) +
    geom_histogram(binwidth = 0.5, 
                   centre = 0.25, 
                   colour = 'black') +
    labs(x = "Valor de FoldX", y = "Hits (N)") +
    scale_fill_manual(name = "Strength", 
                      values = colores, 
                      breaks = names(colores), 
                      drop = FALSE)
  
  # Adición de frecuencias acumuladas al gráfico
  frecAcum <- foldxPlot + 
    stat_bin(aes(y = cumsum(after_stat(count)) / sum(after_stat(count)) * max_plot, fill = NULL), 
             geom = "step", 
             binwidth = 0.5, 
             center = 0.25, 
             linewidth = 0.7)
  
  # Personalización del gráfico
  frecAcum <- frecAcum + 
    tema_plots + 
    theme(legend.position = 'none', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 15))
  
  # Guardado del gráfico
  ggsave(filename = FileOutName, 
         plot = frecAcum, 
         device = "png", 
         dpi = 600, 
         width = 16, 
         height = 10, 
         units = "cm")
}

#Agrupo en estas "funciones" los programas para graficar los datos de Rb y p107 de todas las matrices y que el codigo no sea tan largo. NO SON FUNCIONES. 
datos_plotsP107 <- function(nada){
  
  ##    1GUX_9    ##
  p107_1GUX9 <- bindersData(matrix = paste(p107LxcxePath, '1GUX_9', sep='/'),
                         BenchFile = p107BenchFile)
  
  recallp107_1GUX9 <-  recall_specific(df = p107_1GUX9,#largoDeMatriz = 8,
                                       umbral=5)
  
  p107PLOTrecall_1GUX9 <- plotRecall(df = recallp107_1GUX9,
                                     xintercept = 5,
                                     fileOutName = 'p107PLOTrecall_1GUX9.png')
  
  hist1GUX9 <- histograma(df= p107_1GUX9,
                          FileOutName = 'p107Binders_1GUX9.png',
                          yLims= c(0,10),
                          yBy= 2,
                          labelPosition= 8)
  ##    1GUX_8    ##
  p107_1GUX8 <- bindersData(matrix = paste(p107LxcxePath, '1GUX_8', sep='/'),
                            BenchFile = p107BenchFile)
  
  recallp107_1GUX8 <- recall_specific(df = p107_1GUX8,#largoDeMatriz = 8,
                                      umbral=5)
  
  p107PLOTrecall_1GUX8<- plotRecall(df = recallp107_1GUX8,
                                     xintercept = 5,
                                     fileOutName = 'p107PLOTrecall_1GUX8.png')
  
  hist1GUX8 <- histograma(df= p107_1GUX8,
                         FileOutName = 'p107Binders_1GUX8.png',
                         yLims= c(0,10),
                         yBy= 2,
                         labelPosition= 8)
  ##    2R7G_10   ##
  p107_2R7G10 <- bindersData(matrix = paste(p107E2FPath, '2R7G_10', sep='/'),
                             BenchFile = p107BenchFile)
  
  recallp107_2R7G10 <- recall_specific(df = p107_2R7G10,#largoDeMatriz = 8,
                                      umbral=3)
  
  p107PLOTrecall2R7G10<- plotRecall(df = recallp107_2R7G10,
                                    xintercept = 3,
                                    fileOutName = 'recallp107_2R7G10.png')
  
  hist2R7G10 <- histograma(df= p107_2R7G10,
                           FileOutName = 'p107Binders_2R7G10.png',
                           yLims= c(0,6),
                           yBy= 2,
                           labelPosition= 8)
  
  ##    2R7G_6   ##
  p107_2R7G6 <- bindersData(matrix = paste(p107E2FPath, '2R7G_6', sep='/'),
                            BenchFile = p107BenchFile)
  
  recallp107_2R7G6 <- recall_specific(df = p107_2R7G6,#largoDeMatriz = 8,
                                       umbral=3)
  
  p107PLOTrecall2R7G6<- plotRecall(df = recallp107_2R7G6,
                                    xintercept = 3,
                                    fileOutName = 'recallp107_2R7G6.png')
  
  hist2R7G6 <- histograma(df= p107_2R7G6,
                          FileOutName = 'p107Binders_2R7G6.png',
                          yLims= c(0,6),
                          yBy= 2,
                          labelPosition= 8)
  
  ##    2R7G_5   ##
  
  p107_2R7G5 <- bindersData(matrix = paste(p107E2FPath, '2R7G_5', sep='/'),
                            BenchFile = p107BenchFile)
  
  recallp107_2R7G5 <- recall_specific(df = p107_2R7G5,#largoDeMatriz = 8,
                                      umbral=3)
  
  p107PLOTrecall2R7G5<- plotRecall(df = recallp107_2R7G5,
                                   xintercept = 3,
                                   fileOutName = 'recallp107_2R7G5.png')
  
  
  hist2R7G5 <- histograma(df= p107_2R7G5,
                          FileOutName = 'p107Binders_2R7G5.png',
                          yLims= c(0,6),
                          yBy= 2,
                          labelPosition= 8)
  
  ##    1N4M_9   ##
  
  p107_1N4M9 <- bindersData(matrix = paste(p107E2FPath, '1N4M_9', sep='/'),
                            BenchFile = p107BenchFile)
  
  recallp107_1N4M9 <- recall_specific(df = p107_1N4M9,#largoDeMatriz = 8,
                                      umbral=3)
  
  p107PLOTrecall1N4M9<- plotRecall(df = recallp107_1N4M9,
                                   xintercept = 3,
                                   fileOutName = 'recallp107_1N4M9.png')
  
  hist1N4M9 <- histograma(df= p107_1N4M9,
                          FileOutName = 'p107Binders_1N4M9.png',
                          yLims= c(0,6),
                          yBy= 2,
                          labelPosition= 8)
  
  ##    1N4M_6   ##
  
  p107_1N4M6 <- bindersData(matrix = paste(p107E2FPath, '1N4M_6', sep='/'),
                            BenchFile = p107BenchFile)
  
  recallp107_1N4M6 <- recall_specific(df = p107_1N4M6,#largoDeMatriz = 8,
                                      umbral=3)
  
  p107PLOTrecall1N4M6<- plotRecall(df = recallp107_1N4M6,
                                   xintercept = 3,
                                   fileOutName = 'recallp107_1N4M6.png')
  
  hist1N4M6 <- histograma(df= p107_1N4M6,
                          FileOutName = 'p107Binders_1N4M6.png',
                          yLims= c(0,6),
                          yBy= 2,
                          labelPosition= 8)
  
  ##    1N4M_5   ##
  p107_1N4M5 <- bindersData(matrix = paste(p107E2FPath, '1N4M_5', sep='/'),
                            BenchFile = p107BenchFile)
  
  recallp107_1N4M5<- recall_specific(df = p107_1N4M5,#largoDeMatriz = 8,
                                      umbral=3)
  
  p107PLOTrecall1N4M5<- plotRecall(df = recallp107_1N4M5,
                                   xintercept = 3,
                                   fileOutName = 'recallp107_1N4M5.png')
  
  hist1N4M5 <- histograma(df= p107_1N4M5,
                          FileOutName = 'p107Binders_1N4M5.png',
                          yLims= c(0,6),
                          yBy= 2,
                          labelPosition= 8)
  
  
}
datos_plotsRB <- function(nada){
  
  
  
      ## Retinoblastoma ##
  
  data_1gux9 <- bindersData(matrix = paste(rbLxcxePath, 
                                           '1GUX_9',
                                           sep=''),
                            BenchFile = rbBenchFile)
  
  recall1GUX9 <- recall_specific(df = data_1gux9, #largoDeMatriz = 8,
                                 umbral=5)
  
  
  rb_PLOTrecall_1GUX9<- plotRecall(df = recall1GUX9,
                                   xintercept = 5,
                                   fileOutName = 'rb_PLOTrecall_1GUX9.png')
  
  histograma_1gux9 <- histograma(df= data_1gux9,
                                 FileOutName = 'BindersRB_1GUX_9.png',
                                 yLims= c(0,10),
                                 yBy= 2,
                                 labelPosition= 8)
  
  
  data_1gux8 <- bindersData(matrix = paste(rbLxcxePath, '1GUX_8', sep=''),
                            BenchFile = rbBenchFile)
  
  recall1GUX8 <- recall_specific(df = data_1gux8, #largoDeMatriz = 8,
                                 umbral=5)
  
  rb_PLOTrecall_1GUX8 <- plotRecall(df = recall1GUX8,
                                   xintercept = 5,
                                   fileOutName = 'rb_PLOTrecall_1GUX8.png')
  
  histograma_1gux8 <- histograma(df= data_1gux8,
                                 FileOutName = 'BindersOut_1GUX_8.png',
                                 yLims= c(0,10),
                                 yBy= 2,
                                 labelPosition= 8)
  
  
  data_2r7g10 <- bindersData(matrix = paste(rbE2FPath, '2R7G_10', sep=''),
                            BenchFile =rbBenchFile )
  
  recall2r7g10  <- recall_specific(df = data_2r7g10,#largoDeMatriz = 8,
                                   umbral=3)
  
  rb_PLOTrecall_2r7g10 <- plotRecall(df = recall2r7g10,
                                    xintercept = 3,
                                    fileOutName = 'rb_PLOTrecall2r7g10.png')
  
  histograma_2r7g10 <- histograma(df= data_2r7g10,
                                 FileOutName = 'BindersOut_2R7G_10.png',
                                 yLims= c(0,6),
                                 yBy= 2,
                                 labelPosition= 6)
  
  
  
  data_2r7g6 <- bindersData(matrix = paste(rbE2FPath, '2R7G_6', sep=''),
                             BenchFile = rbBenchFile)
  
  recall2r7g6  <- recall_specific(df = data_2r7g6,#largoDeMatriz = 8,
                                  umbral=3)
  
  rb_PLOTrecall_2r7g6 <- plotRecall(df = recall2r7g6,
                                     xintercept = 3,
                                     fileOutName = 'rb_PLOTrecall2r7g6.png')
  
  histograma_2r7g6 <- histograma(df= data_2r7g6,
                                  FileOutName = 'BindersOut_2R7G_6.png',
                                 yLims= c(0,6),
                                 yBy= 2,
                                 labelPosition= 6)
  
  
  data_2r7g5 <- bindersData(matrix =  paste(rbE2FPath, '2R7G_5', sep=''),
                            BenchFile = rbBenchFile)
  
  recall2r7g5  <- recall_specific(df = data_2r7g5,#largoDeMatriz = 8,
                                  umbral=3)
  
  rb_PLOTrecall_2r7g5 <- plotRecall(df = recall2r7g5,
                                    xintercept = 3,
                                    fileOutName = 'rb_PLOTrecall2r7g5.png')
  
  histograma_2r7g5 <- histograma(df= data_2r7g5,
                                 FileOutName = 'BindersOut_2R7G_5.png',
                                 yLims= c(0,6),
                                 yBy= 2,
                                 labelPosition= 6)
  
  
  data_1n4m9 <- bindersData(matrix = paste(rbE2FPath, '1N4M_9', sep=''),
                            BenchFile = rbBenchFile)
  
  recall1n4m9  <- recall_specific(df = data_1n4m9,#largoDeMatriz = 8,
                                  umbral=3)
  
  rb_PLOTrecall_1n4m9 <- plotRecall(df = recall1n4m9,
                                    xintercept = 3,
                                    fileOutName = 'rb_PLOTrecall1n4m9.png')
  
  histograma_1n4m9 <- histograma(df= data_1n4m9,
                                 FileOutName = 'BindersOut_1N4M_9.png',
                                 yLims= c(0,6),
                                 yBy= 2,
                                 labelPosition= 6)
  
  data_1n4m6 <- bindersData(matrix = paste(rbE2FPath, '1N4M_6', sep=''),
                            BenchFile = rbBenchFile)
  
  recall1n4m6 <- recall_specific(df = data_1n4m6,#largoDeMatriz = 8,
                                 umbral=3)
  
  rb_PLOTrecall_1n4m6 <- plotRecall(df = recall1n4m6,
                                    xintercept = 3,
                                    fileOutName = 'rb_PLOTrecall1n4m6.png')
  
  histograma_1n4m6 <- histograma(df= data_1n4m6,
                                 FileOutName = 'BindersOut_1N4M_6.png',
                                 yLims= c(0,6),
                                 yBy= 2,
                                 labelPosition= 6)
  
  
  data_1n4m5 <- bindersData(matrix = paste(rbE2FPath, '1N4M_5', sep=''),
                            BenchFile = rbBenchFile)
  
  recall1n4m5 <- recall_specific(df = data_1n4m5,#largoDeMatriz = 8,
                                 umbral=3)
  
  
  rb_PLOTrecall_1n4m5 <- plotRecall(df = recall1n4m5,
                                    xintercept = 3,
                                    fileOutName = 'rb_PLOTrecall1n4m5.png')
  
  histograma_1n4m5 <- histograma(df= data_1n4m5,
                                 FileOutName = 'BindersOut_1N4M_5.png',
                                 yLims= c(0,6),
                                 yBy= 2,
                                 labelPosition= 6)
  
}

#Funcion donde guardo las listas de peptidos con valores de FoldX de los hits testeados 
miniDFbinders <- function(df,
                          fileName){
    finalDF <- df[,c("SeqID", 
                     "Sequence", 
                     "FoldX", 
                     "Sub.Sequence",
                     "Strength" )]
    
    write.table(x = finalDF,
                file = paste("C:/Users/Carla/Desktop/", fileName, sep='' ),
                sep = ";",
                append = FALSE,
                quote = F,
                na = "NA",
                row.names = F,
                col.names = T)
    
    return(finalDF)
}

rbLXCXE <- miniDFbinders(df = data_1gux8,
                         fileName = 'BindersRB_LXCXE.csv')

rbE2F <- miniDFbinders(data_1n4m5,
                       fileName = 'BindersRB_E2F.csv')

p107LXCXE <- miniDFbinders(p107_1GUX8,
                           fileName = 'Bindersp107_LXCXE.csv')

p107E2F <- miniDFbinders(p107_1N4M5,
                         fileName = 'Bindersp107_E2F.csv')

#Para chequear algunos datos
bindersE2F <-  benchamarkDF[which((benchamarkDF$Strength == "SP" | benchamarkDF$Strength == "W" |benchamarkDF$Strength == "N") & benchamarkDF$Rb_Class == "Hit" & benchamarkDF$Motif == 'E2F'),]

bindersLxcxe <-  benchamarkDF[which((benchamarkDF$Strength == "SP" | benchamarkDF$Strength == "W" |benchamarkDF$Strength == "N") & benchamarkDF$Rb_Class == "Hit" & (benchamarkDF$Motif == 'L.[CA].E'| benchamarkDF$Motif == 'L.S.E')),]


noMotifBinders <- benchamarkDF[which((benchamarkDF$Strength == "SP" | benchamarkDF$Strength == "W" |benchamarkDF$Strength == "N") & benchamarkDF$Rb_Class == "Hit" & benchamarkDF$Motif == 'No motif'),]