# Carga de bibliotecas necesarias para la visualización y manipulación de datos
library(ggplot2)      # Para crear gráficos
library(svglite)      # Para guardar gráficos en formato SVG
library(readxl)       # Para leer archivos de Excel

# Definición de las rutas de los directorios y archivos de datos
foldXDir <- 'C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/experimentalTP/NuevoScript/'  
rbPath <- 'Rb/ELMfiles/foldxFiles/'
rbBenchmark <- 'RB_Benchmarking.xlsx'
p107Path <- 'p107/ELM/FoldXscan/'
p107Benchmark <- 'p107_Benchmarking.xlsx'
benchmarkPath <- 'C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/experimentalTP/'

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
                             size = 22, 
                             colour = 'black'),
  axis.title = element_text(size = 30),
  axis.ticks = element_line(linewidth = 1, 
                            colour = 'black'),
  panel.border = element_rect(linewidth = 2, 
                              colour = 'black'),
  panel.grid = element_blank(),
  plot.margin = margin(t = 20, 
                       r = 5, 
                       b = 10, 
                       l = 5, 
                       unit = "pt")
)

# Función para procesar datos de FoldX y benchmarking
mismoPepTPTN <- function(matrix, benchFile) {
  
  # Inicialización de un dataframe vacío para almacenar datos de FoldX
  tPdf <- data.frame("SeqID" = NA, "Sequence" = NA, "Matrix" = NA, "Start" = NA, 
                     "End" = NA, "FoldX" = NA, "Min" = NA, "Sub.Sequence" = NA, 
                     "FoldXperRes" = NA, "RelativeMin" = NA, "MoreThan5" = NA, 
                     "RegexMatch" = NA, "RegexPattern" = NA)
  tNdf <- tPdf  # Copia del dataframe vacío para datos de negativos
  
  # Carga del archivo de benchmarking
  benchPaste <- paste(benchmarkPath, benchFile, sep = "")
  benchamarkDF <- read_xlsx(path = benchPaste)
  
  # Carga de archivos FoldX
  filePath <- paste(foldXDir, matrix, sep = "")
  files <- list.files(path = filePath, pattern = "TP.tsv")
  
  # Procesamiento de archivos de True Positives
  for (f in files) {
    archivoFoldX <- paste(filePath, f, sep = '/')
    singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
    
    # Filtrado de datos que cumplen ciertas condiciones
    selectedTP <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$RegexMatch == 'True'), ]
    
    # Selección del mínimo valor de FoldX
    minimo <- min(selectedTP$FoldX)
    selectedTP <- selectedTP[which(selectedTP$FoldX == minimo), ]
    
    # Agregación de los datos seleccionados al dataframe
    tPdf <- rbind(tPdf, selectedTP)
  }
  
  # Eliminación de la primera fila vacía
  tPdf <- tPdf[-1,]
  
  # Procesamiento de archivos de True Negatives
  filesTN <- list.files(path = filePath, pattern = "TN.tsv")
  for (f in filesTN) {
    archivoFoldX <- paste(filePath, f, sep = '/')
    singleF <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
    
    # Unión de los datos de True Negatives con los True Positives
    for (k in 1:nrow(tPdf)) {
      indicesTP <- rownames(tPdf)
      singleF <- singleF[which(rownames(singleF) == indicesTP[k]),]
      tNdf <- rbind(tNdf, singleF)
    }
  }
  
  # Eliminación de la primera fila vacía
  tNdf <- tNdf[-1,]
  
  # Combinación de True Positives y True Negatives en un solo dataframe
  elmInstances <- rbind(tPdf, tNdf)
  elmInstances$Instance <- NA
  
  # Columna que define la clase basada en el benchmarking
  columnaClass <- 'p107_Class'
  
  # Asignación de clases a las instancias basadas en el archivo de benchmarking
  for (i in 1:nrow(elmInstances)) {
    protName <- elmInstances$SeqID[i]
    for (j in 1:nrow(benchamarkDF)) {
      if (benchamarkDF$Name[j] == protName) {
        elmInstances$Instance[i] <- benchamarkDF[[columnaClass]][j]
      }
    }
  }
  
  # Impresión del rango de valores de FoldX
  cat('Rango de valores de ')
  cat(matrix)
  cat('\t')
  cat(range(elmInstances$FoldX))
  cat('\n')
  
  return(elmInstances)  # Devuelve el dataframe de instancias
}

# Función para calcular métricas de recall y especificidad
recall_specific <- function(umbral, df) {
  
  umbralArray <- c(1:10)
  
  # Inicialización de un dataframe para almacenar resultados
  recall_esp <- data.frame('Umbral' = umbralArray, 'Recall' = 0, 
                           'Especificidad' = 0, 'TP_menorAumbral' = 0, 
                           'TN_mayorAumbral' = 0)
  
  # Cálculo de totales de True Positives y True Negatives
  recall_esp$TP_total <- length(df$SeqID[which(df$Instance == "True Positive")]) 
  recall_esp$TN_total <- length(df$SeqID[which(df$Instance == "True Negative")])
  
  # Cálculo de recall y especificidad para cada umbral
  for (i in 1:nrow(recall_esp)) {
    umbral <- recall_esp$Umbral[i]
    recall_esp$TP_menorAumbral[which(recall_esp$Umbral == umbral)] <- length(df$SeqID[which((df$Instance == "True Positive") & df$FoldX <= umbral)])
    recall_esp$Recall[which(recall_esp$Umbral == umbral)] <- (recall_esp$TP_menorAumbral[i] / recall_esp$TP_total[i]) * 100
    recall_esp$TN_mayorAumbral[which(recall_esp$Umbral == umbral)] <- length(df$SeqID[which(df$Instance == "True Negative" & df$FoldX >= umbral)])
    recall_esp$Especificidad[which(recall_esp$Umbral == umbral)] <- recall_esp$TN_mayorAumbral[i] / recall_esp$TN_total[i]
  }
  
  return(recall_esp)  # Devuelve el dataframe de métricas
}

# Función para crear un histograma de FoldX
histograma <- function(df, FileOutName, yLims, xLims, xBy, yBy, labelPosition) {
  
  class <- c('True Positive', 'True Negative')
  df$Instance <- factor(df$Instance, levels = class, labels = class, ordered = TRUE)
  
  # Definición de colores para las instancias
  colores <- c("True Positive" = 'darkgreen', "True Negative" = '#990000ff')
  colores2 <- c("True Positive" = 'lightgreen', "True Negative" = '#e06666ff')
  
  max_plot <- yLims[2] - 0.1 * yLims[2]  # Ajuste para el gráfico
  
  # Creación del gráfico de histograma
  foldxPlot <- ggplot(df, aes(x = FoldX, 
                              fill = Instance)) + 
    scale_y_continuous(expand = c(0, 0), 
                       limits = yLims, 
                       breaks = seq(yLims[1], 
                                    yLims[2], 
                                    by = yBy),
                       sec.axis = sec_axis(trans = ~./max_plot, 
                                           name = "Fr. acumulada", 
                                           breaks = seq(0, 1, by = 0.2), 
                                           labels = scales::percent)) +
    scale_x_continuous(limits = xLims, 
                       breaks = seq(xLims[1], 
                                    xLims[2], 
                                    by = xBy), 
                       expand = c(0, 0)) +
    geom_histogram(binwidth = 0.5, 
                   center = 0.25, 
                   colour = 'black', 
                   position = 'dodge') +
    labs(x = "Valor de FoldX", y = "Instancias (N)") +
    scale_fill_manual(name = 'Instance', 
                      values = colores2, 
                      breaks = names(colores2), 
                      drop = FALSE)
  
  # Adición de frecuencias acumuladas al gráfico
  frecAcum <- foldxPlot + 
    stat_bin(data = subset(df, 
                           Instance == 'True Positive'), 
             aes(y = cumsum(after_stat(count)) / sum(after_stat(count)) * max_plot, col = Instance), 
             geom = "step", 
             binwidth = 0.25, 
             center = 0.125, 
             linewidth = 0.7) +
    stat_bin(data = subset(df, 
                           Instance == 'True Negative'), 
             aes(y = cumsum(after_stat(count)) / sum(after_stat(count)) * max_plot, 
                 col = Instance), 
             geom = "step", 
             binwidth = 0.25, 
             center = 0.125, 
             linewidth = 0.7) +
    scale_color_manual(values = colores, 
                       guide = NULL)
  
  # Personalización del gráfico
  frecAcum <- frecAcum + tema_plots + theme(legend.position = 'none', 
                                            legend.title = element_blank(), 
                                            legend.text = element_text(size = 15))
  
  # Guardado del gráfico
  ggsave(filename = FileOutName, 
         plot = frecAcum, 
         device = "png", 
         dpi = 600, 
         width = 16, 
         height = 12, 
         units = "cm")
}

# Función para trazar el Recall y la Especificidad
plotRecall <- function(df, 
                       xintercept, 
                       fileOutName) {
  
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
                       breaks = seq(0, 100, by = 20), 
                       name = "Recall (%)",
                       sec.axis = sec_axis(trans = ~ . / 100, 
                                           name = "Especificidad", 
                                           breaks = seq(0, 1, by = 0.2))) +
    scale_x_continuous(limits = c(0, 10), 
                       breaks = seq(0, 10, by = 2), 
                       expand = c(0, 0)) + 
    tema_plots
  
  # Guardado del gráfico
  ggsave(filename = fileOutName, plot = plotRecall, device = "png", dpi = 600, width = 16, height = 12, units = "cm")
}



##    1GUX_9    ##

p107_1GUX9 <- mismoPepTPTN(matrix = paste(p107Path, 
                                          '1GUX_9', 
                                          sep=''),
                           benchFile = p107Benchmark)

recallp107_1GUX9 <- recall_specific(df = p107_1GUX9, 
                                    #largoDeMatriz = 9,
                                    umbral=5)

p107_recall1GUX9 <- plotRecall(df = recallp107_1GUX9,
                               xintercept = 5, 
                               fileOutName = 'p107_recall1GUX9.png')
  
hist1GUX9 <- histograma(df= p107_1GUX9,
                        FileOutName = 'p107TP&N_1GUX9.png',
                        yLims= c(0,10),
                        yBy= 2,
                        xLims = c(-2,16),
                        xBy= 2,
                        labelPosition = 10)
##    1GUX_8    ##
p107_1GUX8  <- mismoPepTPTN(matrix = paste(p107Path, '1GUX_8', sep=''), 
                            benchFile = p107Benchmark)

recallp107_1GUX8 <- recall_specific(df = p107_1GUX8,
                                   # largoDeMatriz = 8,
                                    umbral=5)

p107_recall1GUX8 <- plotRecall(df = recallp107_1GUX8,
                               xintercept = 5, 
                               fileOutName = 'p107_recall1GUX8.png')

hist1GUX8 <- histograma(df= p107_1GUX8,
                        FileOutName = 'p107TP&N_1GUX8.png',
                        yLims= c(0,10),
                        yBy= 2,
                        xLims = c(-2,16),
                        xBy= 2,
                        labelPosition = 10)

##    2R7G_10   ##
p107_2R7G10  <- mismoPepTPTN(matrix = paste(p107Path, '2R7G_10', sep=''), 
                             benchFile = p107Benchmark)

recallp107_2R7G10 <- recall_specific(df = p107_2R7G10,
                                     #largoDeMatriz = 10,
                                     umbral=3)

p107_recallPlot2R7G10 <- plotRecall(df = recallp107_2R7G10,
                               xintercept = 3, 
                               fileOutName = 'recallp107_2R7G10.png')

hist2R7G10 <- histograma(df= p107_2R7G10,
                         FileOutName = 'p107TP&N_2R7G10.png',
                         yLims= c(0,6),
                         yBy= 2,
                         xLims = c(-2,12),
                         xBy= 2,
                         labelPosition = 10)
##    2R7G_6   ##
p107_2R7G6  <- mismoPepTPTN(matrix = paste(p107Path, '2R7G_6', sep=''), 
                            benchFile = p107Benchmark)

recallp107_2R7G6 <- recall_specific(df = p107_2R7G6, 
                                   # largoDeMatriz = 6,
                                    umbral=3)

p107_recallPlot2R7G6 <- plotRecall(df = recallp107_2R7G6,
                                    xintercept = 3, 
                                    fileOutName = 'recallp107_2R7G6.png')

hist2R7G6 <- histograma(df= p107_2R7G6,
                        FileOutName = 'p107TP&N_2R7G6.png',
                        yLims= c(0,6),
                        yBy= 2,
                        xLims = c(-2,12),
                        xBy= 2,
                        labelPosition = 10)
##    2R7G_5   ##
p107_2R7G5  <- mismoPepTPTN(matrix = paste(p107Path, '2R7G_5', sep=''), 
                            benchFile = p107Benchmark)
recallp107_2R7G5 <- recall_specific(df = p107_2R7G5, 
                                    #largoDeMatriz = 8,
                                    umbral=3)

p107_recallPlot2R7G5 <- plotRecall(df = recallp107_2R7G5,
                                   xintercept = 3, 
                                   fileOutName = 'recallp107_2R7G5.png')

hist2R7G5 <- histograma(df= p107_2R7G5,
                        FileOutName = 'p107TP&N_2R7G5.png',
                        yLims= c(0,6),
                        yBy= 2,
                        xLims = c(-2,12),
                        xBy= 2,
                        labelPosition = 10)
##    1N4M_9   ##
p107_1N4M9  <- mismoPepTPTN(matrix = paste(p107Path, '1N4M_9', sep=''), 
                            benchFile = p107Benchmark)

recallp107_1N4M9 <- recall_specific(df = p107_1N4M9,
                                    #largoDeMatriz = 9,
                                    umbral=3)

p107_recallPlot1N4M9 <- plotRecall(df = recallp107_1N4M9 ,
                                   xintercept = 3, 
                                   fileOutName = 'recallp107_1N4M9.png')

hist1N4M9 <- histograma(df= p107_1N4M9,
                        FileOutName = 'p107TP&N_1N4M9.png',
                        yLims= c(0,6),
                        yBy= 2,
                        xLims = c(-2,20),
                        xBy= 2,
                        labelPosition = 10)
##    1N4M_6   ##
p107_1N4M6  <- mismoPepTPTN(matrix = paste(p107Path, '1N4M_6', sep=''), 
                            benchFile = p107Benchmark)

recallp107_1N4M6 <- recall_specific(df = p107_1N4M6, 
                                    #largoDeMatriz = 6,
                                    umbral=3)

p107_recallPlot1N4M6 <- plotRecall(df = recallp107_1N4M6 ,
                                   xintercept = 3, 
                                   fileOutName = 'recallp107_1N4M6.png')


hist1N4M6 <- histograma(df= p107_1N4M6,
                        FileOutName = 'p107TP&N_1N4M6.png',
                        yLims= c(0,6),
                        yBy= 2,
                        xLims = c(-2,20),
                        xBy= 2,
                        labelPosition = 10)
##    1N4M_5   ##
p107_1N4M5  <- mismoPepTPTN(matrix = paste(p107Path, '1N4M_5', sep=''), 
                            benchFile = p107Benchmark)
recallp107_1N4M5 <- recall_specific(df = p107_1N4M5,
                                   # largoDeMatriz = 5,
                                    umbral=3)

p107_recallPlot1N4M5 <- plotRecall(df = recallp107_1N4M5 ,
                                   xintercept = 3, 
                                   fileOutName = 'recallp107_1N4M5.png')

hist1N4M5 <- histograma(df= p107_1N4M5,
                        FileOutName = 'p107TP&N_1N4M5.png',
                        yLims= c(0,6),
                        yBy= 2,
                        xLims = c(-2,20),
                        xBy= 2,
                        labelPosition = 10)


##    Retinoblastoma   ##




  ##    1GUX_9    ##

pep1GUX_9RB <- mismoPepTPTN(matrix = paste(rbPath, '1GUX_9', sep=''),
                            benchFile = rbBenchmark)

recall1GUX9RB <- recall_specific(df =pep1GUX_9RB,
                                # largoDeMatriz = 9,
                                 umbral=5)

recallPLOT1GUX9RB <- plotRecall(df = recall1GUX9RB ,
                                   xintercept = 5, 
                                   fileOutName = 'recall1GUX9RB.png')

mismoPEPplot_1gux9RB <- histograma(df= pep1GUX_9RB,
                           FileOutName= 'TP&N_1GUX9.png',
                           yLims= c(0,10),
                           xLims = c(-2,16),
                           xBy= 2,
                           yBy= 2,
                           labelPosition = 10)
##    1GUX_8   ##

samePEP1GUX8 <- mismoPepTPTN(matrix = paste(rbPath, '1GUX_8', sep=''),
                             benchFile = rbBenchmark)

recall1GUX8 <- recall_specific(df =samePEP1GUX8,
                               #largoDeMatriz = 8,
                               umbral=5)

recallPLOT1GUX8rb <- plotRecall(df = recall1GUX8 ,
                                xintercept = 5, 
                                fileOutName = 'recall1GUX8.png')

mismoPEPplot_1gux8 <- histograma(df= samePEP1GUX8,
                                 FileOutName= 'TP&N_1GUX8.png',
                                 yLims= c(0,10),
                                 xLims = c(-2,16),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 10)


  ##    2R7G_10   ##

samePEP_2R7G10 <- mismoPepTPTN(matrix = paste(rbPath, '2R7G_10', sep=''),
                               benchFile = rbBenchmark)

recall2R7G10 <- recall_specific(df =samePEP_2R7G10, 
                                #largoDeMatriz = 10,
                                umbral=3)

recallPLOT2R7G10 <- plotRecall(df = recall2R7G10 ,
                                xintercept = 3, 
                                fileOutName = 'recall2R7G10.png')

mismoPEPplot_2R7G10 <- histograma(df= samePEP_2R7G10,
                                 FileOutName= 'TP&N_2R7G10.png',
                                 yLims= c(0,6),
                                 xLims = c(-2,20),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 5)

##    2R7G_6  ##

samePEP_2R7G6 <- mismoPepTPTN(paste(rbPath, '2R7G_6', sep=''),
                              benchFile = rbBenchmark)

recall2R7G6 <- recall_specific(df =samePEP_2R7G6, 
                               #largoDeMatriz = 6,
                               umbral=3)

recallPLOT2R7G6 <- plotRecall(df = recall2R7G6 ,
                               xintercept = 3, 
                               fileOutName = 'recall2R7G6.png')



mismoPEPplot_2R7G6 <- histograma(df= samePEP_2R7G6,
                                  FileOutName= 'TP&N_2R7G6.png',
                                 yLims= c(0,6),
                                 xLims = c(-2,20),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 5)

##    2R7G_5  ##

samePEP_2R7G5 <- mismoPepTPTN(matrix = paste(rbPath, '2R7G_5', sep=''),
                              benchFile = rbBenchmark)

recall2R7G5 <- recall_specific(df =samePEP_2R7G5, 
                               #largoDeMatriz = 5,
                               umbral=3)


recallPLOT2R7G5 <- plotRecall(df = recall2R7G5 ,
                              xintercept = 3, 
                              fileOutName = 'recall2R7G5.png')

mismoPEPplot_2R7G5 <- histograma(df= samePEP_2R7G5,
                                 FileOutName= 'TP&N_2R7G5.png',
                                 yLims= c(0,6),
                                 xLims = c(-2,20),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 5)

  ##    1N4M_9   ##

samePEP_1n4m9 <- mismoPepTPTN(matrix = paste(rbPath, '1N4M_9', sep=''),
                              benchFile = rbBenchmark)

recall1n4m9 <- recall_specific(df =samePEP_1n4m9,  
                              # largoDeMatriz = 9,
                               umbral=3)


recallPLOT1n4m9 <- plotRecall(df = recall1n4m9 ,
                              xintercept = 3, 
                              fileOutName = 'RBrecall1n4m9.png')

mismoPEPplot_1n4m9 <- histograma(df= samePEP_1n4m9,
                                 FileOutName= 'TP&N_1N4M9.png',
                                 yLims= c(0,6),
                                 xLims = c(-2,20),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 5)


##    1N4M_6   ##

samePEP_1n4m6 <- mismoPepTPTN(matrix = paste(rbPath, '1N4M_6', sep=''),
                              benchFile = rbBenchmark)

recall1n4m6 <- recall_specific(df =samePEP_1n4m6, 
                               #largoDeMatriz = 6,
                               umbral=3)


recallPLOT1n4m6 <- plotRecall(df = recall1n4m6 ,
                              xintercept = 3, 
                              fileOutName = 'RBrecall1n4m6.png')

mismoPEPplot_1n4m6 <- histograma(df= samePEP_1n4m6,
                                 FileOutName= 'TP&N_1N4M6.png',
                                 yLims= c(0,6),
                                 xLims = c(-2,20),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 5)

##    1N4M_5   ##

samePEP_1n4m5 <- mismoPepTPTN(matrix = paste(rbPath, '1N4M_5', sep=''),
                              benchFile = rbBenchmark)

recall1n4m5<- recall_specific(df =samePEP_1n4m5, 
                             # largoDeMatriz = 5,
                              umbral=3)



recallPLOT1n4m5 <- plotRecall(df = recall1n4m5 ,
                              xintercept = 3, 
                              fileOutName = 'RBrecall1n4m5.png')

mismoPEPplot_1n4m5 <- histograma(df= samePEP_1n4m5,
                                 FileOutName= 'TP&N_1N4M5.png',
                                 yLims= c(0,6),
                                 xLims = c(-2,20),
                                 xBy= 2,
                                 yBy= 2,
                                 labelPosition = 5)



