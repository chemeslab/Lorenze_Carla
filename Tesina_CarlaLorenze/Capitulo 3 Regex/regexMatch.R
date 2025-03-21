# Cargar las librerías necesarias
library(ggplot2)
library(ggh4x)

# Definir rutas de archivos
regexDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/Regex/NuevoRegex/Predictions.tsv'
foldXDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/FoldX/'

# Rutas para archivos de Rb
rb1GUX8 <- 'Retinoblastoma/HD2/FoldXScan/LxCxE/1GUX_8'
rb1N4M5 <- 'Retinoblastoma/HD2/FoldXScan/E2F/1N4M_5'
regexRBfile <- 'Rb_Predictions.tsv'

# Rutas para archivos de p107
p1071GUX8 <- 'p107/LxCxE/FoldXScan/1GUX_8'
p1071N4M5 <- 'p107/E2F/FoldXScan/1N4M_5'
regexp107file <- 'p107_Predictions.tsv'

# Definir el tema para los gráficos
tema_plots <- theme_bw() +  
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = 18, margin = margin(t = 0.7, b = 0.7, unit = "line"), colour = 'black'),
    axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 18, colour = 'black'),
    axis.title = element_text(size = 25),
    axis.ticks = element_line(linewidth = 1, colour = 'black'),
    axis.text.y.right = element_text(margin = margin(t = 0.7, r = 0.7, l = 0.7, unit = "line")),
    axis.text.y.left = element_text(margin = margin(r = 0.7, l = 0.7, unit = "line")),
    panel.border = element_rect(linewidth = 2, colour = 'black'),
    panel.grid = element_blank(),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt")
  )

# Función para generar gráficos de torta de pep con SLiM LxCxE
lxcxe_pieChart <- function(regexFile, pieName) {
  # Leer archivo de predicciones
  path <- paste(regexDir, regexFile, sep = '/')
  regexDF <- read.csv(file = path, sep = '\t', header = TRUE)
  
  # Definir columnas lógicas
  logicColumns <- c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 'LxCxE_Acidic_Hidrofobic', 'LxCxE_WC1', 'LxCxE_WC2', 'E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Asignar categorías a los péptidos
  regexDF$Categoria <- NA
  regexDF$Categoria[which(regexDF$LxCxE)] <- 'LxCxE'
  regexDF$Categoria[which(regexDF$LxCxE_Acidic)] <- 'LxCxE_Acidic'
  regexDF$Categoria[which(regexDF$LxCxE_Hidrofobic)] <- 'LxCxE_Hidrofobic'
  regexDF$Categoria[which(regexDF$LxCxE_Acidic_Hidrofobic)] <- 'LxCxE_Acidic_Hidrofobic'
  
  # Contar la frecuencia de cada categoría
  frecuencia_lxcxe <- data.frame(
    'Categoria' = c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 'LxCxE_Acidic_Hidrofobic'),
    'Frecuencia' = c(
      sum(regexDF$Categoria == 'LxCxE'),
      sum(regexDF$Categoria == 'LxCxE_Acidic'),
      sum(regexDF$Categoria == 'LxCxE_Hidrofobic'),
      sum(regexDF$Categoria == 'LxCxE_Acidic_Hidrofobic')
    )
  )
  
  # Definir colores
  colores <- c("LxCxE" = "#c83737", "LxCxE_Acidic" = '#672178', "LxCxE_Hidrofobic" = '#cd87de', "LxCxE_Acidic_Hidrofobic" = '#d38d5f')
  
  # Crear gráfico de torta
  lxcxeFeatures <- ggplot(frecuencia_lxcxe, aes(x = "", y = Frecuencia, fill = Categoria)) +
    geom_col(col = "black", linewidth = 0.5) +
    coord_polar(theta = "y") +
    geom_text(aes(x = 1.8, label = Frecuencia), position = position_stack(vjust = 0.5), size = 9) + 
    guides(fill = guide_legend(ncol = 2)) + 
    scale_fill_manual(values = colores) + 
    tema_plots +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.y = unit(1, 'line'), legend.text = element_text(size = 14), plot.margin = margin(b = 1, unit = "line"), text = element_text(size = 16))
  
  return(regexDF)
}

# Función para generar gráficos de torta de pep con SLiM E2F
e2f_pieChart <- function(regexFile, pieName) {
  # Leer archivo de predicciones
  path <- paste(regexDir, regexFile, sep = '/')
  regexDF <- read.csv(file = path, sep = '\t', header = TRUE)
  
  # Definir columnas lógicas
  logicColumns <- c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 'LxCxE_Acidic_Hidrofobic', 'LxCxE_WC1', 'LxCxE_WC2', 'E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Asignar categorías a los péptidos
  regexDF$Categoria <- NA
  regexDF$Categoria[which(regexDF$E2F)] <- 'E2F'
  regexDF$Categoria[which(regexDF$E2F_Acidic)] <- 'E2F_Acidic'
  regexDF$Categoria[which(regexDF$E2F_Aromatic)] <- 'E2F_Aromatic'
  regexDF$Categoria[which(regexDF$E2F_Acidic_Aromatic)] <- 'E2F_Acidic_Aromatic'
  
  # Contar la frecuencia de cada categoría
  frecuencia_e2f <- data.frame(
    'Categoria' = c('E2F', 'E2F_Acidic', 'E2F_Aromatic', 'E2F_Acidic_Aromatic'),
    'Frecuencia' = c(
      sum(regexDF$Categoria == 'E2F'),
      sum(regexDF$Categoria == 'E2F_Acidic'),
      sum(regexDF$Categoria == 'E2F_Aromatic'),
      sum(regexDF$Categoria == 'E2F_Acidic_Aromatic')
    )
  )
  
  # Definir colores
  colores <- c("E2F_Acidic" = '#672178', "E2F" = '#d38d5f', "E2F_Aromatic" = '#cd87de', "E2F_Acidic_Aromatic" = "#c83737")
  
  # Crear gráfico de torta
  e2fFeatures <- ggplot(frecuencia_e2f, aes(x = "", y = Frecuencia, fill = Categoria)) +
    geom_col(col = "black", linewidth = 0.5) +
    coord_polar(theta = "y") +
    geom_text(aes(x = 1.8, label = Frecuencia), position = position_stack(vjust = 0.5), size = 9) + 
    guides(fill = guide_legend(ncol = 2)) + 
    scale_fill_manual(values = colores) +
    theme_void() +
    theme(legend.position = "bottom", legend.title = element_blank(), legend.spacing.y = unit(1, 'line'), legend.text = element_text(size = 25), plot.margin = margin(b = 1, unit = "line"), text = element_text(size = 20))
  
  return(regexDF)
}

# Función para generar gráficos de puntos de dispersion de pep con SLiM LxCxE
lxcxeJitter <- function(matrixScanFile,
                        regexFile,
                        plotTitle,
                        jitterName){
  
  # Inicializa un data frame vacío para almacenar los datos procesados
  FoldXDF <- data.frame("SeqID"= NA,
                        "Sequence"= NA,
                        "Matrix"= NA,
                        "Start"= NA,
                        "End"= NA,
                        "FoldX"= NA,
                        "Min"= NA,
                        "Sub.Sequence"= NA,
                        "FoldXperRes"= NA,
                        "RelativeMin"= NA,
                        "MoreThan5"= NA,
                        "RegexMatch"= NA,
                        "RegexPattern"= NA)
  
  # Define la ruta de los archivos de FoldX
  filePath <- paste(foldXDir, matrixScanFile, sep="")
  files  <- list.files(path = filePath)
  
  # Carga el archivo de regex que contiene los patrones de SLiMs
  path <- paste(regexDir,regexFile, sep='/')
  regexDF <- read.csv(file = path, sep='\t', header=T)
  
  # Itera sobre los archivos de FoldX y filtra aquellos que empiezan con '2024'
  for (f in files){
    if (startsWith(f, '2024')){
      archivoFoldX <- paste(filePath, f, sep = '/')
      singleFile <- read.csv(archivoFoldX, header= T, sep='\t', skip = 4)
      
      # Filtra solo las entradas con mas de 5 residuos de largo y RegexMatch en TRUE
      df <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$RegexMatch =='True'), ]
      
      # Encuentra el valor mínimo de FoldX
      minimo <- min(df$FoldX)
      
      # Filtra las filas con el valor mínimo de FoldX
      df <- df[which(df$FoldX==minimo),]
      
      # Agrega los datos filtrados al data frame acumulado
      FoldXDF <- rbind(FoldXDF, df)
    }
  }
  
  # Elimina la primera fila de NA
  FoldXDF <- FoldXDF[-1,]
  
  # Agrega una nueva columna en regexDF para almacenar los valores de FoldX
  regexDF$FoldX_value <- 0
  
  # Asigna los valores de FoldX a regexDF utilizando las secuencias coincidentes
  for (i in 1:nrow(FoldXDF)){
    sequence <- FoldXDF$Sequence[i]
    regexDF$FoldX_value[which(regexDF$SEQ == sequence)] <- FoldXDF$FoldX[i]
  }
  
  # Convierte las columnas lógicas en valores booleanos
  logicColumns <- c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 'LxCxE_Acidic_Hidrofobic', 'LxCxE_WC1', 'LxCxE_WC2')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Categorización de los datos en base a la presencia de ciertos patrones
  regexDF$Categoria <- NA
  regexDF$Categoria[which(regexDF$LxCxE)] <- 'L.C.E'
  regexDF$Categoria[which(regexDF$LxCxE_Acidic)] <- '[DE]\nL.C.E'
  regexDF$Categoria[which(regexDF$LxCxE_Hidrofobic)] <- 'L.C.E.\n[WFILYVM]'
  regexDF$Categoria[which(regexDF$LxCxE_Acidic_Hidrofobic)] <- '[DE]\nL.C.E.\n[WFILYVM]'
  
  # Determina si pertenece a una categoría WC
  regexDF$WC <- FALSE
  regexDF$WC[which(regexDF$LxCxE_WC1 | regexDF$LxCxE_WC2)] <- TRUE
  
  # Filtra los hits con el motivo LxCxE
  lxcxeHits <- regexDF[which(regexDF$LxCxE), ]
  
  # Definición de colores para las categorías
  colores <- c("L.C.E"="#ff2a00ff" ,
               '[DE]\nL.C.E'= 'orange', 
               'L.C.E.\n[WFILYVM]'= 'skyblue', 
               '[DE]\nL.C.E.\n[WFILYVM]'= 'darkgreen')
  
  # Definición de etiquetas
  labels <- c('L.C.E', '[DE]\nL.C.E', 'L.C.E.\n[WFILYVM]', '[DE]\nL.C.E.\n[WFILYVM]')
  
  # Convierte la columna de categorías en factor ordenado
  lxcxeHits$Categoria <- factor(lxcxeHits$Categoria, 
                                labels = labels, 
                                ordered = T,
                                levels= labels) 
  
  # Generación del gráfico de dispersión con jitter
  #set.seed para que los puntos se ubiquen siempre en el mismo lugar y no se generen de manera aleatoria cada vez que recreo el grafico
  set.seed(3)
  regexFoldXPlot <- ggplot(data= lxcxeHits,
                           aes(x= Categoria, y= FoldX_value))+
    scale_fill_manual(values = colores, na.value = 'grey80', drop= F, guide = 'none')+
    scale_shape_manual(values = c('TRUE' = 24, 'FALSE'= 21), guide = 'none')+
    geom_jitter(aes(fill= Categoria, shape= WC), 
                col='black', width = 0.3, height = 0, size= 3, alpha= 0.7)+ 
    scale_y_continuous(name= 'Valor de FoldX', limits = c(-2,12), breaks = seq(-2,12, by= 2), expand= c(0, 0))+
    scale_x_discrete(labels=labels)+
    labs(title = plotTitle) + 
    tema_plots+ 
    theme(plot.margin = margin(b=0.7, t= 0.7, r= 0.7, l= 0.7, unit= "line"),
          plot.title = element_text(size= 25))+
    annotate(geom = 'text', x= 'L.C.E', y= 9, label = paste('N = ', nrow(lxcxeHits), sep=''), size= 8)+
    geom_hline(yintercept = 5, linetype= 'dashed', color= '#7393B3', linewidth = 1)
  
  # Guarda el gráfico como imagen PNG
  FileOut <- paste(regexDir, jitterName, sep='/')
  ggsave(filename= FileOut, plot =regexFoldXPlot, device = "png", dpi = 600, width = 20, height = 12, units = "cm")
  
  return(lxcxeHits)
}
# Función e2fJitter: Procesa archivos de matriz y genera un gráfico de dispersión de valores de FoldX
e2fJitter <- function(matrixScanFile, 
                      regexFile,
                      plotTitle, 
                      jitterName) {
  
  # Crear un data frame vacío para almacenar resultados
  FoldXDF <- data.frame("SeqID"= NA,
                        "Sequence"= NA,
                        "Matrix"= NA,
                        "Start"= NA,
                        "End"= NA,
                        "FoldX"= NA,
                        "Min"= NA,
                        "Sub.Sequence"= NA,
                        "FoldXperRes"= NA,
                        "RelativeMin"= NA,
                        "MoreThan5"= NA,
                        "RegexMatch"= NA,
                        "RegexPattern"= NA)
  
  # Construir la ruta del archivo y leer los archivos en el directorio
  filePath <- paste(foldXDir, matrixScanFile, sep="")
  files <- list.files(path = filePath)
  path <- paste(regexDir, regexFile, sep='/')
  regexDF <- read.csv(file = path, sep='\t', header=T)
  
  # Convertir columnas de lógica a tipo booleano
  logicColumns <- c('E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Clasificar categorías basadas en coincidencias de regex
  regexDF$Categoria <- NA
  regexDF$Categoria[which(regexDF$E2F)] <- 'General'
  regexDF$Categoria[which(regexDF$E2F_Acidic)] <- '[DE]P3'
  regexDF$Categoria[which(regexDF$E2F_Aromatic)] <- '[FY]P5'
  regexDF$Categoria[which(regexDF$E2F_Acidic_Aromatic)] <- '[DE]P3 &\n[FY]P5'
  
  # Procesar cada archivo en el directorio
  for (f in files) {
    if (startsWith(f, '2024')) {
      archivoFoldX <- paste(filePath, f, sep = '/')
      singleFile <- read.csv(archivoFoldX, header= T, sep='\t', skip = 4)
      
      # Filtrar filas donde MoreThan5 y RegexMatch son verdaderos
      df <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$RegexMatch =='True'), ]
      
      # Encontrar el valor mínimo de FoldX
      minimo <- min(df$FoldX)
      df <- df[which(df$FoldX == minimo),]
      
      # Agregar resultados al data frame principal
      FoldXDF <- rbind(FoldXDF, df)
    }
  }
  
  # Eliminar la primera fila vacía del data frame
  FoldXDF <- FoldXDF[-1,]
  
  # Inicializar columna para valores de FoldX
  regexDF$FoldX_value <- 0
  
  # Asignar valores de FoldX a las coincidencias en regexDF
  for (i in 1:nrow(FoldXDF)) {
    sequence <- FoldXDF$Sequence[i]
    regexDF$FoldX_value[which(regexDF$SEQ == sequence)] <- FoldXDF$FoldX[i]
  }
  
  # Filtrar coincidencias de E2F
  e2fHits <- regexDF[which(regexDF$E2F), ]
  
  # Definir colores para las categorías en el gráfico
  colores <- c("General"= 'orange', "[DE]P3"= 'skyblue', "[FY]P5"="#ff2a00ff", "[DE]P3 &\n[FY]P5"= 'darkgreen')
  levels <- c("General", "[DE]P3", "[FY]P5", "[DE]P3 &\n[FY]P5")
  
  # Convertir la categoría a un factor ordenado para el gráfico
  e2fHits$Categoria <- factor(e2fHits$Categoria, levels = levels, labels = levels, ordered = T)
  
  # Crear gráfico de dispersión (jitter)
  set.seed(3)
  regexFoldXPlot <- ggplot(data= e2fHits, 
                           aes(x= Categoria, 
                               y= FoldX_value)) +
    geom_jitter(aes(fill = Categoria), 
                width = 0.3, 
                height = 0, 
                size= 3, 
                alpha= 0.7, 
                pch=21) + 
    scale_y_continuous(name= 'Valor de FoldX', 
                       limits = c(-4,10), 
                       breaks = seq(-4,10, by= 2), 
                       expand= c(0, 0)) +
    scale_x_discrete(labels= levels) +
    scale_fill_manual(values = colores, 
                      guide = 'none') + 
    tema_plots + 
    labs(title = plotTitle) + 
    theme(plot.margin = margin(b = 0.7, 
                               t = 0.7, 
                               r = 0.7, 
                               l = 0.7, 
                               unit = "line"),
          plot.title = element_text(size = 25),
          legend.position = "none") +
    annotate(geom = 'text', 
             x= '[DE]P3 &\n[FY]P5', 
             y= 8, 
             label = paste('N = ',
                           nrow(e2fHits),
                           sep=''), 
             size=8) +
    geom_hline(yintercept = 3, 
               linetype= 'dashed', 
               color= '#7393B3', 
               linewidth = 1)
  
  # Guardar el gráfico generado como un archivo PNG
  FileOut <- paste(regexDir, 
                   jitterName, 
                   sep='/')
  
  ggsave(filename= FileOut, 
         plot = regexFoldXPlot, 
         device = "png", 
         dpi = 600, 
         width = 20, 
         height = 12, 
         units = "cm")
}

# Función lxcxeVenn: Realiza análisis de Venn sobre coincidencias de regex
lxcxeVenn <- function(regexFile) {
  
  # Leer el archivo de regex
  path <- paste(regexDir, regexFile, sep='/')
  regexDF <- read.csv(file = path, sep='\t', header=T)
  
  # Convertir columnas de lógica a tipo booleano
  logicColumns <- c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 'LxCxE_Acidic_Hidrofobic', 'LxCxE_WC1', 'LxCxE_WC2', 'E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  regexDF$WC <- FALSE
  regexDF$WC[which(regexDF$LxCxE_WC1 | regexDF$LxCxE_WC2)] <- TRUE
  
  # Contar coincidencias y mostrar resultados
  matchLxcxe <- length(regexDF$ID[which(regexDF$LxCxE)]) 
  cat('N total de péptidos con match L.[CAST].E: ', matchLxcxe, '\n')
  
  # Contar coincidencias específicas
  acidicos <- regexDF$ID[which(regexDF$LxCxE_Acidic)]
  cat('N total de péptidos con match acidico: ', length(acidicos), '\n')
  
  hidrofobicos <- regexDF$ID[which(regexDF$LxCxE_Hidrofobic)]
  cat('N total de péptidos con match hidrofobicos: ', length(hidrofobicos), '\n')
  
  wc <- regexDF$ID[which(regexDF$WC)]
  cat('N total de péptidos con wildcards: ', length(wc), '\n')
  
  # Análisis adicional de coincidencias
  acidicos1 <- acidicos[which(!(acidicos %in% hidrofobicos) & !(acidicos %in% wc))]
  cat('N de péptidos con match en pos ACIDICA pero NO hidrofobica NI WC: ', length(acidicos1), '\n')
  
  acidicos2 <- acidicos[which((acidicos %in% hidrofobicos) & !(acidicos %in% wc))]
  cat('N de péptidos con match en pos ACIDICA e hidrofobica pero SIN WC: ', length(acidicos2), '\n')
  
  acidicos3 <- acidicos[which((acidicos %in% wc) & !(acidicos %in% hidrofobicos))]
  cat('N de péptidos con match en pos ACIDICA, con WC pero SIN hidrofobica: ', length(acidicos3), '\n')
  
  all_acidicos <- acidicos[which((acidicos %in% wc) & (acidicos %in% hidrofobicos))]
  cat('N de péptidos con match en pos ACIDICA, HIDROFOBICA, con WC: ', length(all_acidicos), '\n')
  
  wc1 <- wc[which(!(wc %in% hidrofobicos) & !(wc %in% acidicos))]
  cat('N de péptidos con WC SIN pos hidrofobica NI acida: ', length(wc1), '\n')
  
  wc2 <- wc[which((wc %in% hidrofobicos) & !(wc %in% acidicos))]
  cat('N de péptidos con WC y pos hidrofobica pero SIN pos acida: ', length(wc2), '\n')
  
  hidrofobicos1 <- hidrofobicos[which(!(hidrofobicos %in% acidicos) & !(hidrofobicos %in% wc))]
  cat('N de péptidos con pos hidrofobica pero SIN pos acida NI WC: ', length(hidrofobicos1), '\n')
  
  return(regexDF)
}

# Función e2fVenn: Realiza análisis de Venn para coincidencias E2F
e2fVenn <- function(regexFile) {
  
  path <- paste(regexDir, regexFile, sep='/')
  regexDF <- read.csv(file = path, sep='\t', header=T)
  
  # Convertir columnas de lógica a tipo booleano
  logicColumns <- c('LxCxE', 'LxCxE_Acidic', 'LxCxE_Hidrofobic', 'LxCxE_Acidic_Hidrofobic', 'LxCxE_WC1', 'LxCxE_WC2', 'E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')
  regexDF[, logicColumns] <- lapply(regexDF[, logicColumns], as.logical)
  
  # Filtrar coincidencias de E2F
  e2fs <- regexDF[regexDF$E2F, c("ID", 'E2F', 'E2F_Acidic', 'E2F_Acidic_Aromatic', 'E2F_Aromatic')]
  
  # Contar coincidencias y mostrar resultados
  matchE2F <- regexDF$ID[which(regexDF$E2F)]
  cat('N total de péptidos con match E2F: ', length(matchE2F), '\n')
  
  acidicos <- regexDF$ID[which(regexDF$E2F_Acidic)]
  cat('N total de péptidos con match acidico: ', length(acidicos), '\n')
  
  aromaticos <- regexDF$ID[which(regexDF$E2F_Aromatic)]
  cat('N total de péptidos con match aromatico: ', length(aromaticos), '\n')
  
  acidicosAromaticos <- regexDF$ID[which(regexDF$E2F_Acidic_Aromatic)]
  cat('N total de péptidos con match acidico y aromatico: ', length(acidicosAromaticos), '\n')
  
  # Análisis adicional de coincidencias
  acidicos1 <- acidicos[which(!(acidicos %in% acidicosAromaticos))]
  cat('N total de péptidos con match acidico y NO aromatico: ', length(acidicos1), '\n')
  
  aromaticos1 <- aromaticos[which(!(aromaticos %in% acidicos))]
  cat('N total de péptidos con match aromatico y NO acidico: ', length(aromaticos1), '\n')
  
  e2fCore <- matchE2F[which(!(matchE2F %in% acidicos) & !(matchE2F %in% aromaticos))]
  cat('N total de péptidos con match E2F general, SIN aromatico NI acidico: ', length(e2fCore), '\n')
  
  return(regexDF)
}

rbLxcxe_pieChart <- lxcxe_pieChart(regexFile = regexRBfile,
                                   pieName= 'lxcxeMatch_Rb.png')

rbLxcxe_jitter <- lxcxeJitter(matrixScanFile = rb1GUX8,
                              regexFile = regexRBfile,
                              plotTitle= 'Hits Rb escaneados con 1GUX_8',
                             jitterName='Jitter_RbLxCxE.png')

rbVenn_lxcxe <- lxcxeVenn(regexFile = regexRBfile)

rbVenn_e2f <- e2fVenn(regexFile = regexRBfile )

rbE2F_pieChart <- e2f_pieChart(regexFile =regexRBfile,
                               pieName = 'E2Fmatch_Rb.png' )

rbe2f_jitter <- e2fJitter(matrixScanFile = rb1N4M5,
                          regexFile = regexRBfile ,
                          plotTitle= 'Hits Rb escaneados con 1N4M_5',
                          jitterName = 'regexFoldx_RbE2F.png')





p107Lxcxe_pieChart <- lxcxe_pieChart(regexFile = regexp107file, 
                                     pieName= 'lxcxeMatch_p107.png')

p107Lxcxe_jitter <- lxcxeJitter(matrixScanFile = p1071GUX8,
                                regexFile = regexp107file,
                                plotTitle= 'Hits p107 escaneados con 1GUX_8',
                                jitterName='regexFoldx_p107lxcxe.png')

p107Venn_lxcxe <- lxcxeVenn(regexFile = regexp107file)

p107Venn_e2f <- e2fVenn(regexFile = regexp107file)

p107E2F_pieChart <- e2f_pieChart(regexFile =regexp107file,pieName = 'E2Fmatch_p107.png' )

p107e2f_jitter <- e2fJitter(matrixScanFile = p1071N4M5,
                            regexFile = regexp107file,
                            plotTitle= 'Hits p107 escaneados con 1N4M_5',
                            jitterName = 'regexFoldx_p107E2F.png')
