#Dotplots de TP para FoldX
# Carga de bibliotecas necesarias para la visualización y manipulación de datos
library(ggplot2)      # Para crear gráficos
library(patchwork)    # Para combinar múltiples gráficos
library(dplyr)        # Para manipulación de datos
library(purrr)        # Para programación funcional

# Configuración de la fuente que se usará en los gráficos
windowsFonts('Courier' = windowsFont('Courier New'))

# Definición de las rutas de los directorios que contienen los archivos de datos
FoldXDir <- 'C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/experimentalTP/NuevoScript/'
rbDir <- 'Rb/ELMfiles/foldxFiles/'
p107Dir <- 'p107/ELM/FoldXscan/'

# Personalización del tema para los gráficos
theme_plots_recall <- theme_bw() +
  theme(
    axis.text.x.bottom = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10, colour = 'black'),
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 15, family = 'Courier', colour = 'black'),
    axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.2, size = 10, colour = 'black'),
    axis.title = element_text(size = 15, vjust = -2, hjust = 0.5),
    plot.title = element_text(size = 20, hjust = 0.5),
    plot.margin = margin(t = 1, r = 1, l = 0.5, unit = "line"),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.spacing.y = unit(2, 'line'),
    legend.key.size = unit(1.5, 'line')
  )

# Función para crear gráficos de FoldX
plotTPELM <- function(matrix, fileOutName, Ngroup, ylims, yBy, plotTitle, ultimaRow) {
  
  # Construcción de la ruta del archivo
  filePath <- paste(FoldXDir, matrix, sep = "")
  
  # Listado de archivos FoldX que cumplen con el patrón especificado
  filesFoldX <- list.files(path = filePath, pattern = "TP.tsv")
  
  # Creación de un dataframe vacío para almacenar los datos de FoldX
  foldXDF <- data.frame("SeqID" = NA, "Sequence" = NA, "Matrix" = NA, "Start" = NA, "End" = NA, 
                        "FoldX" = NA, "Min" = NA, "Sub.Sequence" = NA, "FoldXperRes" = NA, 
                        "RelativeMin" = NA, "MoreThan5" = NA, "RegexMatch" = NA, "RegexPattern" = NA)
  
  # Lectura de los archivos FoldX y combinación de los datos en un solo dataframe
  for (f in filesFoldX) {
    if (startsWith(f, '202')) {
      archivoFoldX <- paste(filePath, f, sep = '/')
      singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
      foldXDF <- rbind(foldXDF, singleFile)
    }
  }
  
  # Eliminación de la primera fila vacía
  foldXDF <- foldXDF[-1,]
  
  # Cálculo del punto medio para cada fila
  foldXDF$MidPoint <- rowSums(foldXDF[, c(4:5)]) / 2
  seqUnicas <- unique(foldXDF$Sequence)
  
  # Inicialización de variables para el procesamiento de secuencias
  endTemp <- 0
  endTempPlot <- c()
  
  # Ajuste de puntos medios según las secuencias únicas
  for (i in 1:length(seqUnicas)) {
    seqString <- unlist(strsplit(seqUnicas[i], split = '', fixed = TRUE))
    if (i == 1) {
      endTemp <- length(seqString)
      endTempPlot <- append(x = endTempPlot, values = endTemp)
    } else {
      foldXDF$MidPoint[foldXDF$Sequence == seqUnicas[i]] <- foldXDF$MidPoint[foldXDF$Sequence == seqUnicas[i]] + endTemp
      endTemp <- endTemp + length(seqString)
      endTempPlot <- append(x = endTempPlot, values = endTemp)
    }
  }
  
  # Ajuste para el gráfico
  endTempPlot <- endTempPlot + 0.5
  
  # Conversión de columnas lógicas a valores booleanos
  foldXDF$RegexMatch <- as.logical(foldXDF$RegexMatch)
  foldXDF$MoreThan5 <- as.logical(foldXDF$MoreThan5)
  foldXDF$RelativeMin <- as.logical(foldXDF$RelativeMin)
  foldXDF$Min <- as.logical(foldXDF$Min)
  
  # Clasificación de datos en categorías según condiciones lógicas
  colMin <- 'Min'
  col_Min <- foldXDF[[colMin]]
  foldXDF$Category <- NA
  foldXDF$Category[which(!col_Min & !foldXDF$RegexMatch)] <- 'noMin_noRegex'
  foldXDF$Category[which(col_Min & !foldXDF$RegexMatch)] <- 'MinFoldx'
  foldXDF$Category[which(!col_Min & foldXDF$RegexMatch)] <- 'RegexMatch'
  foldXDF$Category[which(col_Min & foldXDF$RegexMatch)] <- 'MinFoldx_RegexMatch'
  foldXDF$Category[which(!foldXDF$MoreThan5)] <- 'LessThan5Res'
  
  # Creación de un dataframe para la secuencia
  secuencia <- paste(seqUnicas, collapse = '')
  secuencia <- unlist(strsplit(secuencia, split = '', fixed = TRUE))
  
  uniquePep <- unique(foldXDF$Sequence)
  grupos <- 0
  
  dfSeq <- data.frame(Secuencia = secuencia, Res = 1:length(secuencia))
  dfSeq$Sequence <- NA
  for (i in 1:length(seqUnicas)) {
    ss <- ((i - 1) * 18) + 1
    ee <- 18 * i
    dfSeq$Sequence[ss:ee] <- seqUnicas[i]
  }
  
  # Combinación de datos para crear un dataframe final
  df <- merge(x = foldXDF, y = dfSeq, by.x = c("MidPoint", "Sequence"), by.y = c("Res", "Sequence"), all = TRUE)
  
  # Asignación de grupos a las secuencias
  df$Group <- 0
  for (i in 1:length(uniquePep)) {
    df$Group[df$Sequence == uniquePep[i]] <- grupos
    if (i %% Ngroup == 0) {
      grupos <- grupos + 1
    }
  }
  
  # Configuración de los niveles de categorías para el gráfico
  catLevels <- c("noMin_noRegex", "MinFoldx", "RegexMatch", "MinFoldx_RegexMatch", 'LessThan5Res')
  df$Category <- factor(df$Category, labels = catLevels, levels = catLevels, ordered = TRUE)
  
  # Definición de colores para cada categoría
  colores <- c("noMin_noRegex" = "darkgray", 
               "MinFoldx" = 'orange2', 
               "RegexMatch" = 'blue', 
               "MinFoldx_RegexMatch" = 'yellow')
  
  # Creación de gráficos utilizando ggplot2
  pl = df %>%
    split(df$Group) %>%
    map(~ ggplot(data = .x) +
          geom_vline(data = subset(.x, 
                                   (Category == "MinFoldx_RegexMatch" | Category == "RegexMatch")),
                     aes(xintercept = MidPoint), 
                     col = 'darkred') +
          geom_point(stat = 'identity', 
                     color = "black", 
                     pch = 21, 
                     size = 3, 
                     aes_string(x = 'MidPoint', 
                                y = 'FoldX', 
                                fill = 'Category'), 
                     na.rm = TRUE) +
          geom_vline(xintercept = endTempPlot[between(endTempPlot, 
                                                      range(.x$MidPoint)[1], 
                                                      range(.x$MidPoint)[2])]) +
          geom_rect(data = subset(.x, 
                                  Category == "RegexMatch" | Category == "MinFoldx_RegexMatch"),
                    aes(xmin = MidPoint - nchar(Sub.Sequence) / 2, 
                        xmax = MidPoint + nchar(Sub.Sequence) / 2,
                        ymin = ylims[2] - 10, ymax = ylims[2] - 1),
                    col = 'darkred', 
                    linewidth = 0.5, 
                    fill = '#ff000066') +
          theme_plots_recall +
          scale_x_continuous(name = 'Position', 
                             expand = expansion(add = 1), 
                             breaks = seq(1, 740, 
                                          by = 5)) +
          geom_text(mapping = aes(x = MidPoint, 
                                  y = ylims[2] - 5, 
                                  label = Secuencia), 
                    na.rm = TRUE, 
                    family = 'Courier', 
                    fontface = 'bold') +
          scale_y_continuous(name = 'Valor FoldX', 
                             limits = ylims, 
                             breaks = seq(ylims[1], 
                                          ylims[2] - 2, 
                                          by = yBy), 
                             expand = c(0, 0)) +
          scale_fill_manual(name = "Category", 
                            values = colores, 
                            breaks = names(colores), 
                            drop = FALSE))
  
  # Ajuste del último gráfico para que se alinee correctamente
  n = nrow(df[df$Group == max(df$Group),])
  pl[[length(pl)]] = pl[[length(pl)]] + theme(plot.margin = margin(r = 0)) + plot_spacer() + plot_layout(widths = c(n, ultimaRow - n))
  
  # Combinación de todos los gráficos en uno solo
  patch <- wrap_plots(pl, ncol = 1, guides = "collect")
  patch <- patch + plot_annotation(title = plotTitle)
  
  # Definición de dimensiones para guardar el gráfico
  hh <- if(length(pl) == 1) 15 else length(pl) * 5
  ww <- 40
  
  # Guardado del gráfico 
  savePlot <- 'Rb/ELMfiles/DotPlots/'
  fileOut <- paste(FoldXDir, savePlot, fileOutName, sep = '')
  
  ggsave(filename = fileOut, plot = patch, device = "png", width = ww, height = hh, units = "cm", dpi = 300, limitsize = FALSE)
  
  # Devolver el dataframe.
  # return(foldXDF) 
}


#########   Retinoblastoma  ############

rb1GUX_9 <- plotTPELM(matrix = paste(rbDir, '1GUX_9', sep=''),
                       fileOutName = paste('/', "Rb_1GUX9.png",sep=''),
                       ylims = c(-4,70),
                       yBy= 16,
                       Ngroup= 4,
                      ultimaRow= 70,
                       plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1GUX_9")


rb1GUX_8 <- plotTPELM(matrix = paste(rbDir, '1GUX_8', sep=''),
                      fileOutName =paste('/', "Rb_1GUX8.png",sep=''),
                      ylims = c(-4,70),
                      yBy= 16,
                      Ngroup= 4,
                      ultimaRow= 100,
                      plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1GUX_8")

rb2R7G_10 <- plotTPELM(matrix = paste(rbDir, '2R7G_10', sep=''),
                      fileOutName =paste('/', "Rb_2R7G10.png",sep=''),
                      ylims = c(-4,40),
                      yBy= 10,
                      Ngroup= 3,
                      plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 2R7G_10")

rb2R7G_6 <- plotTPELM(matrix = paste(rbDir, '2R7G_6', sep=''),
                       fileOutName =paste('/', "Rb_2R7G6.png",sep=''),
                       ylims = c(-4,40),
                       yBy= 10,
                       Ngroup= 3,
                       plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 2R7G_6")

rb2R7G_5 <- plotTPELM(matrix = paste(rbDir, '2R7G_5', sep=''),
                       fileOutName =paste('/', "Rb_2R7G5.png",sep=''),
                       ylims = c(-4,40),
                       yBy= 10,
                       Ngroup= 3,
                      ultimaRow= 70,
                       plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 2R7G_5")

rb1N4M_9 <- plotTPELM(matrix = paste(rbDir, '1N4M_9', sep=''),
                       fileOutName =paste('/', "Rb_1N4M9.png",sep=''),
                       ylims = c(-4,40),
                       yBy= 10,
                       Ngroup= 3,
                      ultimaRow= 70,
                       plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1N4M_9")

rb1N4M_6 <- plotTPELM(matrix = paste(rbDir, '1N4M_6', sep=''),
                      fileOutName =paste('/', "Rb_1N4M6.png",sep=''),
                      ylims = c(-4,40),
                      yBy= 10,
                      Ngroup= 3,
                      plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1N4M_6")

rb1N4M_5 <- plotTPELM(matrix = paste(rbDir, '1N4M_5', sep=''),
                      fileOutName =paste('/', "Rb_1N4M5.png",sep=''),
                      ylims = c(-4,40),
                      yBy= 10,
                      Ngroup= 3,
                      ultimaRow= 70,
                      plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1N4M_5")




    ###   p107 ###

p1071GUX_9 <- plotTPELM(matrix = paste(p107Dir, '1GUX_9', sep=''),
                      fileOutName = paste('/', "p107_1GUX9.png",sep=''),
                      ylims = c(-4,70),
                      yBy= 16,
                      Ngroup= 4,
                      ultimaRow= 82,
                      plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1GUX_9")

p1071GUX_8 <- plotTPELM(matrix = paste(p107Dir, '1GUX_8', sep=''),
                        fileOutName = paste('/', "p107_1GUX8.png",sep=''),
                        ylims = c(-4,70),
                        yBy= 16,
                        Ngroup= 4,
                        ultimaRow= 120,
                        plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1GUX_8")

p1072R7G_10 <- plotTPELM(matrix = paste(p107Dir, '2R7G_10', sep=''),
                        fileOutName = paste('/', "p107_2R7G10.png",sep=''),
                        ylims = c(-4,40),
                        yBy= 10,
                        Ngroup= 2,
                        ultimaRow= 35,
                        plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 2R7G_10")

p1072R7G_6 <- plotTPELM(matrix = paste(p107Dir, '2R7G_6', sep=''),
                         fileOutName = paste('/', "p107_2R7G6.png",sep=''),
                         ylims = c(-4,40),
                         yBy= 10,
                         Ngroup= 2,
                        ultimaRow= 35,
                         plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 2R7G_6")

p1072R7G_5 <- plotTPELM(matrix = paste(p107Dir, '2R7G_5', sep=''),
                        fileOutName = paste('/', "p107_2R7G5.png",sep=''),
                        ylims = c(-4,40),
                        yBy= 10,
                        Ngroup= 2,
                        ultimaRow= 35,
                        plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 2R7G_5")

p1071N4M_9 <- plotTPELM(matrix = paste(p107Dir, '1N4M_9', sep=''),
                        fileOutName = paste('/', "p107_1N4M9.png",sep=''),
                        ylims = c(-4,40),
                        yBy= 10,
                        Ngroup= 2,
                        ultimaRow= 42,
                        plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1N4M_9")

p1071N4M_6 <- plotTPELM(matrix = paste(p107Dir, '1N4M_6', sep=''),
                        fileOutName = paste('/', "p107_1N4M6.png",sep=''),
                        ylims = c(-4,40),
                        yBy= 10,
                        Ngroup= 2,
                        ultimaRow= 62,
                        plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1N4M_6")

p1071N4M_5 <- plotTPELM(matrix = paste(p107Dir, '1N4M_5', sep=''),
                        fileOutName = paste('/', "p107_1N4M5.png",sep=''),
                        ylims = c(-4,40),
                        yBy= 10,
                        Ngroup= 2,
                        ultimaRow= 35,
                        plotTitle ="Valores minimos FoldX para interactores reportados en ELM escaneados con matriz 1N4M_5")
