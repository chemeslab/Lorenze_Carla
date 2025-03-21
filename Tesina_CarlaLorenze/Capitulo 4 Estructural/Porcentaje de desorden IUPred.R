#Graficos del porcentaje de score promedio de prediccion del desorden (IUPred, longrun) de los peptidos que fueron hit en propd para cada pocket como carnada. 

#Se corrio IUPred para las secuencias y te tomaron los scores de los 16 residuos que componen el peptido unicamente. Se tomo cuantos residuos obtuvieron un score mayor a 0,4 para estimar el porcentaje de desorden del peptido. 

# Cargar la librería ggplot2 para la visualización de datos
library(ggplot2)

# Definir la ruta del directorio donde se encuentran los archivos de IUPred
iupredDir <- "C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/Estructural/IUPred scores_Propd"

# Definir un tema de gráfico personalizado para mejorar la presentación
plotTheme <- theme(axis.text.x = element_text(angle= 0,
                                              hjust = 0.5,
                                              vjust = 0.5,
                                              size = 20,
                                              margin= margin(t=0.7, b=0.7, unit= "line"),
                                              colour = 'black'),
                   axis.text.y = element_text(angle= 0,
                                              hjust = 0.5,
                                              vjust = 0.5,
                                              size = 20,
                                              colour = 'black'),
                   axis.title = element_text(size = 25),
                   axis.ticks = element_line(linewidth = 1, colour = 'black'),
                   axis.text.y.right = element_text(margin= margin(r= 0.7,l=0.7, unit= "line")),
                   axis.text.y.left = element_text(margin= margin(r= 0.7,l=0.7, unit= "line")),
                   panel.border = element_rect(size=2, colour = 'black'), 
                   panel.grid= element_blank())

# Función para crear un histograma del porcentaje de desorden en un archivo de pocket
pepDisorder <- function(pocketFile, histColor, fileName) {
  
  # Construir la ruta completa del archivo de pocket
  pocketPath <- paste(iupredDir, pocketFile, sep="/")
  
  # Leer el archivo CSV en un DataFrame
  pocketDF <- read.csv(pocketPath, sep = ";", header = TRUE)
  
  # Calcular el porcentaje de péptidos con más del 50% de residuos desordenados
  porcentaje <- round((length(pocketDF$Accession_number[which(pocketDF$Disorder_. >= 50)]) * 100) / (length(pocketDF$Accession_number)), 0)
  cat('Npep con > 50% de sus residuos desordenados: ')
  cat(length(pocketDF$Accession_number[which(pocketDF$Disorder_. >= 50)]))
  
  # Crear el histograma del porcentaje de desorden
  plot <- ggplot(pocketDF, mapping = aes(x = Disorder_.)) + 
    geom_histogram(fill = histColor, 
                   colour = 'black', 
                   binwidth = 5, 
                   center = 2.5) + 
    theme_bw() + 
    labs(x = "Porcentaje de desorden", y = "Fr. absoluta") +
    scale_x_continuous(limits = c(0, 100), 
                       breaks = seq(0, 100, by = 20),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 95),
                       breaks = seq(0, 95, by = 20),
                       expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~./max_plot,
                                           name = "Fr. acumulada",
                                           breaks = seq(0, 1, by = 0.2),
                                           labels = scales::percent)) +
    plotTheme
  
  # Establecer el valor máximo para el eje secundario del gráfico
  max_plot <- 90
  
  # Añadir la frecuencia acumulada al gráfico
  plot_frecAcum <- plot + 
    stat_bin(aes(y = cumsum(after_stat(count)) / sum(after_stat(count)) * max_plot),
             geom = "step",
             binwidth = 5,
             center = 2.5,
             linewidth = 0.7) +
    annotate('rect', 
             xmin = 50,
             xmax = 100, 
             ymin = 0, 
             ymax = 95, 
             col = 'black',
             fill = histColor,
             alpha = 0.3) +
    annotate('text',
             x = 75, 
             y = 80, 
             label = paste(porcentaje, '%', sep=''),
             size = 10)
  
  # Guardar el gráfico como un archivo PNG
  ggsave(filename = fileName,
         plot = plot_frecAcum,
         device = "png",
         dpi = 600,
         width = 16,
         height = 12,
         units = "cm")
}

#Porcentaje de desorden Rb

rbPlot <- pepDisorder(pocketFile = '20230403_120140_IUPred_scores_Rb_HD2.csv',
                        histColor ='#0087a7ff',
                        fileName= "IUPred_Rb.png")


#Porcentaje de desorden p107

p107Plot <- pepDisorder(pocketFile = '20230403_120140_IUPred_scores_p107.csv',
                        histColor ='#C75C70',
                        fileName= "IUPred_p107.png")
