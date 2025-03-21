## Structural analysis.

#Alphafold accessibility for pocket.
#Plots de distribuicion de peptidos de acuerdo a su score de accesibilidad estructural. 
#Como input recibe un archivo propd, toma la columna de alphafold_accessibility y plotea frec absoluta de pep + curva de frecuencia acumulada. 

# Cargar la librería ggplot2 para la visualización de datos
library(ggplot2)

# Definir la ruta del directorio donde se encuentran los archivos de las pocket
pocketDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/'

# Definir un tema de gráfico personalizado para mejorar la presentación
tema_plots <- theme(axis.text.x = element_text(angle= 0,
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
                    axis.text.y.right = element_text(margin= margin(r= 0.7, l=0.7, unit= "line")),
                    axis.text.y.left = element_text(margin= margin(r= 0.7, l=0.7, unit= "line")),
                    panel.border = element_rect(size=2, colour = 'black'), 
                    panel.grid= element_blank())

# Definir un valor máximo para el eje secundario
max_alpha <- 30

# Función para crear un gráfico de accesibilidad usando datos de AlphaFold
alphaFoldPlot <- function(pocketFile, colorFill, fileName) {
  
  # Construir la ruta completa del archivo pocket
  file <- paste(pocketDir, pocketFile, sep='')
  # Leer el archivo CSV en un DataFrame
  pocketDF <- read.csv(file = file, sep = ";", header = TRUE)
  
  # Crear un nuevo DataFrame con la accesibilidad de AlphaFold
  aFoldDF <- data.frame("AlphaFold" = pocketDF$alphafold_accessibility)
  # Reemplazar valores negativos por 0
  aFoldDF$PosVal <- replace(aFoldDF$AlphaFold, which(aFoldDF$AlphaFold < 0), 0)
  
  # Calcular el porcentaje de péptidos con RSA mayor a 0.4
  porcentaje <- round((length(aFoldDF$PosVal[which(aFoldDF$PosVal > 0.4)]) * 100) / (length(aFoldDF$PosVal)), 0)
  cat('Porcentaje de péptidos con RSA > 0,4: ')
  cat((length(aFoldDF$PosVal[which(aFoldDF$PosVal > 0.4)])))
  
  # Crear el gráfico de accesibilidad
  alphafoldPlot <- ggplot(data = aFoldDF, aes(x = PosVal)) +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.2),
                       n.breaks = 4,
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 32), 
                       breaks = seq(0, 30, by = 10),
                       expand = c(0, 0),
                       sec.axis = sec_axis(trans = ~./max_alpha,
                                           name = "Fr. acumulada",
                                           breaks = seq(0, 1, by = 0.2),
                                           labels = scales::percent)) +
    annotate('rect', 
             xmin = 0.4,
             xmax = 1, 
             ymin = 0, 
             ymax = 32, 
             col = 'black',
             fill = colorFill,
             alpha = 0.3) +
    geom_histogram(binwidth = 0.02,
                   center = 0.01,
                   fill = colorFill, 
                   colour = 'black') +
    labs(x = "Accesibilidad", 
         y = "Fr. absoluta") +
    annotate('text',
             x = 0.65, 
             y = 25, 
             label = paste(porcentaje, '%', sep = ''),
             size = 10)
  
  # Añadir la frecuencia acumulada al gráfico
  frecAcum <- alphafoldPlot +  
    stat_bin(aes(y = cumsum(after_stat(count)) / sum(after_stat(count)) * max_alpha),
             geom = "step",
             binwidth = 0.02,
             center = 0.01,
             linewidth = 0.7) 
  
  # Aplicar tema y estilo al gráfico
  frecAcum <- frecAcum + theme_bw() + tema_plots
  
  # Guardar el gráfico como un archivo PNG
  ggsave(filename = fileName,
         plot = frecAcum,
         device = "png",
         dpi = 600,
         width = 16,
         height = 12,
         units = "cm")
}

#File Rb con datos Davey + Ylva HD2
Rb_afPlot <- alphaFoldPlot(pocketFile= 'ProP-PD_Rb_HD2.csv',
                           colorFill= '#0087a7ff',
                           fileName = 'RSA_Rb.png')


#File p107 con datos Davey + Ylva

p107_afPlot <- alphaFoldPlot(pocketFile= 'ProP-PD_p107_HD2.csv',
                           colorFill= '#C75C70',
                           fileName = 'RSA_p107.png')





# pocket <- read.csv(file = 'C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/IlvaVSDavey/ProPD_files/Propd_ListaDefinitivaDePeptidos/ProP-PD_p107.csv', sep = ";",header = TRUE)
