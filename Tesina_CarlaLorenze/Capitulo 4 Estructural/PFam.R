#Plots de las tres pocket que muestra cuantos residuos de los peptidos que fueron hit, pertenecen a dominios Pfam
#Como input toma el archivo de lista de peptidos propd y se grafica la columna de 'Numero de residuos' incluida en la lista de datos. 
library(ggplot2)  # Carga la librería ggplot2 para la generación de gráficos.

# Define el directorio donde se encuentran los archivos con datos de pockets.
pocketDir <-'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/'

# Definición de un tema personalizado para los gráficos, ajustando fuentes, colores y márgenes.
tema_plots <- theme(
  axis.text.x = element_text(angle= 0, hjust = 0.5, vjust = 0.5, size = 20,
                             margin= margin(t=0.7, b=0.7, unit= "line"), colour = 'black'),
  axis.text.y = element_text(angle= 0, hjust = 0.5, vjust = 0.5, size = 20, colour = 'black'),
  axis.title = element_text(size = 25),
  axis.ticks = element_line(linewidth = 1, colour = 'black'),
  axis.text.y.right = element_text(margin= margin(r= 0.7,l=0.7, unit= "line")),
  axis.text.y.left = element_text(margin= margin(r= 0.7,l=0.7, unit= "line")),
  panel.border = element_rect(linewidth=2, colour = 'black'),
  panel.grid= element_blank()
)

# Definición de un valor máximo para la escala de los gráficos.
max_plot <- 198

# Función para generar y visualizar un histograma de residuos dentro de dominios Pfam.
plotPfam <- function(pocketFile, colorFill, fileName){
  
  # Construcción de la ruta completa del archivo a leer.
  archivo <- paste(pocketDir, pocketFile, sep='')
  
  # Lectura del archivo CSV con separador de punto y coma (;).
  pocketDF <- read.csv(file= archivo, sep = ";", header = TRUE)
  
  # Cálculo del porcentaje de péptidos con 8 o menos residuos dentro de dominios Pfam.
  porcentaje <- round((length(pocketDF$UniqueID[which(pocketDF$PFAM_residues <= 8)]) * 100) /
                        length(pocketDF$UniqueID), 0)
  
  # Impresión en consola del número de péptidos con 8 o menos residuos en dominios Pfam.
  cat('Npep con 8 o menos residuos solapados con dominios Pfam: ')
  cat(length(pocketDF$UniqueID[which(pocketDF$PFAM_residues <= 8)]))
  
  # Creación del histograma base.
  plot <- ggplot(pocketDF, aes(x= PFAM_residues)) +
    
    # Configuración del eje Y con límites y una segunda escala de frecuencia acumulada.
    scale_y_continuous(limits = c(0, 220),
                       breaks = seq(0, 202, by = 50),
                       expand = c(0,0),
                       sec.axis = sec_axis(trans= ~./max_plot,
                                           name = "Fr. acumulada",
                                           breaks = seq(0, 1, by= 0.2),
                                           labels= scales::percent))+
    
    # Configuración del eje X con límites y marcas cada 2 unidades.
    scale_x_continuous(limits = c(-1, 17),
                       breaks = seq(0, 17, by= 2), 
                       expand = c(0,0))+
    
    # Añade una región sombreada para resaltar los valores de interés (<=8 residuos).
    annotate('rect', xmin=-1, xmax=8, ymin=0, ymax=220, col='black',
             fill=colorFill, alpha=0.3) +
    
    # Creación del histograma con binwidth=1.
    geom_histogram(binwidth = 1, fill= colorFill, colour= 'black')+
    
    # Etiquetas de los ejes.
    labs(x= "Res. dentro de dominios Pfam", y= "Fr. absoluta")+
    
    # Añadir anotación del porcentaje de péptidos en la zona resaltada.
    annotate('text', x= 4, y= 170 , label= paste(porcentaje, '%', sep=''), size = 10)
  
  # Creación de un nuevo gráfico que agrega la frecuencia acumulada al histograma.
  plotfrecAcum <- plot + 
    stat_bin(aes(y=cumsum(after_stat(count))/sum(after_stat(count)) * max_plot),
             geom="step",
             binwidth = 2,
             center = 1,
             linewidth= 0.7)
  
  # Aplicación del tema definido previamente.
  plotfrecAcum <- plotfrecAcum + theme_bw() + tema_plots
  
  # Guardado del gráfico como archivo PNG.
  # ggsave(filename= fileName,
  #        plot = plotfrecAcum,
  #        device = "png",
  #        dpi = 600,
  #        width = 16,
  #        height = 12,
  #        units = "cm")
  
  return(pocketDF)  # Retorna el data frame leído para su posible uso posterior.
}



pfamRb <- plotPfam(pocketFile = 'ProP-PD_Rb_HD2.csv',
                   colorFill='#0087a7ff',
                   fileName= 'Pfam_Rb.png')


pfamp107 <- plotPfam(pocketFile = 'ProP-PD_p107_HD2.csv',
                     colorFill = '#C75C70',
                     fileName= 'Pfam_p107.png' )

length(pfamRb$PFAM_residues[which(pfamRb$PFAM_residues <= 8)]) #215
length(pfamp107$PFAM_residues[which(pfamp107$PFAM_residues <= 8)]) #73
