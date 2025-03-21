 #Funcion para hacer los histogramas stackeados, que despues no fueron presentados en la tesina de esta manera.  
stackHist <- function(regexFile,
                      foldXFilePath,
                      FileOut){
  
  regexPath <- paste(regexDir,regexFile, sep='')
  regexDF <- read.csv(file= regexPath, sep=';', header=T)
  
  foldXPath <- paste(FoldXDir,foldXFilePath, sep='')
  
  variantes <- data.frame('PepID'= regexDF$unique_id, "RegexMatch"= regexDF$Regex)
  variantes$Group <- NA
  
  
  for (i in 1:nrow(variantes)){
    
    if (startsWith(prefix = '(?=([IL].C.E.', x= variantes$RegexMatch[i])){
      variantes$Group[i] <- "LxCxE"
      
    }else if (startsWith(prefix = '(?=([IL].A.E.', x= variantes$RegexMatch[i])){
      variantes$Group[i] <- "LxAxE"
      
    }else if (startsWith(prefix = '(?=([IL].S.E.', x= variantes$RegexMatch[i])){
      variantes$Group[i] <- "LxSxE"
      
    }else if (startsWith(prefix = '(?=([IL].T.E.', x= variantes$RegexMatch[i])){
      variantes$Group[i] <- "LxTxE"
    } 
  }
  
  #Borro las filas con regex match de E2F
  variantes <- variantes[!is.na(variantes$Group),]

  # cat('variantesDF:')
  # cat(variantes$PepID[2])
  
  foldXDF <- data.frame("SeqID"= NA,
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
  
  files  <- list.files(path = foldXPath)

  # Dataframe para peptidos CON match de regex
  for (f in files){
    archivoFoldX <- paste(foldXPath, f, sep = '/')
    singleFile <- read.csv(archivoFoldX, header= T, sep='\t', skip = 4)

    df <- singleFile[which(singleFile$RegexMatch== 'True'), ]
    minimo <- min(df$FoldXperRes)

    df <- df[which(df$FoldXperRes==minimo),]

    foldXDF <- rbind(foldXDF, df)
  }
  foldXDF <- foldXDF[-1,]
  # #Check
  # cat(length(unique(foldXDF$SeqID)))
  # cat('\n')
  # cat(range(foldXDF$FoldX))

  foldXDF$VarianteLxCxE <- NA

   for (i in 1:nrow(variantes)){
     pepID <- variantes$PepID[i]
  
     for (j in 1:nrow(foldXDF)){
  
       if (foldXDF$SeqID[j] == pepID){
         foldXDF$VarianteLxCxE[j] <- variantes$Group[i]
       }
     }
   }
   

  
  
 max_plotRb <- 12

 catLevels <- c("LxCxE","LxAxE","LxSxE","LxTxE")
 foldXDF$VarianteLxCxE <- factor(foldXDF$VarianteLxCxE, labels =catLevels , levels= catLevels, ordered = T)
 
 colores <- c("LxCxE" = "#d45500","LxAxE"= '#00aa88',"LxSxE"= '#0044aa',"LxTxE"= '#e9afdd' )
 
 foldxPlot <- ggplot(foldXDF,
                     aes(x= FoldX,
                         fill = VarianteLxCxE))+
   scale_y_continuous(limits = c(0, 13),
                      breaks = seq(0,12, by = 3),
                      expand = c(0,0),
                      sec.axis = sec_axis(trans= ~./max_plotRb,
                                          name = "Fr. acumulada",
                                          labels= scales::percent))+
   scale_x_continuous(limits = c(-6,12),
                      breaks = seq(-6, 12, by= 3),
                      expand = c(0,0))+
   geom_histogram(binwidth = 0.5,
                  center= 0.25,
                  colour= 'black')+
   labs(x= "Valor de FoldX",
        y= "Fr. absoluta")+
   geom_vline(xintercept = 3,
              linetype= 'dashed',
              color= 'black',
              linewidth=1)+
   scale_fill_manual(name= "VarianteLxCxE", 
                     values = colores, 
                     breaks = names(colores),
                     drop= F)
 
 frecAcum <- foldxPlot +
   stat_bin(aes(y=cumsum(after_stat(count))/sum(after_stat(count))*max_plotRb, fill= NULL),
              geom="step",
              binwidth = 0.5,
              center =0.25,
              linewidth= 0.7)
   
   
   frecAcum <- frecAcum+ tema_plots + 
     theme(legend.position = "bottom",
           legend.title = element_blank(),
           legend.spacing.y = unit(2,'line'),
           legend.text = element_text(size = 15))
   
  ggsave(filename= FileOut,
         plot =frecAcum,
         device = "png",
         dpi = 600,
         width = 16,
         height = 10,
         units = "cm")
 #return(foldXDF)
  
}

## Distribución de peptidos (hits Rb hd2, Rb hd+hd2 y p107) de acuerdo a estabilidad energetica escaneados con matrices FoldX 1GUX, 2R7G, 1N4M y 4YOS

# Toma como input archivos escaneados con matrices FoldX y plotea los histogramas.

# Carga de bibliotecas necesarias para la visualización
library(ggplot2)      # Para crear gráficos
library(svglite)      # Para guardar gráficos en formato SVG

# Definición de las rutas de los directorios donde se encuentran los archivos de FoldX
FoldXDir <- 'C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/FoldX/'
rbE2FPath <- 'Retinoblastoma/HD2/FoldXScan/E2F'
rbLxcxePath <- 'Retinoblastoma/HD2/FoldXScan/LxCxE'
rbSRpath <- 'Retinoblastoma/HD2/FoldXScan/SinRegex'
p107e2fPath <- 'p107/E2F/FoldXScan'
p107lxcxePath <- 'p107/LxCxE/FoldXScan'
p107SRpath <- 'p107/SinRegex/FoldXScan'

# Configuración del tema para los gráficos
tema_plots <- theme_bw() + theme(
  axis.text.x = element_text(angle = 0, 
                             hjust = 0.5, 
                             vjust = 0.5, 
                             size = 20, 
                             colour = 'black'),
  axis.text.y = element_text(angle = 0, 
                             hjust = 0.5, 
                             vjust = 0.5, 
                             size = 20, 
                             colour = 'black'),
  axis.title = element_text(size = 25),
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

# Función para procesar datos de péptidos que coinciden con una expresión regular
matchPep <- function(FoldXFilePath, umbralMotivo) {
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
  
  # Ruta del archivo de FoldX
  filePath <- paste(FoldXDir, FoldXFilePath, sep = "/")
  files <- list.files(path = filePath)
  
  # Procesamiento de los archivos de FoldX para péptidos con coincidencia de regex
  for (f in files) {
    archivoFoldX <- paste(filePath, f, sep = '/')
    singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
    
    # Filtrado de datos que cumplen la condición de coincidencia de regex
    df <- singleFile[which(singleFile$RegexMatch == 'True'), ]
    
    # Selección del mínimo valor de FoldX
    minimo <- min(df$FoldX)
    df <- df[which(df$FoldX == minimo), ]
    
    # Agregación de los datos seleccionados al dataframe
    FoldXDF <- rbind(FoldXDF, df)
  }
  
  # Eliminación de la primera fila vacía
  FoldXDF <- FoldXDF[-1,]
  
  # Cálculo de la cantidad y porcentaje de péptidos únicos con FoldX menor o igual al umbral
  nPep <- length(unique(FoldXDF$SeqID[which(FoldXDF$FoldX <= umbralMotivo)]))
  porcPep <- nPep / nrow(FoldXDF) * 100
  
  # Impresión de resultados
  cat(range(FoldXDF$FoldX))
      cat('\n')
      cat('Cantidad de péptidos únicos con FoldX <= ')
      cat(umbralMotivo)
      cat(' = ')
      cat(nPep)
      cat('\n')
      cat('\n')
      cat('Porcentaje: ')
      cat(porcPep)
      cat('%')
      cat('\n')
      
      return(FoldXDF)  # Devuelve el dataframe de FoldX
}

# Función para procesar datos de péptidos sin coincidencia de expresión regular
nonMatchPep <- function(FoldXFilePath, umbralMotivo) {
  
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
  
  # Ruta del archivo de FoldX
  filePath <- paste(FoldXDir, FoldXFilePath, sep = "/")
  files <- list.files(path = filePath)
  
  # Procesamiento de los archivos de FoldX para péptidos sin coincidencia de regex
  for (f in files) {
    archivoFoldX <- paste(filePath, f, sep = '/')
    singleFile <- read.csv(archivoFoldX, header = TRUE, sep = '\t', skip = 4)
    
    # Filtrado de datos que cumplen las condiciones
    df <- singleFile[which(singleFile$MoreThan5 == 'True' & singleFile$Min == 'True'), ]
    
    # Selección del mínimo valor de FoldX
    minimo <- min(df$FoldX)
    df <- df[which(df$FoldX == minimo), ]
    
    # Agregación de los datos seleccionados al dataframe
    FoldXDF <- rbind(FoldXDF, df)
  }
  
  # Eliminación de la primera fila vacía
  FoldXDF <- FoldXDF[-1,]
  
  # Cálculo de la cantidad de péptidos únicos con FoldX mayor o igual al umbral
  nPep <- length(unique(FoldXDF$SeqID[which(FoldXDF$FoldX >= umbralMotivo)]))
  
  # Impresión de resultados
  cat(range(FoldXDF$FoldX))
      cat('\n')
      cat('Cantidad de péptidos únicos con FoldX >= ')
      cat(umbralMotivo)
      cat(' = ')
      cat(nPep)
      cat('\n')
      
      return(FoldXDF)  # Devuelve el dataframe de FoldX
}

# Función para crear un histograma de FoldX
histograma <- function(df, 
                       FileOutName, 
                       yLims, 
                       xintercept, 
                       yBy, 
                       umbralMotivo, 
                       pocketColor) {
  
  max_plot <- yLims[2] - 0.1 * yLims[2]  # Ajuste para el gráfico
  
  # Creación del gráfico de histograma
  foldxPlot <- ggplot(df, 
                      aes(x = FoldX)) + 
    scale_y_continuous(expand = c(0, 0), 
                       limits = yLims, 
                       breaks = seq(yLims[1], 
                                    yLims[2], 
                                    by = yBy), 
                       sec.axis = sec_axis(trans = ~./max_plot, 
                                           name = "Fr. acumulada", 
                                           breaks = seq(0, 1, by = 0.2), 
                                           labels = scales::percent)) +
    geom_vline(xintercept = xintercept, 
               linetype = 'dashed', 
               color = '#7393B3', 
               linewidth = 1) +
    scale_x_continuous(limits = c(-2, 16), 
                       breaks = seq(-2, 16, by = 2), 
                       expand = c(0, 0)) +
    geom_histogram(binwidth = 0.25, center = 0.125, 
                   fill = pocketColor, 
                   colour = 'black') +
    labs(x = "Valor de FoldX", 
         y = "Fr. absoluta")
  
  # Adición de frecuencias acumuladas al gráfico
  frecAcum <- foldxPlot + 
    stat_bin(aes(y = cumsum(after_stat(count)) / sum(after_stat(count)) * max_plot), 
             geom = "step", 
             binwidth = 0.25, 
             center = 0.125, 
             linewidth = 0.7)
  
  frecAcum <- frecAcum + 
    theme_bw() + 
    tema_plots
  
  # Guardado del gráfico
  ggsave(filename = FileOutName, 
         plot = frecAcum, 
         device = "png", 
         dpi = 600, 
         width = 16, 
         height = 9, 
         units = "cm")
  
  return(frecAcum)  # Devuelve el gráfico
}

#Uso estas "funciones" para agrupar los programas de cada pocket y de cada matriz para que el codigo no sea tan largo. NO SON FUNCIONES
datosRB <- function(nada){
  
  ## Péptidos Rb escaneados con 1GUX
  
  
  lxcxe_1GUX9 <- matchPep(umbralMotivo= 5,
                          FoldXFilePath= paste(rbLxcxePath, '1GUX_9', sep= '/'))
  
  lxcxe_1GUX8 <- matchPep(umbralMotivo= 5,
                          FoldXFilePath=paste(rbLxcxePath, '1GUX_8', sep= '/'))
  
  sR_1GUX9 <- nonMatchPep(umbralMotivo= 5,
                          FoldXFilePath=paste(rbSRpath, '1GUX_9', sep= '/'))
  
  sR_1GUX8 <- nonMatchPep(umbralMotivo= 5,
                          FoldXFilePath= paste(rbSRpath, '1GUX_8', sep= '/'))
  rownames(sR_1GUX8) <- NULL
  sR_1GUX8 <- sR_1GUX8[-82,]
  
  
  
  E2F_2R7G10 <- matchPep(FoldXFilePath= paste(rbE2FPath, '2R7G_10', sep= '/'),
                         umbralMotivo= 3)
  
  E2F_2R7G6<- matchPep(FoldXFilePath= paste(rbE2FPath, '2R7G_6', sep= '/'),
                       umbralMotivo= 3)
  
  E2F_2R7G5<- matchPep(FoldXFilePath= paste(rbE2FPath, '2R7G_5', sep= '/'),
                       umbralMotivo= 3)
  
  sR_2R7G10<- nonMatchPep(FoldXFilePath= paste(rbSRpath, '2R7G_10', sep= '/'),
                          umbralMotivo= 3)
  
  sR_2R7G6<- nonMatchPep(FoldXFilePath= paste(rbSRpath, '2R7G_6', sep= '/'),
                         umbralMotivo= 3)
  
  sR_2R7G5<- nonMatchPep(FoldXFilePath= paste(rbSRpath, '2R7G_5', sep= '/'),
                         umbralMotivo= 3)
  rownames(sR_2R7G5) <- NULL
  sR_2R7G5 <- sR_2R7G5[-c(136),]
  
  
  
  E2F_1N4M9 <- matchPep(FoldXFilePath= paste(rbE2FPath, '1N4M_9', sep= '/'),
                        umbralMotivo= 3)
  
  E2F_1n4m6<- matchPep(FoldXFilePath= paste(rbE2FPath, '1N4M_6', sep= '/'),
                       umbralMotivo= 3)
  
  E2F_1n4m5<- matchPep(FoldXFilePath= paste(rbE2FPath, '1N4M_5', sep= '/'),
                       umbralMotivo= 3)
  
  sR_1n4m9<- nonMatchPep(FoldXFilePath= paste(rbSRpath, '1N4M_9', sep= '/'),
                         umbralMotivo= 3)
  
  sR_1n4m6<- nonMatchPep(FoldXFilePath= paste(rbSRpath, '1N4M_6', sep= '/'),
                         umbralMotivo= 3)
  
  sR_1n4m5<- nonMatchPep(FoldXFilePath= paste(rbSRpath, '1N4M_5', sep= '/'),
                         umbralMotivo= 3)
}
histogramas_Retinoblastoma <- function(nada){
  

lxcxe_1GUX9_histograma <- histograma(df= lxcxe_1GUX9,
                                     FileOutName = 'Rblxcxe_1GUX9.png',
                                     yLims= c(0,14),
                                     yBy= 2,
                                     xintercept= 5,
                                     umbralMotivo= 5,
                                     pocketColor ='#0087a7ff')
lxcxe_1GUX9_histograma


lxcxe_1GUX8_histograma <- histograma(df= lxcxe_1GUX8,
                                     FileOutName = 'Rblxcxe_1GUX8.png',
                                     yLims= c(0,14),
                                     yBy= 2,
                                     xintercept= 5,
                                     umbralMotivo= 5,
                                     pocketColor ='#0087a7ff')
lxcxe_1GUX8_histograma


##--- SIN REGEX ----    ##
sR_1GUX9_histograma <- histograma(df= sR_1GUX9,
                                  FileOutName = 'RbsR_1GUX9.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 5,
                                  umbralMotivo= 5,
                                  pocketColor ='#0087a7ff')




sR_1GUX8_histograma <- histograma(df= sR_1GUX8,
                                  FileOutName = 'RbsR_1GUX8.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 5,
                                  umbralMotivo= 5,
                                  pocketColor ='#0087a7ff')



## Péptidos ** Rb ** CON match E2F escaneados con 2R7G y variantes


E2F_2R7G10_histograma <- histograma(df= E2F_2R7G10,
                                    FileOutName = 'RbE2F2_2r7g10.png',
                                    yLims= c(0,14),
                                    yBy= 2,
                                    xintercept= 3,
                                    umbralMotivo= 3,
                                    pocketColor ='#0087a7ff')



E2F_2R7G6_histograma <- histograma(df= E2F_2R7G6,
                                    FileOutName = 'RbE2F2_2r7g6.png',
                                   yLims= c(0,14),
                                   yBy= 2,
                                   xintercept= 3,
                                   umbralMotivo= 3,
                                    pocketColor ='#0087a7ff')



E2F_2R7G5_histograma <- histograma(df= E2F_2R7G5,
                                   FileOutName = 'RbE2F2_2r7g5.png',
                                   yLims= c(0,14),
                                   yBy= 2,
                                   xintercept= 3,
                                   umbralMotivo= 3,
                                   pocketColor ='#0087a7ff')



##--- SIN REGEX ----    ##


sR_2R7G10_histograma <- histograma(df= sR_2R7G10,
                                   FileOutName = 'RbsR_2r7g10.png',
                                   yLims= c(0,14),
                                   yBy= 2,
                                   xintercept= 3,
                                   umbralMotivo= 3,
                                   pocketColor ='#0087a7ff')



sR_2R7G6_histograma <- histograma(df= sR_2R7G6,
                                   FileOutName = 'RbsR_2r7g6.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 3,
                                  umbralMotivo= 3,
                                   pocketColor ='#0087a7ff')



sR_2R7G5_histograma <- histograma(df= sR_2R7G5,
                                  FileOutName = 'RbsR_2r7g5.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 3,
                                  umbralMotivo= 3,
                                  pocketColor ='#0087a7ff')




##      Péptidos Rb CON match E2F escaneados con 1N4M y variantes     ##

E2F_1N4M9_histograma <- histograma(df= E2F_1N4M9,
                                    FileOutName = 'RbE2F2_1N4M9.png',
                                   yLims= c(0,14),
                                   yBy= 2,
                                   xintercept= 3,
                                   umbralMotivo= 3,
                                    pocketColor ='#0087a7ff')



E2F_1n4m6_histograma <- histograma(df= E2F_1n4m6,
                                   FileOutName = 'RbE2F2_1n4m6.png',
                                   yLims= c(0,14),
                                   yBy= 2,
                                   xintercept= 3,
                                   umbralMotivo= 3,
                                   pocketColor ='#0087a7ff')




E2F_1n4m5_histograma <- histograma(df= E2F_1n4m5,
                                   FileOutName = 'RbE2F2_1n4m5.png',
                                   yLims= c(0,14),
                                   yBy= 2,
                                   xintercept= 3,
                                   umbralMotivo= 3,
                                   pocketColor ='#0087a7ff')



    ## SIN REGEX     ##


sR_1n4m9_histograma <- histograma(df= sR_1n4m9,
                                   FileOutName = 'RBsR_1n4m9.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 3,
                                  umbralMotivo= 3,
                                   pocketColor ='#0087a7ff')



sR_1n4m6_histograma <- histograma(df= sR_1n4m6,
                                  FileOutName = 'RcsR_1n4m6.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 3,
                                  umbralMotivo= 3,
                                  pocketColor ='#0087a7ff')


sR_1n4m5_histograma <- histograma(df= sR_1n4m5,
                                  FileOutName = 'RcsR_1n4m5.png',
                                  yLims= c(0,14),
                                  yBy= 2,
                                  xintercept= 3,
                                  umbralMotivo= 3,
                                  pocketColor ='#0087a7ff')


}

##   p107   ##

datosP107 <- function(nada){
  
  
  p107E2F_2R7G10 <- matchPep(FoldXFilePath= paste(p107e2fPath, '2R7G_10', sep= '/'),
                             umbralMotivo= 3)
  
  
  
  p107E2F_2R7G6 <- matchPep(FoldXFilePath= paste(p107e2fPath, '2R7G_6', sep= '/'),
                            umbralMotivo= 3)
  

  
  p107E2F_2R7G5<- matchPep(FoldXFilePath= paste(p107e2fPath, '2R7G_5', sep= '/'),
                           umbralMotivo= 3)
  
  
  ##--- SIN REGEX ----    ##
  
  
  p107sR_2R7G10 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '2R7G_10', sep= '/'),
                               umbralMotivo= 3)
  
  
  
  p107sR_2R7G6 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '2R7G_6', sep= '/'),
                              umbralMotivo= 3)
  
  
  
  
  p107sR_2R7G5 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '2R7G_5', sep= '/'),
                              umbralMotivo= 3)
  rownames(p107sR_2R7G5) <- NULL
  p107sR_2R7G5 <- p107sR_2R7G5[-10,]
  
  
  
  
  ## Péptidos p107 CON escaneados con 1N4M
  
  p107E2F_1n4m9 <- matchPep(FoldXFilePath= paste(p107e2fPath, '1N4M_9', sep= '/'),
                            umbralMotivo= 3)
  
  
  
  p107E2F_1n4m6 <- matchPep(FoldXFilePath= paste(p107e2fPath, '1N4M_6', sep= '/'),
                            umbralMotivo= 3)
  
  
  p107E2F_1n4m5 <- matchPep(FoldXFilePath= paste(p107e2fPath, '1N4M_5', sep= '/'),
                            umbralMotivo=3)
  
  
  # ##---     SIN REGEX     ---##
  
  
  p107sR_1n4m9 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '1N4M_9', sep= '/'),
                              umbralMotivo= 3)
  
  
  p107sR_1n4m6 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '1N4M_6', sep= '/'),
                              umbralMotivo=3)
  rownames(p107sR_1n4m6) <- NULL
  p107sR_1n4m6 <- p107sR_1n4m6[-22,]
  
  
  p107sR_1n4m5 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '1N4M_5', sep= '/'),
                              umbralMotivo= 3)
  
  
  # 
  # ## Péptidos p107  escaneados con 1GUX
  
  p107lxcxe_1gux9 <- matchPep(FoldXFilePath= paste(p107lxcxePath, '1GUX_9', sep= '/'),
                              umbralMotivo= 5)
  
  
  
  p107lxcxe_1gux8 <- matchPep(FoldXFilePath= paste(p107lxcxePath, '1GUX_8', sep= '/'),
                              umbralMotivo= 5)
  
  
  # ##--- SIN REGEX ----    ##
  
  p107sR_1GUX9 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '1GUX_9', sep= '/'),
                              umbralMotivo= 5)
  
  
  p107sR_1GUX8 <- nonMatchPep(FoldXFilePath= paste(p107SRpath, '1GUX_8', sep= '/'),
                              umbralMotivo= 5)
  
}
histogramas_p107 <- function(nadaDeNada){
  
  p107E2F_2R7G10_histograma <- histograma(df= p107E2F_2R7G10,
                                          FileOutName = 'p107E2F2_2r7g10.png',
                                          yLims= c(0,6),
                                          yBy= 2,
                                          xintercept= 3,
                                          umbralMotivo= 3,
                                          pocketColor ='#C75C70')
  
  
  p107E2F_2R7G6_histograma <- histograma(df= p107E2F_2R7G6,
                                         FileOutName = 'p107E2F2_2r7g6.png',
                                         yLims= c(0,6),
                                         yBy= 2,
                                         xintercept= 3,
                                         umbralMotivo= 3,
                                         pocketColor ='#C75C70')
  
  p107E2F_2R7G5_histograma <- histograma(df= p107E2F_2R7G5,
                                         FileOutName = 'p107E2F2_2r7g5.png',
                                         yLims= c(0,6),
                                         yBy= 2,
                                         xintercept= 3,
                                         umbralMotivo= 3,
                                         pocketColor ='#C75C70')

  
  p107sR_2R7G10_histograma <- histograma(df= p107sR_2R7G10,
                                         FileOutName = 'p107sR_2r7g10.png',
                                         yLims= c(0,6),
                                         yBy= 2,
                                         xintercept= 3,
                                         umbralMotivo= 3,
                                         pocketColor ='#C75C70')
 
  p107sR_2R7G6_histograma <- histograma(df= p107sR_2R7G6,
                                        FileOutName = 'p107sR_2r7g6.png',
                                        yLims= c(0,6),
                                        yBy= 2,
                                        xintercept= 3,
                                        umbralMotivo= 3,
                                        pocketColor ='#C75C70')
  
  
  p107sR_2R7G5_histograma <- histograma(df= p107sR_2R7G5,
                                        FileOutName = 'p107sR_2r7g5.png',
                                        yLims= c(0,6),
                                        yBy= 2,
                                        xintercept= 3,
                                        umbralMotivo= 3,
                                        pocketColor ='#C75C70')
  
  
  p107E2F_1n4m9_histograma <- histograma(df= p107E2F_1n4m9,
                                         FileOutName = 'p107E2F2_1n4m9.png',
                                         yLims= c(0,6),
                                         yBy= 2,
                                         xintercept= 3,
                                         umbralMotivo= 3,
                                         pocketColor ='#C75C70')
  
  
  p107E2F_1n4m6_histograma <- histograma(df= p107E2F_1n4m6,
                                         FileOutName = 'p107E2F2_1n4m6.png',
                                         yLims= c(0,6),
                                         yBy= 2,
                                         xintercept= 3,
                                         umbralMotivo= 3,
                                         pocketColor ='#C75C70')
  
  
  p107E2F_1n4m5_histograma <- histograma(df= p107E2F_1n4m5,
                                         FileOutName = 'p107E2F2_1n4m5.png',
                                         yLims= c(0,6),
                                         yBy= 2,
                                         xintercept= 3,
                                         umbralMotivo= 3,
                                         pocketColor ='#C75C70')
  
  
  
  p107sR_1n4m9_histograma <- histograma(df= p107sR_1n4m9,
                                        FileOutName = 'p107sR_1n4m9.png',
                                        yLims= c(0,6),
                                        yBy= 2,
                                        xintercept= 3,
                                        umbralMotivo= 3,
                                        pocketColor ='#C75C70')
  
  
  p107sR_1n4m6_histograma <- histograma(df= p107sR_1n4m6,
                                        FileOutName = 'p107sR_1n4m6.png',
                                        yLims= c(0,6),
                                        yBy= 2,
                                        xintercept= 3,
                                        umbralMotivo= 3,
                                        pocketColor ='#C75C70')
  
  
  p107sR_1n4m5_histograma <- histograma(df= p107sR_1n4m5,
                                        FileOutName = 'p107sR_1n4m5.png',
                                        yLims= c(0,6),
                                        yBy= 2,
                                        xintercept= 3,
                                        umbralMotivo= 3,
                                        pocketColor ='#C75C70')
  
  
  p107lxcxe_1gux9_histograma <- histograma(df= p107lxcxe_1gux9,
                                           FileOutName = 'p107lxcxe_1gux9.png',
                                           yLims= c(0,8),
                                           yBy= 2,
                                           xintercept= 5,
                                           umbralMotivo= 5,
                                           pocketColor ='#C75C70')

  
  p107lxcxe_1gux8_histograma <- histograma(df= p107lxcxe_1gux8,
                                           FileOutName = 'p107lxcxe_1gux8.png',
                                           yLims= c(0,8),
                                           yBy= 2,
                                           xintercept= 5,
                                           umbralMotivo= 5,
                                           pocketColor ='#C75C70')

  
  p107sR_1GUX8_histograma <- histograma(df= p107sR_1GUX8,
                                        FileOutName = 'p107sR_1GUX8.png',
                                        yLims= c(0,8),
                                        yBy= 2,
                                        xintercept= 5,
                                        umbralMotivo= 5,
                                        pocketColor ='#C75C70')
  
  

  
  
  p107sR_1GUX9_histograma <- histograma(df= p107sR_1GUX9,
                                        FileOutName = 'p107sR_1GUX9.png',
                                        yLims= c(0,8),
                                        yBy= 2,
                                        xintercept= 5,
                                        umbralMotivo= 5,
                                        pocketColor ='#C75C70')
  
  

  
}
