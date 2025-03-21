library(ggplot2) #Cargo libreria para graficar
library(dplyr) #Cargo libreria para manipular datos

# Cargo el archivo HD2
hd2File <- "C:/Users/Carla/Desktop/Tesina/GitHub/phageD/Data/files_all/HD2 umbral ELM 050-intact-proteomic sin isoformas/20230308_140221_HD2_050.csv"
hd2df <- read.csv(file= hd2File, sep=';', header=T)
# Cargo el archivo Rb
rbFile <- "C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/ProP-PD_Rb_HD2.csv"

# Cargo el archivo p107
p107File <- "C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/ProP-PD_p107_HD2.csv"

# Cargo el archivo p130
p130File <- "C:/Users/Carla/Desktop/Propd_ListaDefinitivaDePeptidos/ProP-PD_p130_HD2.csv"
# p130_df <- read.csv(file= p130File, sep=';', header=T)

# p130_df <- p130_df %>%
#   left_join(hd2df %>% select(UniqueID, 
#                              LxCxE, 
#                              E2F_mimic, 
#                              LxSxE), 
#             by = 'UniqueID') 


#Defino el conjunto de factores que voy a usar para graficar
twoMotifFactor <- c("LxCxE", "E2F")
allMotifFactor <- c("LxCxE", "E2F", "LxSxE")

tema_plots <- theme_bw()+ 
  theme(axis.text.x = element_text(angle= 0,
                                               hjust = 0.5,
                                               vjust = 0.5,
                                               size = 32,
                                               margin= margin(t=0.7, b=0.7, unit= "line"),
                                               colour = 'black'),
                    axis.text.y = element_text(angle= 0,
                                               hjust = 0.5,
                                               vjust = 0.5,
                                               size = 30,
                                               colour = 'black'),
                    axis.title = element_text(size = 35),
                    axis.ticks = element_line(linewidth = 1, colour = 'black'),
                    axis.text.y.right = element_text(margin= margin(r= 0.5,
                                                                    l=0.5, 
                                                                    unit= "line")),
                    axis.text.y.left = element_text(margin= margin(r= 0.5,l=0.5, unit= "line")),
                    panel.border = element_rect(linewidth =2, colour = 'black'),
                    panel.grid= element_blank())


##  Creo un DF de TP de ELM para cada set de datos
# Función para calcular la frecuencia de proteínas con SLiMs en un archivo de datos
frecuencyTPELM <- function(dataFile){
  
 readFile <- read.csv(file= dataFile, sep=';', header=T)
  # readFile <- p130_df
 # Crear un nuevo data frame con columnas relevantes
 ELM_DF <- data.frame('Accession_n'= readFile$Accession,
                      "Unique_ID"= readFile$UniqueID,
                      "LxCxE"= readFile$LxCxE, 
                      "E2F_mimic"= readFile$E2F_mimic, 
                      "LxSxE"= readFile$LxSxE)
 
 # Convertir valores en booleanos y luego en etiquetas de texto
 ELM_DF$LxCxE <- as.logical(ELM_DF$LxCxE)
 ELM_DF$E2F_mimic <- as.logical(ELM_DF$E2F_mimic)
 ELM_DF$LxSxE <- as.logical(ELM_DF$LxSxE)
 
 ELM_DF$LxCxE[which(ELM_DF$LxCxE)] <- "LxCxE"
 ELM_DF$E2F_mimic[which(ELM_DF$E2F_mimic)] <- "E2F"
 ELM_DF$LxSxE[which(ELM_DF$LxSxE)] <- "LxSxE"
 
 # Crear una columna combinada con los motivos presentes
 ELM_DF$Motif <- paste(ELM_DF$LxCxE, ELM_DF$E2F_mimic, ELM_DF$LxSxE, sep=",")
 ELM_DF$Motif <- gsub(pattern=",FALSE",replacement= '',x=ELM_DF$Motif)
 ELM_DF$Motif <- gsub(pattern="FALSE,",replacement= '',x=ELM_DF$Motif)
 ELM_DF$Motif <- gsub(pattern="FALSE",replacement='Not in ELM',x=ELM_DF$Motif)
 
 # Calcular la frecuencia de proteínas para cada motivo
 lxcxeProt <- length(unique(ELM_DF$Accession_n[which(ELM_DF$Motif == 'LxCxE')]))
 e2fProt <-length(unique(ELM_DF$Accession_n[which(ELM_DF$Motif == 'E2F')]))
 lxsxeProt <-length(unique(ELM_DF$Accession_n[which(ELM_DF$Motif == 'LxSxE')]))
 
 # Crear un data frame con los resultados
 frecDF <- data.frame('Motivo'= c('LxCxE', 'E2F', "LxSxE"), 
                      'Frecuencia' = c(lxcxeProt, e2fProt, lxsxeProt))
  
  return(frecDF)
}


# Función para graficar la frecuencia de SLiMs
pocketRecallPlot <- function(hd2DF,
                         pocketDF,
                         factor,
                         pocketColor,
                         #pocketTitle,
                         fileOut){
  
hd2DF$Motivo <- factor(x= hd2DF$Motivo, labels = factor, levels = factor, ordered = T)

pocketDF$Motivo<- factor(x= pocketDF$Motivo, labels = factor, levels = factor, ordered = T)


RecallPlot <-  ggplot(hd2DF,
                    aes (x= Motivo, 
                         y = Frecuencia)) +
  geom_histogram(color= 'black',
                 fill = 'grey',
                 stat= 'identity')+
  geom_histogram(data = pocketDF, 
                 aes(y = Frecuencia),
                 color = 'black',
                 fill = pocketColor,
                 position= "stack",
                 stat= 'identity')+
  geom_text(data = pocketDF,
            aes(label = Frecuencia),
            vjust = -0.1,    #Para el plot Rb
            # vjust = -0.5,    #Para el plot p107
            size= 10) +
  geom_text(data = hd2DF,
            aes(label = Frecuencia),
            vjust = -0.5,
            size= 10) +
  scale_x_discrete(drop= FALSE) +
  scale_y_continuous(limits = c(0,14), 
                     breaks = seq(0, 12, 
                                  by= 2),
                     expand = c(0,0))+
  labs(#title= pocketTitle,
       x= "SLiM",
       y ="Frecuencia absoluta")+ tema_plots

ggsave(filename = fileOut ,
       plot = RecallPlot,
       device = "png",
       dpi = 300,
       width = 18,
       height = 16,
       units = "cm")
# return(RecallPlot)
}

# Calcular frecuencia de SLiMs en diferentes conjuntos de datos
hd2TPELM <- frecuencyTPELM(hd2File)
rbTPELM <- frecuencyTPELM(rbFile)
p107TPELM <- frecuencyTPELM(p107File)
p130TPELM <- frecuencyTPELM(p130File)

# Generar gráficos de comparación
p107Plot <- pocketRecallPlot(hd2DF = hd2TPELM, pocketDF = p107TPELM, factor = allMotifFactor, pocketColor = '#C75C70', fileOut = "ELMrecallp107.png")
p130Plot <- pocketRecallPlot(hd2DF = hd2TPELM, pocketDF = p130TPELM, factor = allMotifFactor, pocketColor = '#75a270ff', fileOut = "ELMrecallp130.png")

# Excluir el tercer motivo para el plot de Rb
hd2TPELM <- hd2TPELM[-3,]
rbTPELM <- rbTPELM[-3,]

rbPlot <- pocketRecallPlot(hd2DF = hd2TPELM,
                           pocketDF = rbTPELM,
                           factor = twoMotifFactor,
                           pocketColor = '#31538fff',
                          # pocketTitle = "Retinoblastoma",
                           fileOut = "ELMrecallRb.png")
