#SCRIPT DE ILUSTRACIÓN PARA EL TRABAJO 
#Castillo Moine, M. A., Balzarini, M. G.. 2019. Gestión de bases de datos espacio-temporales. 
#Agriscientia. 
#Descargar los datos de ilustración del siguiente enlace 

#https://agroubaar-my.sharepoint.com/:f:/g/personal/mcastillomoine_agro_uba_ar/Em7ca8y1BFBAkRWCm_ug7NIBZgiOUmepvrWQvKn0j2HTaA?e=JD5hCy

#Descomprimir la carpeta (por ejemplo en D:\\Google Drive\\Agrsc-GeDETIS) 

#Setear el directorio de trabajo 
setwd('D:\\Google Drive\\Agrsc-GeDETIS')
#Librerias y otros####
library(stringr)

library(bigmemory)
options(bigmemory.allow.dimnames=TRUE)
library(raster) 
rasterOptions(progress='text', timer=T, todisk=T, default=FALSE)

library(reshape2)
library(foreach)
library(doSNOW)

library(dtwclust)
library(kohonen)

library(tibbletime) # Future of tidy time series analysis
library(tidyquant)  # Loads tidyverse, tq_get()
library(reshape2)
library(tidyr)
library(cleangeo)
library(magrittr)
library(rgdal)
library(sf)
library(sp)
library(zoo)
library(sabre)
library(lwgeom)
library(tidyverse)
library(DataCombine)

#Especificar un directorio de salida para los procesamientos. 
dir_out='D:\\Google Drive\\Agrsc-GeDETIS' 

#Procesamiento de las bandas 1, 2 y atributos de calidad del producto MOD13Q1
#Vector de fechas####
folder=paste0(dir_out, '/MOD09Q1_Surf_Ref_8Days_250m_v6/QA_bits')
files=sort(grep(list.files(path=folder, all.files = FALSE, full.names = T), pattern = 'tif', inv=F, value=T)) #Se puede especificar otro patrón en caso de dobles coincidencias/dobles extensiones de archivos
Z=substr(files, 86, 93)
Z=gsub('_', '', Z)
Z=strptime(Z, format="%Y%j")
Z=as.Date(Z)
#Banda de calidad####
#El criterio de calidad es eliminar aquellos pixeles que no tengan calidad optima 
#Se realiza la extracción con MODISrts u otro paquete... 
folder='./MOD09Q1_Surf_Ref_8Days_250m_v6/QA_bits'
files=sort(grep(list.files(path=folder, all.files = FALSE, full.names = T), pattern = 'tif', inv=F, value=T)) #Se puede especificar otro patrón en caso de dobles coincidencias/dobles extensiones de archivos
qc=stack(files)
names(qc)<-Z
setZ(qc, Z, name='time')
#Crear un valor de máscara: 
#Se usa en todas las demás capas la siguiente expresión
#data=mask(data, qc, maskvalue=0, inverse=T)#Reemplaza por NA cualquier valor que no sea óptimo

####Variable 1: 250m Surface Reflectance Band 1 (620-670 nm)####
beginCluster()
folder='./MOD09Q1_Surf_Ref_8Days_250m_v6/b1_Red/'
files=sort(grep(list.files(path=folder, all.files = FALSE, full.names = T), pattern = 'tif', inv=F, value=T)) #Se puede especificar otro patrón en caso de dobles coincidencias/dobles extensiones de archivos
MOD09Q1_B1=stack(files)
names(MOD09Q1_B1)<-Z
setZ(MOD09Q1_B1, Z, name='time')

funcion <- function(x){ calc(x, function(x) { x[x>16000] <- NA; x[x<(-100)]<-NA; x[x==(-28672)]<-NA; return(x)})} #hecho
MOD09Q1_B1 <-clusterR(MOD09Q1_B1, fun=funcion, filename='MOD09Q1_B1', overwrite=T)#proceso lento, tener paciencia
MOD09Q1_B1=mask(MOD09Q1_B1, qc, maskvalue=0, inverse=T)
endCluster()

beginCluster()
nas=clusterR(MOD09Q1_B1, is.na)#identifica los NA en el stack; forma un nuevo stack
mm=clusterR(nas, fun=mean)#genera una máscara, donde las celdas que son siempre NA (ej. límite de la capa) valen 1
mask=clusterR(MOD09Q1_B1, fun=mask, args = list(maskvalue=1, updatevalue=0, updateNA=T), filename = masked)
MOD09Q1_B1=approxNA(mask)#Aparentemente no se puede aplicar cluster
MOD09Q1_B1=setZ(MOD09Q1_B1, Z)
names(MOD09Q1_B1)=Z
#Aplicar el factor de corrección
funcion <- function(x) x*0.0001
MOD09Q1_B1=clusterR(MOD09Q1_B1, calc, args=list(fun=funcion), filename = 'MOD09Q1_B1', overwrite=T) #Aplicar cluster en apariencia no mejora el tiempo de procesado... 
endCluster() #Cierra los cluster usados en este paso para evitar errores de sistema
removeTmpFiles(h=0.1)#Remueve los archivos termporales más antiguos que 0.1 horas

####Variable 2: 250m Surface Reflectance Band 2 (841-876 nm)####
folder='./MOD09Q1_Surf_Ref_8Days_250m_v6/b2_NIR/'
files=sort(grep(list.files(path=folder, all.files = FALSE, full.names = T), pattern = 'tif', inv=F, value=T)) #Se puede especificar otro patrón en caso de dobles coincidencias/dobles extensiones de archivos
MOD09Q1_B2=stack(files)
names(MOD09Q1_B2)<-Z
setZ(MOD09Q1_B2, Z, name='time')

beginCluster()
funcion <- function(x){ calc(x, function(x) { x[x>16000] <- NA; x[x<(-100)]<-NA;  x[x==(-28672)]<-NA; return(x)})} #hecho
MOD09Q1_B2 <-clusterR(MOD09Q1_B2, fun=funcion, filename='MOD09Q1_B2', overwrite=T)#hecho, lento pero hecho
MOD09Q1_B2=mask(MOD09Q1_B2, qc, maskvalue=0, inverse=T)
endCluster()

beginCluster()
nas=clusterR(MOD09Q1_B2, is.na)#identifica los NA en el stack; forma un nuevo stack
mm=clusterR(nas, fun=mean)#genera una máscara, donde las celdas que son siempre NA (ej. límite de la capa) valen 1
plot(mm)
#enmascarar todos los valores que son 1
mask=mask(MOD09Q1_B2, mm, maskvalue=1, updatevalue=-28672, updateNA=T)#-0,3 lo uso para indicar valores NA absurdos
#Ahora está apta para usar con approxNA
MOD09Q1_B2=approxNA(mask, filename = 'MOD09Q1_B2', overwrite=T)#Aparentemente no se puede aplicar cluster
names(MOD09Q1_B2)<-Z
setZ(MOD09Q1_B2, Z, name='time')

#Aplicar el factor de corrección
funcion <- function(x) x*0.0001
MOD09Q1_B2=clusterR(MOD09Q1_B2, calc, args=list(fun=funcion), filename = 'MOD09Q1_B2', overwrite=T) #Aplicar cluster en apariencia no mejora el tiempo de procesado... 
MOD09Q1_B2=calc(MOD09Q1_B2, function(x) { x=x*0.0001; return(x)}, filename='MOD09Q1_B2', overwrite=T)
endCluster() #Cierra los cluster usados en este paso para evitar errores de sistema
removeTmpFiles(h=0.1)#Remueve los archivos termporales más antiguos que 0.1 horas

#Cálculo del NDVI####
beginCluster()
NDVI=(MOD09Q1_B2-MOD09Q1_B1)/(MOD09Q1_B2+MOD09Q1_B1)
NDVI[NDVI<(-1)]=-1
NDVI[NDVI>1]=1
NDVI[is.na(NDVI)]=-1
endCluster()



#Capas auxiliares de información####

Mascara_Cordoba_soloborde <- readOGR(".\\SIG\\Mascara_Cordoba_soloborde.shp")
Mascara_Cordoba_soloborde=clgeo_Clean(Mascara_Cordoba_soloborde)
plot(Mascara_Cordoba_soloborde)

#En otra sesion de trabajo, se pueden volver a cargar los datos de entrada con 
#NDVI=stack('./NDVI')

#Para realizar la comparación de clasificaciones usaremos el 
#Mapa de Coberturas de Nivel 2 (IDECOR, 2018) recategorizado
m=raster('./SIG/Nivel2_Completo_31Agosto2018_segm_30m_EPSG_5345.tif')
n=m#genera una copia en memoria
n[n==1]=100
n[n==2]=100
n[n==3]=100
n[n==4]=200
n[n==5]=200
n[n==6]=300
n[n==7]=300
n[n==8]=300
n[n==9]=400
n[n==10]=400
n[n==11]=400
n[n==12]=500
n[n==13]=500
n[n==14]=500
n[n==15]=500
n[n==16]=600
n[n==17]=600
n[n==18]=700
n[n==19]=700
n[n==20]=800
n[n==21]=800
n[n==22]=900
plot(n)
writeRaster(n, '.\\SIG\\Nivel2_Completo_31Agosto2018_segm_30m_EPSG_5345_9cat.tif', 
            overwrite=T)
#En QGis se generó rápidamente una capa vectorial de ese mapa, recortada al área de estudio
l_cover_i=st_read('./SIG/Nivel2_Completo_31Agosto2018_segm_30m_EPSG_5345_9cat_mascara.shp')
st_crs(l_cover_i) = 5345#EPSG for the SRC of the layer
l_cover_i=st_make_valid(l_cover_i)
plot(l_cover_i)
#También se generó una capa raster, recortada al área de estudio
l_cover_r=raster('.\\SIG\\Nivel2_Completo_31Agosto2018_segm_30m_EPSG_5345_9cat_mascara.tif')
levelplot(l_cover_r)

#Ejemplo para generar una capa de máscara para zonas urbanas (no correr)
# z_u=raster('D:\\SIG\\Capitulo 1 y paper agriscientia\\Nivel2_Completo_31Agosto2018_segm_30m_EPSG_5345_parazonaurbana_masc.tif')
# z_u=rasterToPolygons(z_u, fun=function(x){x==500}, dissolve = T)
# writeOGR(z_u, dsn = dir_out, layer = 'zonas_urbanas_sin_inf_vial2', 
#          driver="ESRI Shapefile", overwrite_layer=T)

#La capa que usaremos de base se es la de Unidades de la Vegetación de la Argetina
#(Oyarzabal et al., 2018)
uv=readOGR('.\\UnidadesVegetacionArg_EPSG_5345_recorte.shp')
uv=clgeo_Clean(uv)
uv_i=intersect(uv, Mascara_Cordoba_soloborde)#formato OGR 
plot(uv_i)
uv_i_sf=st_as_sf(uv_i)


#Funciones para preprocesar las series de tiempo de cada pixel####
prepare_time_series_stack=function(main='Spatio-temporal classification', 
                                   stack=stack(NULL),
                                   dir_out= tempdir()){#En este ejemplo entra una lista de archivos
  require(signal)
  matriz=stack
  matriz=bigmemory::big.matrix(nrow = ncell(stack),
                               ncol = nlayers(stack),
                               backingfile = paste(main, '_raw_data', sep=''),
                               backingpath = dir_out)
  flush(matriz)
  matriz2=bigmemory::big.matrix(nrow = ncell(stack), 
                                ncol = ncell(stack), 
                                backingfile = paste(main, round(sum(runif(3))), sep=''),
                                backingpath = tempdir())
  pb=txtProgressBar(min = 0, max = nlayers(stack), initial = 0, style = 3)
  for(i in 1:nlayers(stack)){
    v=as.vector(stack[[i]])
    matriz[,i]=v
    setTxtProgressBar(pb, i)
  }
  print('Fin de armado de la matriz')
  i_i=seq(1, nrow(matriz), by=50000)
  i_f=c(i_i[2:length(i_i)]-1, nrow(matriz))
  pb=txtProgressBar(min = 0, max = length(i_i), initial = 0, style = 3)
  cl <- makeCluster()#Cuidado! La cantidad de núcleos que puede ocupar depende de la RAM disponible
  registerDoSNOW(cl)
  for(i in 1:length(i_i)){
    subset=sub.big.matrix(matriz, firstRow = i_i[i], lastRow = i_f[i])
    subset=as.matrix(subset)
    subset2=subset
    subset2[,]<-0
    #Aquí opera con la matriz subset y guardar ese resultado en la matriz original o una copia
    s=foreach(j=1:nrow(subset), .combine = 'rbind') %dopar% {
      #Formar la time serie y filtrarla
      time_serie=ts(subset[j,], frequency = 46)
      ts_sg=signal::sgolayfilt(time_serie, p = 3, n = 9, m = 0, ts = 1)
    }
    #Grabar cada parte en la big.matrix en disco 
    matriz2[i_i[i]:i_f[i],]=s
    setTxtProgressBar(pb, i)
  }
  flush(matriz2)
  deepcopy(x = matriz, backingfile = paste(main, '_raw_data', sep=''), backingpath = dir_out)
  deepcopy(x = matriz2, backingfile = paste(main, '_filtered_data', sep=''), backingpath = dir_out)
  stopCluster(cl)
  return(list('Raw data'=matriz, 'Filtered data'=matriz2))
}
data_NDVI_procesada=prepare_time_series_stack(main = 'NDVI', stack = NDVI_m, 
                                              dir_out = dir_out)
#Función que predice el grupo de pertenencia a un SOM a partir de un modelo SOM
predict_km_bigmem=function(big_matrix, model){
  require(DSL)
  require(clue)
  i_i=seq(1, nrow(big_matrix), by=50000)
  i_f=c(i_i[2:length(i_i)]-1, nrow(big_matrix))
  data_pred=list(NULL)
  data_pred=as.DList(data_pred)
  for(i in 1:length(i_i)){
    subset=sub.big.matrix(big_matrix, firstRow = i_i[i], lastRow = i_f[i])
    subset=as.matrix(subset)
    data_pred[[i]]=cl_predict(model, subset, type="class_ids")
  }
  return(unlist(data_pred))
}

#Grafico de dos series de tiempo####
o=as.vector(data_NDVI_procesada$`Raw data`[153548,])
f=as.vector(data_NDVI_procesada$`Filtered data`[153548,])
st=cbind(z=as.character(z), o=o, f=f)
#rownames(st)<-z
st=as.tibble(st)
st=mutate(st, z=as.Date(st$z), o=as.numeric(st$o), f=as.numeric(st$f))
ggplot()+
  scale_linetype_identity()+
  geom_line(data = st,
            aes(z, o), linetype='dashed', color='grey70')+
  geom_line(data = st,
            aes(z, f))+
  ylab('NDVI (desestac.)')+
  xlab('Fecha')+
  expand_limits(y=c(0,0.8))+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent", linetype=0, size=0) # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , axis.line = element_line(color = "grey80")
    , legend.title=element_blank()
  )

#Muestreo aleatorio####
#Para determinar el tamaño mínimo muestral de cada grupo se puede aplicar la fórmula de Cochran (en Oloffson et al. 1977)
c=0.975#percentil para construir el intervalo de confianza para prueba bilateral de hipótesis (1-[alfa/2])
O=0.5#precisión global esperada de encontrar una determinada clase para un pixel (valor de probabilidad que responde a la pregunta: ese pixel pertenece o no a esa clase de interés?)
d=0.025#Equivalente a alfa/2
n=round((qnorm(c)^2*O*(1-O))/d^2)

#A los fines de poder hacer más comparables las clasificaciones se usarán las mismas semillas
#Muestreo de las matrices 
set.seed(10)
#s=sample(1:nrow(data_NDVI_procesada$`Filtered data`), n)#muestrear una cantidad de series de tiempo lo suficientemente grande 
#saveRDS(s, paste(dir_out, 's', sep='/'))
s=readRDS(paste(dir_out, 's', sep = '/'))
data_train_matrix_originales <- as.matrix(data_NDVI_procesada$`Raw data`[s,]) #Crear un subset del RasterStack con las series de tiempo que serán usadas como conjunto de entrenamiento para la SOM.
data_train_matrix_filtradas <- as.matrix(data_NDVI_procesada$`Filtered data`[s,]) #Crear un subset de la BigMatrix con las series de tiempo que serán usadas como conjunto de entrenamiento para la SOM.

#Clasificación no supervisada mediante kmeans####
#Para decidir el número de grupos hay varias estrategias. Una posible considerando
#que las series de tiempo muestreadas al azar son independientes entre sí
#es aplicar la regla del codo. 
#(tomado de https://www.r-bloggers.com/finding-optimal-number-of-clusters/)
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(data_train_matrix_originales, 
                                 k, nstart=50,iter.max = 500 )$tot.withinss})
wss
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Numero de clusters K",
     ylab="Suma de cuadradados dentro de cluster total")


ggplot()+
  aes(x=1:k.max, y=wss)+
  geom_point()+
  geom_line()+
  xlab("Numero de clusters K")+
  ylab("Suma de cuadradados dentro de cluster total" )+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    #, panel.grid.major = element_blank() # get rid of major grid
    #, panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent", linetype=0, size=0) # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , axis.line = element_line(color = "grey80")
    , legend.title=element_blank()
  )

#(a) Generar un conjunto de entrenamiento muestreando las series de tiempo para n pixeles (output del paso anterior)
#(b) Generar un modelo de clasificación basado en el conjunto de entrenamiento con 8 grupos (un grupo es para los datos de máscara)
#(c) Usar el modelo para predecir la pertenecia de las series de tiempo de los pixeles restantes a un grupo

#(b) y (c) series de tiempo crudas####
#Por kmeans
pmt=proc.time()
model_petit_NDVI_originales_kmeans=kmeans(data_train_matrix_originales, 8, 
                                          nstart = 50, iter.max = 500)
#Si se usa el Rasterstack para predecir da error:
#predichos=predict(object=NDVI, model = model_petit_NDVI_originales_km)#Raster-based method

#El método anterior produce un error dado que solo admite como input matrices. En 
#su lugar usmos la función predict_km_bigmem()
predichos=predict_km_bigmem(data_NDVI_procesada$'Raw data', model_petit_NDVI_originales_kmeans)
proc.time()-pmt

ras=raster(NDVI)#Se puede usar como "molde" cualquier imagen raster con la misma extrensión y resolución
ras=setValues(ras, as.vector(predichos))
levelplot(ras)
sha=rasterToPolygons(ras, dissolve = T)
#plot(sha)
writeOGR(sha, dsn = dir_out, layer = 'model_NDVI_originales_kmeans_7G', 
         driver="ESRI Shapefile", overwrite_layer=T)
writeRaster(ras, paste(dir_out, 'model_NDVI_originales_kmeans_7G', sep='/'), format='GTiff', overwrite=T)

#(b) y (c) series de tiempo filtradas####
#Por kmeans
model_petit_NDVI_originales_kmeans=kmeans(data_train_matrix_filtradas, 8, nstart = 50, iter.max = 500)
predichos=predict_km_bigmem(NDVI_SG2, model_petit_NDVI_originales_kmeans)
proc.time()-pmt

ras=raster(NDVI)#Se puede usar como "molde" cualquier imagen raster con la misma extrensión y resolución
ras=setValues(ras, as.vector(predichos))
plot(ras, main='model_NDVI_filtradas_kmeans_7G_SG_masc')
sha=rasterToPolygons(ras, dissolve = T)
#plot(sha)
writeOGR(sha, dsn = dir_out, layer = 'model_NDVI_filtradas_kmeans_7G_SG_masc', 
         driver="ESRI Shapefile", overwrite_layer=T)
writeRaster(ras, paste(dir_out, 'model_NDVI_filtradas_kmeans_7G_SG_masc', sep='/'), 
            format='GTiff', overwrite=T)


#Ilustrando la función de filtrado de las series de tiempo con una serie de tiempo####
o=as.vector(data_NDVI_procesada$`Raw data`[153548,])
f=as.vector(data_NDVI_procesada$`Filtered data`[153548,])
st=cbind(z=as.character(z), o=o, f=f)
#rownames(st)<-z
st=as.tibble(st)
st=mutate(st, z=as.Date(st$z), o=as.numeric(st$o), f=as.numeric(st$f))
ggplot()+
  scale_linetype_identity()+
  geom_line(data = st,
            aes(z, o), linetype='dashed', color='grey70')+
  geom_line(data = st,
            aes(z, f))+
  ylab('NDVI (desestac.)')+
  xlab('Fecha')+
  expand_limits(y=c(0,0.8))+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent", linetype=0, size=0) # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , axis.line = element_line(color = "grey80")
    , legend.title=element_blank()
  )


#Gráficos y medidas de calidad de los mapas obtenidos####
#Muestreo de series de tiempo clasificadas para construir descripciones de los grupos formados####


#Mapa de las series de tiempo sin depurar####
mapa_crudas_r=raster('.\\model_NDVI_originales_kmeans_7G.tif')
mapa_crudas_r=mask(mapa_crudas_r, Mascara_Cordoba_soloborde)
mapa_crudas_s=readOGR('.\\model_NDVI_originales_kmeans_7G.shp')
mapa_crudas_s=cleangeo::clgeo_Clean(mapa_crudas_s)
mapa_crudas_i=intersect(mapa_crudas_s, Mascara_Cordoba_soloborde)

#Hacer el muestreo estratificado 
#Para determinar el tamaño de muestra usaremos la formula 13 de Oloffson 2014 (Good practices for estimating area and assessing accuracy of land change)
#usando un S(O) (error estandar de la precisión general) de 0.01  y considerando una 
#precisión de usuario de 0.50 (esto es, la probabilidad de acertarle a la verdadera 
#clase es azarosa) por lo que Si=sqrt(0.5*(1-0.5))=0.5 (esto es posible solo bajo la sufecha de que nuestra clasificacion en cierta)
a=area(mapa_crudas_i)
a=sapply(a, function(x) x/sum(a))*0.5
n=(sum(a)/0.01)^2

mapa_crudas_r_sample=sampleStratified(mapa_crudas_r, n, xy=T, sp=T)

series_crudas_sample=data_NDVI_procesada$`Raw data`[mapa_crudas_r_sample$cell,]
series_crudas_sample=as.tibble(series_crudas_sample)
dim(series_crudas_sample)
names(series_crudas_sample)=as.character(z)
series_crudas_sample=mutate(series_crudas_sample, grupo=mapa_crudas_r_sample$model_NDVI_originales_kmeans_7G)
series_crudas_sample.large <- series_crudas_sample %>%
  gather(fecha, valor, `2000-02-18`:`2018-02-18`)
series_crudas_sample.large=mutate(series_crudas_sample.large, 
                                  grupo=as.factor(series_crudas_sample.large$grupo), 
                                  fecha=as.Date(series_crudas_sample.large$fecha))
series_crudas_sample.large=series_crudas_sample.large %>% 
  group_by(grupo, fecha) %>%
  summarise(media=mean(valor), 
            q1=quantile(valor, 0.25), mediana=median(valor), q2=quantile(valor, 0.75))
#Grafica todas las series juntas:
ggplot()+
  geom_line(data = series_crudas_sample.large, 
            mapping = aes(x = fecha, y = media, color=grupo))+
  geom_ribbon(data = series_crudas_sample.large, 
              aes(x=fecha, 
                  ymin=q1, 
                  ymax=q2,
                  fill = grupo), alpha = 0.3, linetype=0)+
  ggtitle('crudas 7G')+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

#grafico por grupo individual, ir sumando una tabla a rbind por grupo que se quiere graficar
series_crudas_sample.large_split=split(series_crudas_sample.large, series_crudas_sample.large$grupo)
#Si a las tablas que queremos graficar las re-combinamos en una sola obtenemos las leyendas de  manera automática
series.combinadas=rbind(series_crudas_sample.large_split$`4`, 
                        series_crudas_sample.large_split$`6`, 
                        series_crudas_sample.large_split$`7`)
ggplot() +
  geom_ribbon(data = series.combinadas,
              aes(x=fecha,
                  ymin=q1,
                  ymax=q2,
                  fill=grupo), alpha = 1, linetype=0)+
  scale_fill_manual(values = c('#8fb775', '#dbcf72', '#e7c5c5'))+
  geom_line(data = series.combinadas,
            aes(fecha, mediana, group=grupo)) +
  scale_color_manual(values = c('#8fb775', '#dbcf72', '#e7c5c5'))+
  ylab('NDVI')+
  xlab('')+
  expand_limits(y=c(0,0.8))+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent", linetype=0, size=0) # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , axis.line.y = element_line(color = "grey80")
    , axis.line.x = element_blank()
    , axis.text.x = element_blank()
    , axis.ticks.x = element_blank()
    , legend.title=element_blank()
  )

#
#Recortar con la capa de borde de Córdoba... O leer el objeto que ya está cortado y transformarlo a sf
NDVI.c_sf_masked=st_as_sf(mapa_crudas_i)
NDVI.c_sf_masked=st_make_valid(NDVI.c_sf_masked)
st_crs(NDVI.c_sf_masked) = 5345#EPSG for the SRC of the layer

#Entonces ahpra lo que hay que comparar es NDVI.f_sf_masked con landsat_comp
mcu=mapcurves_calc(x=NDVI.c_sf_masked, x_name = layer, 
                   y=l_cover_i, y_name=DN)


v_m=vmeasure_calc(x=NDVI.c_sf_masked, x_name = layer, 
                  y=l_cover_i, y_name=DN)

st_write(v_m$map1, dsn=paste(dir_out, 'NDVI_crudas_kmeans_vm_l_cover.shp', sep='/'))
st_write(v_m$map2, dsn=paste(dir_out, 'NDVI_crudas_kmeans_vm_2_l_cover.shp', sep='/'))

new=c(Mapa='NDVI_crudas_kmeans_vm_l_cover.shp', mcu_ref_map=mcu$ref_map, mcu_gof=mcu$gof, v_m_v_measure=v_m$v_measure, 
      v_m_homogeneity=v_m$homogeneity, v_m_completeness=v_m$completeness)
results=InsertRow(results, new)

#Mapa de las series de tiempo filtradas####
mapa_filtrada7_r=raster('D:\\SIG\\Capitulo 1 y paper agriscientia\\model_NDVI_filtradas_kmeans_7G_SG_masc.tif')
mapa_filtrada7_r=mask(mapa_filtrada7_r, Mascara_Cordoba_soloborde)
mapa_filtrada7_s=readOGR('D:\\SIG\\Capitulo 1 y paper agriscientia\\model_NDVI_filtradas_kmeans_7G_SG_masc.shp')
mapa_filtrada7_s=cleangeo::clgeo_Clean(mapa_filtrada7_s)
mapa_filtrada7_i=intersect(mapa_filtrada7_s, Mascara_Cordoba_soloborde)
plot(mapa_filtrada7_i)

#GRáficos de las series de tiempo
a=area(mapa_filtrada7_i)
a=sapply(a, function(x) x/sum(a))*0.5
n=(sum(a)/0.01)^2
mapa_filtradas_r_sample=sampleStratified(mapa_filtrada7_r, n, xy=T, sp=T)

series_filtradas_sample=data_NDVI_procesada$`Filtered data`[mapa_filtradas_r_sample$cell,]
series_filtradas_sample=as.tibble(series_filtradas_sample)
dim(series_filtradas_sample)
names(series_filtradas_sample)=as.character(z)

series_filtradas_sample=mutate(series_filtradas_sample, 
                               grupo=
                                 mapa_filtradas_r_sample$model_NDVI_originales_kmeans_7G_SG_masc)
series_filtradas_sample.large <- series_filtradas_sample %>%
  gather(fecha, valor, `2000-02-18`:`2018-02-18`)
series_filtradas_sample.large=mutate(series_filtradas_sample.large, 
                                     grupo=as.factor(series_filtradas_sample.large$grupo), 
                                     fecha=as.Date(series_filtradas_sample.large$fecha))


series_filtradas_sample.large=series_filtradas_sample.large %>% 
  group_by(grupo, fecha) %>%
  summarise(media=mean(valor), 
            q1=quantile(valor, 0.25), mediana=median(valor), q2=quantile(valor, 0.75))

#grafico por grupo individual filtradas#### 
series_filtradas_sample.large_split=split(series_filtradas_sample.large, series_filtradas_sample.large$grupo)
ggplot() +
  geom_line(data = series_filtradas_sample.large_split$`8`,
            aes(fecha, mediana)) +
  geom_ribbon(data = series_filtradas_sample.large_split$`8`,
              aes(x=fecha,
                  ymin=q1,
                  ymax=q2), fill='grey70', alpha = 0.5, linetype=0)+
  ggtitle('Grupo 7')+
  ylab('NDVI')+
  xlab('Fecha')+
  expand_limits(y=c(0,0.8))+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , axis.line = element_line(color = "grey80")
  )


ggplot() +
  geom_line(data = series_filtradas_sample.large_split$`1`,
            aes(fecha, mediana), color='#8fb775') +
  geom_line(data = series_filtradas_sample.large_split$`2`,
            aes(fecha, mediana), color='#dbcf72') +
  geom_line(data = series_filtradas_sample.large_split$`3`,
            aes(fecha, mediana), color='#e7c5c5') +
  geom_ribbon(data = series_filtradas_sample.large_split$`1`,
              aes(x=fecha,
                  ymin=q1,
                  ymax=q2), fill='#8fb775', alpha = 0.5, linetype=0)+
  geom_ribbon(data = series_filtradas_sample.large_split$`2`,
              aes(x=fecha,
                  ymin=q1,
                  ymax=q2), fill='#dbcf72', alpha = 0.5, linetype=0)+
  geom_ribbon(data = series_filtradas_sample.large_split$`3`,
              aes(x=fecha,
                  ymin=q1,
                  ymax=q2), fill='#e7c5c5', alpha = 0.5, linetype=0)+
  ggtitle('filtradas ns cleaned kmeans G 1-2-3')+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent") # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  )

#Si a las tablas las combinamos en una sola obtenemos las leyendas
series.combinadas=rbind(series_filtradas_sample.large_split$`1`, 
                        series_filtradas_sample.large_split$`2`, 
                        series_filtradas_sample.large_split$`3`)
ggplot() +
  geom_ribbon(data = series.combinadas,
              aes(x=fecha,
                  ymin=q1,
                  ymax=q2,
                  fill=grupo), alpha = 1, linetype=0)+
  scale_fill_manual(values = c('#8fb775', '#dbcf72', '#e7c5c5'))+
  geom_line(data = series.combinadas,
            aes(fecha, mediana, group=grupo)) +
  scale_color_manual(values = c('#8fb775', '#dbcf72', '#e7c5c5'))+
  ylab('NDVI (desestac.)')+
  xlab('Fecha')+
  expand_limits(y=c(0,0.8))+
  theme(
    panel.background = element_rect(fill = "transparent") # bg of the panel
    , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
    , panel.grid.major = element_blank() # get rid of major grid
    , panel.grid.minor = element_blank() # get rid of minor grid
    , legend.background = element_rect(fill = "transparent", linetype=0, size=0) # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    , axis.line = element_line(color = "grey80")
    , legend.title=element_blank()
  )


#Convertir a sf
mapa_filtrada7_st=st_read('D:\\SIG\\Capitulo 1 y paper agriscientia\\model_NDVI_originales_kmeans_7G_SG_masc.shp')
Mascara_Cordoba_soloborde_sf=st_as_sf(Mascara_Cordoba_soloborde)
st_crs(mapa_filtrada7_st) = 5345#EPSG for the SRC of the layer
mapa_filtrada7_st=st_make_valid(mapa_filtrada7_st)
mapa_filtrada7_st_i=st_intersection(mapa_filtrada7_st, Mascara_Cordoba_soloborde_sf)
plot(mapa_filtrada7_st_i)

#Entonces ahora lo que hay que comparar es NDVI.f_sf_masked con uv_i
mcu=mapcurves_calc(x=mapa_filtrada7_st_i, x_name = layer, 
                   y=uv_i_sf, y_name=NOMFISONOM)




v_m=vmeasure_calc(x=mapa_filtrada7_st_i, x_name = layer, 
                  y=uv_i_sf, y_name=NOMFISONOM)

st_write(v_m$map1, dsn=paste(dir_out, 'NDVI_filtrada7_kmeans_vm_uvi.shp', sep='/'))
st_write(v_m$map2, dsn=paste(dir_out, 'NDVI_filtrada7_kmeans_vm_2_uvi.shp', sep='/'))

new=c(Mapa='NDVI_crudas6_kmeans_vm_uvi.shp', mcu_ref_map=mcu$ref_map, mcu_gof=mcu$gof, v_m_v_measure=v_m$v_measure, 
      v_m_homogeneity=v_m$homogeneity, v_m_completeness=v_m$completeness)
results=InsertRow(results, new)
