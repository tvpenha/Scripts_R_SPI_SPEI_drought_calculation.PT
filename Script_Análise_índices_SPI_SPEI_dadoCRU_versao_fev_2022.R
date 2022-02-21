# -------------------------------------------------------------------------	#
#
#
#
#					SCRIPT PARA ANÁLISE DE SECA UTILIZANDO OS ÍNDICES SPI e SPEI 							
#
#
#
#
# Criado:		28-09-2021						
# Modificado:	21-02-2022		  
#																			
# Autor:															
#				Thales Vaz Penha - thales.penha@usp.br 					
#				
#																			
# Software: R version 4.0.4 (2021-02-15)
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# -------------------------------------------------------------------------	#
#
# Dado de entrada: 
#
# SPEI-CRU CRU TS 4.05 dataset
#
# Fonte: Cimate Research Unit (CRU)
# Resolução espacial: 0.5°
# Janela temporal do dado: 1901-2020
# Formato: NetCDF
# Referências:
# https://spei.csic.es/database.html
# https://spei.csic.es/spei_database/#map_name=spei01#map_position=1415
# https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/
# https://zenodo.org/record/834462#.YOzJzOhKg2w
# VICENTE-SERRANO, Sergio M. et al. A new global 0.5 gridded dataset (1901-2006)of a multiscalar drought index: comparison with current drought index datasets based on the Palmer Drought Severity Index. Journal of Hydrometeorology, v. 11, n. 4, p. 1033-1043, 2010.
# BEGUERÍA, Santiago; VICENTE-SERRANO, Sergio M.; ANGULO-MARTÍNEZ, Marta. A multiscalar global drought dataset: the SPEIbase: a new gridded product for the analysis of drought variability and impacts. Bulletin of the American Meteorological Society, v. 91, n. 10, p. 1351-1354, 2010.
#
# -------------------------------------------------------------------------	#

############################################################################################################


rm(list=ls()) # remove da memória global os dados anteriormente trabalhados

memory.limit (9999999999) # expande o limite de memória de armazenamento e processamento de dados no R


setwd("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/")


############################################################################################################


# ----- Bibliotecas exigidas ------ # caso não tenha as bibliotecas instaladas, utilize o comando "install.packages" no lugar de "library" e depois utilize o comando "library"

library('raster') ; library('rgdal') ; library('ncdf4') ; library('utils')
library('sp') ; library('RNetCDF') ; library('SPEI') ; library('rasterVis')


############################################################################################################


# ----- Definir área de estudo
# 
# Limites da região Northeastern South America (NES) do IPCC-WGI
# Fonte: https://github.com/SantanderMetGroup/ATLAS/tree/v1.6/reference-regions
# Referência: https://essd.copernicus.org/articles/6/2959/2020/essd-6-2959-2020.html
# 
# Limites do Semiárido Brasileiro
# Fonte: http://antigo.sudene.gov.br/delimitacao-do-semiarido
#
# -------------------------------------------------------------------------	#


# Opção 1: Abrir shapefile com os limites da área de estudo
NES = readOGR("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/Shapefile/IPCC-WGI-NES_region_v4.shp") # Região NES do IPCC

SAB = readOGR("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/Shapefile/LIM_Semiarido_Brasileiro.shp") # limites do Semiárido Brasileiro

# Opção 2: Definir o quadrangular que engloba a área de estudo | aqui corresponde aos limites da região NES
coords = matrix(c(-34.0, -20.0,
                  -50.0, -20.0,
                  -50.0, 0.0,
                  -34.0, 0.0), 
                  ncol = 2, byrow = TRUE)

NES = Polygon(coords)
NES = SpatialPolygons(list(Polygons(list(NES), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(NES, axes = TRUE)


###############################################################################################################################
#
#
# CONHECENDO ASPECTOS BÁSICOS DOS ÍNDICES SPI e SPEI 
#
#
###############################################################################################################################


# Abrir o arquivo de entrada NetCDF SPI - 6 meses

spi6 <- brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/spi6_cru_1990_2019.nc")

# Visualizar as características do arquivo de entrada NetCDF SPI
print(spi6)


#-----------------------------------------------------------------------------------------------------------#


# Abrir o arquivo de entrada NetCDF SPEI - 6 meses

spei6 <- brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/spei6_cru_1990_2019.nc")

# Visualizar as características do arquivo de entrada NetCDF SPEI
print(spei6)


###############################################################################################################################
#
#
# VISUALIZAÇÃO EM MAPA DO SPEI PARA A REGIÃO NES EM UM ANO ESPECÍFICO
#
#
###############################################################################################################################


# ----- Apresentação em mapa 

# recorte temproal para o ano de 2019

spi6_2019 = subset(spi6, 349:360) #recorte para 2019


spei6_2019 = subset(spei6, 349:360) #recorte para 2019

# atualzação dos nomes de cada mês

names(spi6_2019) <- c("Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")


names(spei6_2019) <- c("Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")

# estatística geral

summary(spi6_2019) # estatística básica


summary(spei6_2019) # estatística básica

# visualização dos mapas

cuts=c(2.0,1.5,1.0,0.5,0,-0.5,-1.0,-1.5,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow","white","light blue","blue","dark blue")) #define as cores da legenda

# Mapa SPI
plot(spi6_2019, breaks=cuts, col = pal(9)) #gera o mapa com as cores e legendas pré-definidas


# Mapa SPEI
plot(spei6_2019, breaks=cuts, col = pal(9)) #gera o mapa com as cores e legendas pré-definidas


# Exporta raster como tif 

writeRaster(spi6_2019,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPEI_6_seca_2019.tif')


writeRaster(spei6_2019,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPEI_6_seca_2019.tif')

###############################################################################################################################
#
#
# GRÁFICO SPEI DE SÉRIE TEMPORAL PARA A REGIÃO NES PARA O PERÍODO DE 1990-2019
#
#
###############################################################################################################################


# ----- Análise em gráfico

#transformação do dado em dataframe

spi6_df = as.data.frame(spi6)
head(spi6_df)


spei6_df = as.data.frame(spei6)
head(spei6_df)


#definindo as datas da sequencia dos dados

dates <- seq(as.Date("1990/1/1"), by = "month", length.out = 360)
head(dates)
tail(dates)


# cálculo de média de todos os pontos de grade para gerar os gráficos

spi6_df_mean = as.data.frame(colMeans(spi6_df, na.rm = T))


spei6_df_mean = as.data.frame(colMeans(spei6_df, na.rm = T))


# substitui valores infinitos por NA

spi6_df_mean=do.call(data.frame,lapply
              (spi6_df_mean, function(value) replace (value, is.infinite(value),NA)))


spei6_df_mean=do.call(data.frame,lapply
                       (spei6_df_mean, function(value) replace (value, is.infinite(value),NA)))

# renomeando colunas

names(spi6_df_mean) <- "spi_6"


names(spei6_df_mean) <- "spei_6"

# Adicionando as datas

spi6_df_mean$Dates <-seq(as.Date("1990/1/1"), by = "month", length.out = 360)


spei6_df_mean$Dates <-seq(as.Date("1990/1/1"), by = "month", length.out = 360)


# ----- VISUALIZAÇÃO EM GRAFICO - SPI


library(ggplot2)
library(scales)
# Basic barplot
ggplot(data=spi6_df_mean, aes(x=Dates, y=spi_6)) +
  geom_bar(aes(fill = spi_6 < 0), stat = "identity") + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("red", "blue"))+
  geom_hline(yintercept=-0.5, linetype="dashed", color = "black")+
  labs(x = "Month/Year", y = "SPI", title = "Região NES SPI - escala 6 meses - Série temporal 1990-2019") +
  scale_y_continuous(limits = c(-2.0, 2.0))+
  scale_x_date(date_breaks = "2 year",
               labels = date_format("%Y"))




# ----- VISUALIZAÇÃO EM GRAFICO - SPEI

library(ggplot2)
library(scales)
# Basic barplot
ggplot(data=spei6_df_mean, aes(x=Dates, y=spei_6)) +
  geom_bar(aes(fill = spei_6 < 0), stat = "identity") + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("red", "blue"))+
  geom_hline(yintercept=-0.5, linetype="dashed", color = "black")+
  labs(x = "Month/Year", y = "SPEI", title = "Região NES SPEI - escala 6 meses - Série temporal 1990-2019") +
  scale_y_continuous(limits = c(-2.0, 2.0))+
  scale_x_date(date_breaks = "2 year",
               labels = date_format("%Y"))




# ----- Teste estatístico para análise de tendência na série temporal

# Utilizaremos o teste não paramétrico de Mann-Kendall
# H0 (hipótese nula) = não há tendência no dado
# H1 (hipótese alternativa) = há tendência no dado

library(Kendall)
# performar o teste de Mann-Kendall

#SPI
MK_test_spi = MannKendall(as.vector(spi6_df_mean$spi_6)) # considerar como significante estatítsticamente quando o p-valor equivale a 0.05 

# SPEI
MK_test_spei = MannKendall(as.vector(spei6_df_mean$spei_6)) # considerar como significante estatítsticamente quando o p-valor equivale a 0.05 


# ----- VISUALIZAÇÃO EM GRAFICO com linha de tendência (linear) - SPI

library(ggplot2)
library(scales)
# Basic barplot
ggplot(data=spi6_df_mean, aes(x=Dates, y=spi_6)) +
  geom_bar(aes(fill = spi_6 < 0), stat = "identity") + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("red", "blue"))+
  labs(x = "Month/Year", y = "SPI", title = "Curva de tendência para Região NES SPI - escala 6 meses - Série temporal 1990-2019") +
  geom_smooth(method = 'lm') +
  scale_y_continuous(limits = c(-2.0, 2.0))+
  scale_x_date(date_breaks = "2 year",
               labels = date_format("%Y"))





# ----- VISUALIZAÇÃO EM GRAFICO com linha de tendência (linear) - SPEI

library(ggplot2)
library(scales)
# Basic barplot
ggplot(data=spei6_df_mean, aes(x=Dates, y=spei_6)) +
  geom_bar(aes(fill = spei_6 < 0), stat = "identity") + scale_fill_manual(guide = FALSE, breaks = c(TRUE, FALSE), values=c("red", "blue"))+
  labs(x = "Month/Year", y = "SPEI", title = "Curva de tendência para Região NES SPEI - escala 6 meses - Série temporal 1990-2019") +
  geom_smooth(method = "lm")+ 
  scale_y_continuous(limits = c(-2.0, 2.0))+
  scale_x_date(date_breaks = "2 year",
               labels = date_format("%Y"))



###############################################################################################################################
#
#
# ANÁLISE DE SECA COM SPI PARA A REGIÃO NES NA SÉRIE HISTÓRICA POR MÊS DO ANO
#
#
###############################################################################################################################

# Definindo eventos de seca na série (1990-2019) quando o limiar de SPI < -0.5

spi6_dry = spi6

spi6_dry[spi6_dry[] > -0.5] <- NA # atribui valores nulos aos valores superiores ao limiar de seca do SPI, ou seja, representam valores normais ou úmidos


# gerar o valor médio de seca (SPI < - 0.5) por mês do ano
spi6_dry_month <- stackApply(spi6_dry, c( rep(seq(1:6), times = 30)), mean, na.rm = T ) 

# atualzação dos nomes de cada mês
names(spi6_dry_month) <- c("Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")

# estatística geral
summary(spi6_dry_month) 

# visualização dos mapas
cuts=c(0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow")) #define as cores da legenda

plot(spi6_dry_month, breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas


# Visualização em gráfico de barras

# transformação para data frame

spi6_dry_month_df = as.data.frame(spi6_dry_month)

# cálculo de média de todos os pontos de grade para gerar os gráficos

spi6_dry_month_df = as.data.frame(colMeans(spi6_dry_month_df, na.rm = T))

# renomeando colunas

names(spi6_dry_month_df) <- "spi_6"

#definindo as datas da sequencia dos dados

meses <- c("Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")

# Adicionando as meses do ano

spi6_dry_month_df$meses <- meses

# Visualização gráfica
ggplot(spi6_dry_month_df, aes(y = spi_6, x = meses, fill = meses)) +
  geom_bar(stat = "identity")



# Visualização mais elegante dos mapas 

library(rasterVis)
library(ggmap)

# Tranformacao dos limites da área de estudo para visualização
SAB_plot = fortify(SAB) #transforma o shapefile em dataframe para o gráfico


gplot(spi6_dry_month) + 
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient(low = "red", high = "yellow", na.value = NA) +
  coord_equal()+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))+
  geom_polygon(data = SAB_plot, aes(x = long, y = lat, group = group), colour = alpha("black", 1/2), fill = NA)+
    labs(x = "Latitude", y = "Longitude", fill="spi",
       title = "Região do Semiárido Brasileiro - SPI - escala 6 meses - Série temporal 1990-2019")
  

#plot(spi6_dry_month)


# Exporta raster como tif 
writeRaster(spi6_dry_month,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPI_6_seca_Mes_1990-2019.tif')


###############################################################################################################################
#
#
# ANÁLISE DE SECA COM SPEI PARA A REGIÃO NES POR SECA RECORRENTE
#
#
###############################################################################################################################

# recorte temporal para os anos de 2010-2019

spei6_2010 = subset(spei6, 241:252) #recorte para 2010
spei6_2011 = subset(spei6, 253:264) #recorte para 2011
spei6_2012 = subset(spei6, 265:276) #recorte para 2012
spei6_2013 = subset(spei6, 277:288) #recorte para 2013
spei6_2014 = subset(spei6, 289:300) #recorte para 2014
spei6_2015 = subset(spei6, 301:312) #recorte para 2015
spei6_2016 = subset(spei6, 313:324) #recorte para 2016
spei6_2017 = subset(spei6, 325:336) #recorte para 2017
spei6_2018 = subset(spei6, 337:348) #recorte para 2018
spei6_2019 = subset(spei6, 349:360) #recorte para 2019

# Definindo eventos de seca na série (2010-2019) quando o limiar de spei < -0.5
# reclassificando para valores binário (valores de seca = 1)

reclass <- cbind(from = c(-Inf,-0.5), to = c(-0.5, Inf), becomes = c(1, 0))

spei6_2010 <- reclassify(spei6_2010, reclass)
spei6_2011 <- reclassify(spei6_2011, reclass)
spei6_2012 <- reclassify(spei6_2012, reclass)
spei6_2013 <- reclassify(spei6_2013, reclass)
spei6_2014 <- reclassify(spei6_2014, reclass)
spei6_2015 <- reclassify(spei6_2015, reclass)
spei6_2016 <- reclassify(spei6_2016, reclass)
spei6_2017 <- reclassify(spei6_2017, reclass)
spei6_2018 <- reclassify(spei6_2018, reclass)
spei6_2019 <- reclassify(spei6_2019, reclass)

# contagem de recorrência da seca entre 2010-2019

spei6_dry_2010_2019 = (spei6_2010 + spei6_2011 + spei6_2012 + spei6_2013 + spei6_2014 +
                        spei6_2015 + spei6_2016 + spei6_2017 + spei6_2018 + spei6_2019)

# atualização dos nomes de cada mês
names(spei6_dry_2010_2019) <- c("Janeiro", "Fevereiro", "Marco", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")


# Mapas de recorrência de seca 
# (quantidade de vezes que o pixel apresentou spei < -0.5 entre 2010-2019)
# visualização dos mapas
cuts=c(0,1,2,3,4,5,6,7,8) #define os intervalos da legenda
pal <- colorRampPalette(c("grey","yellow","orange","red")) #define as cores da legenda

plot(spei6_dry_2010_2019, breaks=cuts, col = pal(8)) #gera o mapa com as cores e legendas pré-definidas


# Mapas mais elegantes

gplot(spei6_dry_2010_2019) + 
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint= 4, na.value = NA) +
  coord_equal()+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))+
  geom_polygon(data = SAB_plot, aes(x = long, y = lat, group = group), colour = alpha("black", 1/2), fill = NA)+
  labs(x = "Latitude", y = "Longitude", fill="Anos",
       title = "Região do Semiárido Brasileiro - recorrência de seca - spei - Série temporal 2010-2019")



# Exporta raster como tif 
writeRaster(spei6_dry_2010_2019,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_spei_6_seca_2010-2019.tif')


###############################################################################################################################
#
#
# ANÁLISE DE SECA COM SPI PARA A REGIÃO NES NA SÉRIE HISTÓRICA POR ESTAÇÃO DO ANO
#
#
###############################################################################################################################


#---------- Mapa - Verão


NES_summer_2015 <- subset(spi6, 310:315) # Out-Mar para 2015/2016
NES_summer_2016 <- subset(spi6, 322:327) # Out-Mar para 2016/2017
NES_summer_2017 <- subset(spi6, 334:339) # Out-Mar para 2017/2018
NES_summer_2018 <- subset(spi6, 346:351) # Out-Mar para 2018/2019


NES_summer_2015_mean <- calc(NES_summer_2015, fun = mean)
NES_summer_2016_mean <- calc(NES_summer_2016, fun = mean)
NES_summer_2017_mean <- calc(NES_summer_2017, fun = mean)
NES_summer_2018_mean <- calc(NES_summer_2018, fun = mean)

# Somente eventos de seca
NES_summer_2015_mean[NES_summer_2015_mean[] > -0.5 ] = NA
NES_summer_2016_mean[NES_summer_2016_mean[] > -0.5 ] = NA
NES_summer_2017_mean[NES_summer_2017_mean[] > -0.5 ] = NA
NES_summer_2018_mean[NES_summer_2018_mean[] > -0.5 ] = NA

NES_summer_2015_2019_map = stack(NES_summer_2015_mean, NES_summer_2016_mean, NES_summer_2017_mean, NES_summer_2018_mean)  

names(NES_summer_2015_2019_map)<- c("2015-2016 Summer (Oct-Mar)", "2016-2017 Summer(Oct-Mar)", "2017-2018 Summer(Oct-Mar)","2018-2019 Summer(Oct-Mar)")

# Visualização do mapa

cuts=c(0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow")) #define as cores da legenda

plot(NES_summer_2015_2019_map, breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas



# Visualização dos mapas individualmente sobre o Semiárido brasileiro
plot(NES_summer_2015_2019_map$X2015.2016.Summer..Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_summer_2015_2019_map$X2016.2017.Summer.Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_summer_2015_2019_map$X2017.2018.Summer.Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_summer_2015_2019_map$X2018.2019.Summer.Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)


# Exporta raster como tif 
writeRaster(NES_summer_2015_2019_map,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_spi_dry_summer_2015-2019.tif')


######################################################################################################################################################


#---------- Mapa - Inverno


NES_winter_2015 <- subset(spi6, 304:309) # abr-set para 2015
NES_winter_2016 <- subset(spi6, 316:321) # abr-set para 2016
NES_winter_2017 <- subset(spi6, 328:333) # abr-set para 2017
NES_winter_2018 <- subset(spi6, 340:345) # abr-set para 2018
NES_winter_2019 <- subset(spi6, 352:357) # abr-set para 2019

NES_winter_2015_mean <- calc(NES_winter_2015, fun = mean)
NES_winter_2016_mean <- calc(NES_winter_2016, fun = mean)
NES_winter_2017_mean <- calc(NES_winter_2017, fun = mean)
NES_winter_2018_mean <- calc(NES_winter_2018, fun = mean)
NES_winter_2019_mean <- calc(NES_winter_2019, fun = mean)

# Somente eventos de seca
NES_winter_2015_mean[NES_winter_2015_mean[] > -0.5 ] = NA
NES_winter_2016_mean[NES_winter_2016_mean[] > -0.5 ] = NA
NES_winter_2017_mean[NES_winter_2017_mean[] > -0.5 ] = NA
NES_winter_2018_mean[NES_winter_2018_mean[] > -0.5 ] = NA
NES_winter_2019_mean[NES_winter_2019_mean[] > -0.5 ] = NA

NES_winter_2015_2019_map = stack(NES_winter_2015_mean, NES_winter_2016_mean, NES_winter_2017_mean, NES_winter_2018_mean, NES_winter_2019_mean)  

names(NES_winter_2015_2019_map)<- c("2015 winter (abr-set)", "2016 winter (abr-set)", "2017 winter (abr-set)","2018 winter (abr-set)", "2019 winter (abr-set)")

# Visualização do mapa

cuts=c(0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow")) #define as cores da legenda

plot(NES_winter_2015_2019_map, breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas



# Visualização dos mapas individualmente sobre o Semiárido brasileiro
plot(NES_winter_2015_2019_map$X2015.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2016.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2017.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2018.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2019.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)


# Exporta raster como tif 
writeRaster(NES_winter_2015_2019_map,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_spi_dry_winter_2015-2019.tif')





###############################################################################################################################
#
#
# ANÁLISE DE SECA COM SPEI PARA A REGIÃO NES NA SÉRIE HISTÓRICA POR MÊS DO ANO
#
#
###############################################################################################################################

# Definindo eventos de seca na série (1990-2019) quando o limiar de SPEI < -0.5

spei6_dry = spei6

spei6_dry[spei6_dry[] > -0.5] <- NA # atribui valores nulos aos valores superiores ao limiar de seca do SPEI, ou seja, representam valores normais ou úmidos


# gerar o valor médio de seca (SPEI < - 0.5) por mês do ano
spei6_dry_month <- stackApply(spei6_dry, c( rep(seq(1:6), times = 30)), mean, na.rm = T ) 

# atualzação dos nomes de cada mês
names(spei6_dry_month) <- c("Janeiro", "Fevereiro", "Março", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")

# estatística geral
summary(spei6_dry_month) 

# visualização dos mapas
cuts=c(0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow")) #define as cores da legenda

plot(spei6_dry_month, breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas


# Visualização dos mapas mais elegantes

library(rasterVis)
library(ggmap)

# Transformação dos limites da áres de estudo para visualização
SAB_plot = fortify(SAB) #transforma o shapefile em dataframe para o gráfico


gplot(spei6_dry_month) + 
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient(low = "red", high = "yellow", na.value = NA) +
  coord_equal()+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))+
  geom_polygon(data = SAB_plot, aes(x = long, y = lat, group = group), colour = alpha("black", 1/2), fill = NA)+
  labs(x = "Latitude", y = "Longitude", fill="SPEI",
       title = "Região do Semiárido Brasileiro - SPEI - escala 6 meses - Série temporal 1990-2019")


#plot(spei6_dry_month)


# Exporta raster como tif 
writeRaster(spei6_dry_month,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPEI_6_seca_Mes_1990-2019.tif')


###############################################################################################################################
#
#
# ANÁLISE DE SECA COM SPEI PARA A REGIÃO NES POR SECA RECORRENTE
#
#
###############################################################################################################################

# recorte temporal para os anos de 2010-2019

spei6_2010 = subset(spei6, 241:252) #recorte para 2010
spei6_2011 = subset(spei6, 253:264) #recorte para 2011
spei6_2012 = subset(spei6, 265:276) #recorte para 2012
spei6_2013 = subset(spei6, 277:288) #recorte para 2013
spei6_2014 = subset(spei6, 289:300) #recorte para 2014
spei6_2015 = subset(spei6, 301:312) #recorte para 2015
spei6_2016 = subset(spei6, 313:324) #recorte para 2016
spei6_2017 = subset(spei6, 325:336) #recorte para 2017
spei6_2018 = subset(spei6, 337:348) #recorte para 2018
spei6_2019 = subset(spei6, 349:360) #recorte para 2019

# Definindo eventos de seca na série (2010-2019) quando o limiar de SPEI < -0.5
# reclassificando para valores binário (valores de seca = 1)

reclass <- cbind(from = c(-Inf,-0.5), to = c(-0.5, Inf), becomes = c(1, 0))

spei6_2010 <- reclassify(spei6_2010, reclass)
spei6_2011 <- reclassify(spei6_2011, reclass)
spei6_2012 <- reclassify(spei6_2012, reclass)
spei6_2013 <- reclassify(spei6_2013, reclass)
spei6_2014 <- reclassify(spei6_2014, reclass)
spei6_2015 <- reclassify(spei6_2015, reclass)
spei6_2016 <- reclassify(spei6_2016, reclass)
spei6_2017 <- reclassify(spei6_2017, reclass)
spei6_2018 <- reclassify(spei6_2018, reclass)
spei6_2019 <- reclassify(spei6_2019, reclass)

# contagem de recorrência da seca entre 2010-2019

spei6_dry_2010_2019 = (spei6_2010 + spei6_2011 + spei6_2012 + spei6_2013 + spei6_2014 +
                          spei6_2015 + spei6_2016 + spei6_2017 + spei6_2018 + spei6_2019)

# atualização dos nomes de cada mês
names(spei6_dry_2010_2019) <- c("Janeiro", "Fevereiro", "Marco", "Abril", "Maio", "Junho", "Julho", "Agosto", "Setembro", "Outubro", "Novembro", "Dezembro")


# Mapas de recorrência de seca (quantidade de vezes que o pixel apresentou SPEI < -0.5 entre 2010-2019)
# visualização dos mapas
cuts=c(0,1,2,3,4,5,6,7,8) #define os intervalos da legenda
pal <- colorRampPalette(c("grey","yellow","orange","red")) #define as cores da legenda

plot(spei6_dry_2010_2019, breaks=cuts, col = pal(8)) #gera o mapa com as cores e legendas pré-definidas


# Mapas mais elegantes

gplot(spei6_dry_2010_2019) + 
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  scale_fill_gradient2(low = "white", mid = "yellow", high = "red", midpoint= 4, na.value = NA) +
  coord_equal()+
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))+
  geom_polygon(data = BR_plot, aes(x = long, y = lat, group = group), colour = alpha("black", 1/2), fill = NA)+
  labs(x = "Latitude", y = "Longitude", fill="Anos",
       title = "Região do Semiárido Brasileiro - recorrência de seca - SPEI - Série temporal 2010-2019")



# Exporta raster como tif 
writeRaster(spei6_dry_2010_2019,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPEI_6_seca_2010-2019.tif')


###############################################################################################################################
#
#
# ANÁLISE DE SECA COM SPEI PARA A REGIÃO NES NA SÉRIE HISTÓRICA POR ESTAÇÃO DO ANO
#
#
###############################################################################################################################


#---------- Mapa - Verão


NES_summer_2015 <- subset(spei6, 310:315) # Out-Mar para 2015/2016
NES_summer_2016 <- subset(spei6, 322:327) # Out-Mar para 2016/2017
NES_summer_2017 <- subset(spei6, 334:339) # Out-Mar para 2017/2018
NES_summer_2018 <- subset(spei6, 346:351) # Out-Mar para 2018/2019


NES_summer_2015_mean <- calc(NES_summer_2015, fun = mean)
NES_summer_2016_mean <- calc(NES_summer_2016, fun = mean)
NES_summer_2017_mean <- calc(NES_summer_2017, fun = mean)
NES_summer_2018_mean <- calc(NES_summer_2018, fun = mean)

# Somente eventos de seca
NES_summer_2015_mean[NES_summer_2015_mean[] > -0.5 ] = NA
NES_summer_2016_mean[NES_summer_2016_mean[] > -0.5 ] = NA
NES_summer_2017_mean[NES_summer_2017_mean[] > -0.5 ] = NA
NES_summer_2018_mean[NES_summer_2018_mean[] > -0.5 ] = NA

NES_summer_2015_2019_map = stack(NES_summer_2015_mean, NES_summer_2016_mean, NES_summer_2017_mean, NES_summer_2018_mean)  

names(NES_summer_2015_2019_map)<- c("2015-2016 Summer (Oct-Mar)", "2016-2017 Summer(Oct-Mar)", "2017-2018 Summer(Oct-Mar)","2018-2019 Summer(Oct-Mar)")

# Visualização do mapa

cuts=c(0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow")) #define as cores da legenda

plot(NES_summer_2015_2019_map, breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas



# Visualização dos mapas individualmente sobre o Semiárido brasileiro
plot(NES_summer_2015_2019_map$X2015.2016.Summer..Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_summer_2015_2019_map$X2016.2017.Summer.Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_summer_2015_2019_map$X2017.2018.Summer.Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_summer_2015_2019_map$X2018.2019.Summer.Oct.Mar., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)


# Exporta raster como tif 
writeRaster(NES_summer_2015_2019_map,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPEI_dry_summer_2015-2019.tif')


######################################################################################################################################################


#---------- Mapa - Inverno


NES_winter_2015 <- subset(spei6, 304:309) # abr-set para 2015
NES_winter_2016 <- subset(spei6, 316:321) # abr-set para 2016
NES_winter_2017 <- subset(spei6, 328:333) # abr-set para 2017
NES_winter_2018 <- subset(spei6, 340:345) # abr-set para 2018
NES_winter_2019 <- subset(spei6, 352:357) # abr-set para 2019

NES_winter_2015_mean <- calc(NES_winter_2015, fun = mean)
NES_winter_2016_mean <- calc(NES_winter_2016, fun = mean)
NES_winter_2017_mean <- calc(NES_winter_2017, fun = mean)
NES_winter_2018_mean <- calc(NES_winter_2018, fun = mean)
NES_winter_2019_mean <- calc(NES_winter_2019, fun = mean)

# Somente eventos de seca
NES_winter_2015_mean[NES_winter_2015_mean[] > -0.5 ] = NA
NES_winter_2016_mean[NES_winter_2016_mean[] > -0.5 ] = NA
NES_winter_2017_mean[NES_winter_2017_mean[] > -0.5 ] = NA
NES_winter_2018_mean[NES_winter_2018_mean[] > -0.5 ] = NA
NES_winter_2019_mean[NES_winter_2019_mean[] > -0.5 ] = NA

NES_winter_2015_2019_map = stack(NES_winter_2015_mean, NES_winter_2016_mean, NES_winter_2017_mean, NES_winter_2018_mean, NES_winter_2019_mean)  

names(NES_winter_2015_2019_map)<- c("2015 winter (abr-set)", "2016 winter (abr-set)", "2017 winter (abr-set)","2018 winter (abr-set)", "2019 winter (abr-set)")

# Visualização do mapa

cuts=c(0,-0.2,-0.4,-0.6,-0.8,-1.0,-1.2,-1.4,-1.6,-1.8,-2.0) #define os intervalos da legenda
pal <- colorRampPalette(c("red","orange","yellow")) #define as cores da legenda

plot(NES_winter_2015_2019_map, breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas



# Visualização dos mapas individualmente sobre o Semiárido brasileiro
plot(NES_winter_2015_2019_map$X2015.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2016.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2017.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2018.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)

plot(NES_winter_2015_2019_map$X2019.winter..abr.set., breaks=cuts, col = pal(11)) #gera o mapa com as cores e legendas pré-definidas
plot(SAB, bg="transparent", add=TRUE)


# Exporta raster como tif 
writeRaster(NES_winter_2015_2019_map,'SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/NES_SPEI_dry_winter_2015-2019.tif')




#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#
#
#
# Autor:	Thales Vaz Penha | thales.penha@usp.br														
#				 					
#				
#																			
# Software: R version 4.0.4 (2021-02-15)
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)
#
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################


