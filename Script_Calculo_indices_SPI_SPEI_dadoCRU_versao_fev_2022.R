# -------------------------------------------------------------------------	#
#
#
#					SCRIPT PARA CÁLCULO DOS ÍNDICES DE SECA 
#
#           STANDARD PRECIPITATION INDEX (SPI)
#                 
#   STANDARD PRECIPITATION-EVAPORTRANSPIRATION INDEX (SPEI)  
#		
#
#
# Criado:		28-09-2021						
# Modificado:	21-02-2022		  
#		
#
# Autor:	Thales Vaz Penha | thales.penha@usp.br														
#				 					
#				
#																			
# Software: R version 4.0.4 (2021-02-15)
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# -------------------------------------------------------------------------	#
#
# Dado de entrada: 
#
# CRU TS 4.05 dataset
#
# Intituição: University of East Aglia Cimatic Research Unit (CRU)
# Fonte: https://catalogue.ceda.ac.uk/uuid/c26a65020a5e4b80b20018f148556681 | https://crudata.uea.ac.uk/cru/data/hrg/ | https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/
# Resolução espacial: 0.5°
# Janela temporal do dado: 1901-2020
# Formato: NetCDF
# Referências: Harris, I., Osborn, T.J., Jones, P. et al. Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. Sci Data 7, 109 (2020). https://doi.org/10.1038/s41597-020-0453-3
#
#
# -------------------------------------------------------------------------	#

############################################################################################################


rm(list=ls()) # remove da memória global os dados anteriormente trabalhados


memory.limit (9999999999) # expande o limite de memória de armazenamento e processamento de dados no R


setwd("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/") # Define o diretório de área de trabalho


############################################################################################################


# ----- Bibliotecas exigidas ------ # caso não tenha as bibliotecas instaladas, utilize o comando "install.packages" no lugar de "library" e depois utilize o comando "library"

library('raster') ; library('rgdal') ; library('ncdf4') ; library('utils')
library('sp') ; library('RNetCDF') ; library('SPEI') ; library('rasterVis')


############################################################################################################
#
#
# -------------------------------------------------------------------------	#
#
#
# ----- Definição da área de estudo
# 
# Limites da região Northeastern South America (NES) do IPCC-WGI
# Fonte: https://github.com/SantanderMetGroup/ATLAS/tree/v1.6/reference-regions
# Referência: https://essd.copernicus.org/articles/6/2959/2020/essd-6-2959-2020.html
# 
# Limites do Semiárido Brasileiro
# Fonte: http://antigo.sudene.gov.br/delimitacao-do-semiarido
#
#
# -------------------------------------------------------------------------	#


# Opção 1: Abrir shapefile com os limites da área de estudo

NES = readOGR("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/Shapefile/IPCC-WGI-NES_region_v4.shp")

# Opção 2: Definir o quadrangular que engloba a área de estudo

coords = matrix(c(-34.0, -20.0,
                  -50.0, -20.0,
                  -50.0, 0.0,
                  -34.0, 0.0), 
                  ncol = 2, byrow = TRUE)

NES = Polygon(coords)
NES = SpatialPolygons(list(Polygons(list(NES), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(NES, axes = TRUE)

############################################################################################################


# conhecendo aspectos básicos do dado de entrada de precipitação e suas características

# Abrir o arquivo de entrada NetCDF de precipitação do CRU

pr <- nc_open("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/NetCDF/cru_ts4.05.1901.2020.pre.dat.nc")

# Visualizar as características do arquivo de entrada NetCDF de precipitação do CRU

print(pr)

# obter informações da dimensão temporal do dado de entrada de precipitação

time <- ncvar_get(pr, "time") # 1440 intervalos (meses) entre 1900 e 2020
time <- as.vector(time)

# verificando a data de início da série

tunits <- ncatt_get(pr,"time","units")
tunits


############################################################################################################
#
#           Recorte espaço-temporal dos dados de entrada para os limites da área de estudo
#
#
#               PRECIPITAÇÃO - dado de entrada necessário para o cálculo do SPI e SPEI
#
#                                                e 
#
#           EVAPOTRANSPIRAÇÃO POTENCIAL - dado de entrada necessário para o cálculo do SPEI
#
#
############################################################################################################


# PRECIPITAçÃO  - CRU TS 4.05 dataset (1901-2020) - https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/cruts.210305643.v4.05/pre/


############################################################################################################


# abre e transforma o arquivo NetCDF de precipitação em Raster para manipulação

precipitation = raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/NetCDF/cru_ts4.05.1901.2020.pre.dat.nc")

# Verificação da estrutura do dado

str(precipitation)


# Identificação das datas de início e fim da série do dado

precipitation_date <- getZ(precipitation)            # extração da informação de data
grep("1990-01-16", precipitation_date)    # valor [1069] corresponde a data "01-01-1990"
grep("2019-6-16", precipitation_date)    # valor [1428] corresponde ao final da série de 30 anos "01-6-2019" 

# Recorte temporal do dado de entrada. Em seguida, recuperação da informação de data.
# precipitation_30y é a variável de precipitação com a série temporal de 30 anos

precipitation_30y <- subset(precipitation, 1069:1428)          # recorte temporal para 1990-2019
precipitation_30y@z$Date <- precipitation@z$Date[1069:1428]    # acrescentando o campo $Date de volta ao campo @z 


# Start the clock | contagem do tempo de processamento
ptm <- proc.time()

# lopping com recorte espacial do dado de entrada
for(i in 1:360){
  message(paste("reading layer", i))
  # Preparo do raster layer para o ano i
  precipitation_30y_NES <- precipitation_30y  
  #Recorte espacial com os limites da área de estudo
  precipitation_30y_NES <- crop(precipitation_30y_NES, NES)
  #Recorte espacial com extração das informações da grade raster nos limites da área de estudo
  precipitation_30y_NES <- mask(precipitation_30y_NES, NES)
  gc()
}

# Stop the clock | tempo decorrido
proc.time() - ptm

# renomeamento da variável cálculada

prec <- precipitation_30y_NES

#######################################################################################################
#
#
# Transforma o novo dado recortado espaçialmente e temporalmente em arquivo NETCDF para facilitar a manipulação nos cálculos dos índices de seca
#
#
########################################################################################################

# nome de saída do arquivo NetCDF a ser salvo

filename <- "cru_ts4.05.1990.2019.precipitation.nc"

# Informações de Longitude e Latitude 

xvals <- unique(values(init(prec, "x"))) # extrai os valores de longitude
yvals <- unique(values(init(prec, "y"))) # extrai os valores de latitude
nx <- length(xvals) # verifica o tamanho dos valores de longitude
ny <- length(yvals) # verifica o tamanho dos valores de latitude
lon <- ncdim_def("longitude", "degrees_east", xvals) # define a componente longitude do netcdf
lat <- ncdim_def("latitude", "degrees_south", yvals) # define a componente latitude do netcdf

# Missing value a ser usado

mv <- -9999

# Componente Tempo

time <- ncdim_def(name = "Time", 
                  units = "months", 
                  vals = 1:360, # entrar com toda a série desejada, neste caso, 30 anos = 360 meses #
                  unlim = TRUE,
                  longname = "Month_of_year")

# Definição das variáveis da componente de precipitação

var_prec <- ncvar_def(name = "precipitation",
                      units = "mm/month",
                      dim = list(lon, lat, time),
                      longname = "Monthly_Total_Precipitation",
                      missval = mv,
                      compression = 9)

# Adicionando as variáveis ao arquivo netcdf

ncout <- nc_create(filename, var_prec, force_v4 = TRUE)
print(paste("The file has", ncout$nvars,"variables"))
print(paste("The file has", ncout$ndim,"dimensions"))

# adicionando algumas informações globais sobre o dado

ncatt_put(ncout, 0, "CRU TS 4.05 precipitation", "Modified using Thales Vaz Penha (2022) code ")
ncatt_put(ncout, 0, "Source Data obtained from","British Atmospheric Data Centre, RAL, UK")
ncatt_put(ncout, 0, "References Information", "The original data is available at http://badc.nerc.ac.uk/data/cru/")
ncatt_put(ncout, 0, "Created on", date())


# Start the clock | contagem do tempo de processamento
ptm <- proc.time()

# Colocando os valores de precipitação ao arquivo netcdf
# é preciso fazer um loop através das camadas para adicionar os valores de precipitação e casar com o valor de tempo corretamente

for (i in 1:nlayers(prec)) { 
  #message("Processing layer ", i, " of ", nlayers(prec))
  ncvar_put(nc = ncout, 
            varid = "precipitation", 
            vals = values(prec[[i]]), 
            start = c(1, 1, i), 
            count = c(-1, -1, 1))
}
nc_close(ncout)

# Stop the clock | tempo decorrido
proc.time() - ptm



############################################################################################################


# EVAPOTRANSPIRAÇÂO POTENCIAL  - CRU TS 4.05 dataset (1901-2020) - https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/cruts.210305643.v4.05/pet/


############################################################################################################


# abre e transforma o arquivo NetCDF de evapipitação em Raster para manipulação

evapotranspiration = raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/NetCDF/cru_ts4.05.1901.2020.pet.dat.nc")

# Verificação da estrutura do dado

str(evapotranspiration)


# Identificação das datas de início e fim da série do dado

evapotranspiration_date <- getZ(evapotranspiration)            # extração da informação de data
grep("1990-01-16", evapotranspiration_date)    # valor [1069] corresponde a data "01-01-1990"
grep("2019-6-16", evapotranspiration_date)    # valor [1428] corresponde ao final da série de 30 anos "01-6-2019"  

# Recorte temporal do dado de entrada. Em seguida, recuperação da informação de data.
# evapotranspiration_30y é a variável de evapotranspiração potencial com a série temporal de 30 anos

evapotranspiration_30y <- subset(evapotranspiration, 1069:1428)          # recorte temporal para 1990-2019
evapotranspiration_30y@z$Date <- evapotranspiration@z$Date[1069:1428]    # acrescentando o campo $Date de volta ao campo @z 



# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# lopping com recorte espacial do dado de entrada

for(i in 1:360){
  message(paste("reading layer", i))
  # Preparo do raster layer para o ano i
  evapotranspiration_30y_NES <- evapotranspiration_30y  
  #Recorte espacial com os limites da área de estudo
  evapotranspiration_30y_NES <- crop(evapotranspiration_30y_NES, NES)
  #Recorte espacial com extração das informações da grade raster nos limites da área de estudo
  evapotranspiration_30y_NES <- mask(evapotranspiration_30y_NES, NES)
  gc()
}

# Stop the clock | tempo decorrido
proc.time() - ptm


# renomeamento da variável cálculada

evap <- evapotranspiration_30y_NES

#######################################################################################################
#
#
# Transforma o novo dado recortado espaçialmente e temporalmente em arquivo NETCDF para facilitar a manipulação nos cálculos dos índices de seca
#
#
########################################################################################################

# nome de saída do arquivo NetCDF a ser salvo

filename <- "cru_ts4.05.1990.2019.evapotranspiration.nc"

# Informações de Longitude e Latitude 

xvals <- unique(values(init(evap, "x"))) # extrai os valores de longitude
yvals <- unique(values(init(evap, "y"))) # extrai os valores de latitude
nx <- length(xvals) # verifica o tamanho dos valores de longitude
ny <- length(yvals) # verifica o tamanho dos valores de latitude
lon <- ncdim_def("longitude", "degrees_east", xvals) # define a componente longitude do netcdf
lat <- ncdim_def("latitude", "degrees_south", yvals) # define a componente latitude do netcdf


# Missing value a ser usado

mv <- -9999

# Componente Tempo

time <- ncdim_def(name = "Time", 
                  units = "months", 
                  vals = 1:360, # entrar com toda a série desejada, neste caso, 30 anos = 360 meses #
                  unlim = TRUE,
                  longname = "Month_of_year")

# Definição das variáveis da componente de evapotrasnpiração potencial

var_evap <- ncvar_def(name = "evapotranspiration",
                      units = "mm/month",
                      dim = list(lon, lat, time),
                      longname = "Monthly_Total_evapotranspiration",
                      missval = mv,
                      compression = 9)

# Adicionando as variáveis ao arquivo netcdf
ncout <- nc_create(filename, var_evap, force_v4 = TRUE)
print(paste("The file has", ncout$nvars,"variables"))
print(paste("The file has", ncout$ndim,"dimensions"))

# adicionando algumas informações globais sobre o dado

ncatt_put(ncout, 0, "CRU TS 4.05 evapotranspiration", "Modified using Thales Vaz Penha (2022) code ")
ncatt_put(ncout, 0, "Source Data obtained from","British Atmospheric Data Centre, RAL, UK")
ncatt_put(ncout, 0, "References Information", "The original data is available at http://badc.nerc.ac.uk/data/cru/")
ncatt_put(ncout, 0, "Created on", date())


# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# Colocando os valores de evapotranspiração potencial ao arquivo netcdf
# é preciso fazer um loop através das camadas para adicionar os valores de evapotranspiração potencial e casar com o valor de tempo corretamente

for (i in 1:nlayers(evap)) { 
  #message("Processing layer ", i, " of ", nlayers(evap))
  ncvar_put(nc = ncout, 
            varid = "evapotranspiration", 
            vals = values(evap[[i]]), 
            start = c(1, 1, i), 
            count = c(-1, -1, 1))
}
nc_close(ncout)

# Stop the clock | tempo decorrido
proc.time() - ptm



############################################################################################################
#
#
# Manipulação do dado de precipitação para calcular o índice de seca SPI
#
#
############################################################################################################


# Abrindo o arquivo raster (NETCDF) com as informações de precipitação mensal

precip <- nc_open("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.precipitation.nc")


# Extração das variáveis de precipitação para transformação de raster para matriz array


lon <- length(ncvar_get(precip, "longitude")) # verifica o número de longitudes
lat <- length(ncvar_get(precip, "latitude")) # verifica o número de latitudes
lon1 <- ncvar_get(precip, "longitude") # extrai os valores de longitudes
lat1 <- ncvar_get(precip, "latitude") # extrai os valores de latitudes
ppt  <- ncvar_get(precip, "precipitation") # extrai os valores de precipitação
#ppt  <- ppt[ , ,1069:1428] # 1990-2020 # Cliamtologia dos últimos 30 anos] # long:40 lat:32 prec:360
Anual = 360/12 #anos na série temporal

############################################################################################################

#Redimensionando o dado de entrada para uma matriz array

precip <- sapply(1:dim(ppt)[3], function(x)t(ppt[,,x]))

############################################################################################################
#
#
# Cálculo do SPI com escala de 6 meses (SPI-6)
#
#
############################################################################################################

#Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# construção do array de três dimensões

spi_6 <- array(list(),(lon*lat))

for (i in 1:(lon*lat)) {
  spi_6[[i]] <- SPEI::spi(precip[i,], scale=6, na.rm=TRUE)
}

# Stop the clock | tempo decorrido
proc.time() - ptm

#############################################################################################################

#Retorno ao formato array 

sapply(spi_6, '[[',2 )->matriz_ppt 
ppt_6 <- array(aperm(matriz_ppt, c(2,1),c(40,32,360)));spi_c <- array(t(ppt_6), dim=c(40,32,360))


#############################################################################################################

#Salvando o índice cálculado em formato NETCDF

for(i in 1:360) { 
  nam <- paste("SPI", i, sep = "")
  assign(nam,raster((spi_c[ , ,i]), xmn=min(lon1), xmx=max(lon1), ymn=min(lat1), ymx=max(lat1), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")) )
}

# agrega todos os raster SPI calculados em um único raster
CRU_spi <- stack(mget(paste0("SPI", 1:360)))

# compatibiliza a resolusão espacial do dado SPI com o de entrada
CRU_spi <- resample(CRU_spi, precip)


#remove os SPEI em demasia
rm(list=ls(pattern="SPI"))

# Exporta o resultado do índice em formato NETCDF
outfile <- "spi6_cru_1990_2019.nc"
crs(CRU_spi) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(CRU_spi, outfile, overwrite=TRUE, format="CDF", varname="SPI", varunit="units",longname="SPI CRU", xname="lon", yname="lat")


############################################################################################################
#
#            SPEI
#
############################################################################################################


# Manipulação do dado de precipitação e evapotranspiração potencial para calcular o índice de seca SPEI


# Abrindo o arquivo raster (NETCDF) com as informações de precipitação mensal

precip <- raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.precipitation.nc")

# Abrindo o arquivo raster (NETCDF) com as informações de evapotranspiração potencial mensal

evap.transp <- raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.evapotranspiration.nc")



# Transformações de unidades 

precip # unidade: mm/mês

evap.transp # unidade: mm/dia


#compatibilizando as unidades

evap.transp = (evap.transp * 30) # mm/dia para mm/mês



# Cálculo do balanço hídrico climático com base no saldo de precipitação menos evapotranspiração


BHC = (precip - evap.transp)



# nome de saída do arquivo NetCDF a ser salvo

filename <- "cru_ts4.05.1990.2019.climatic_water_balance.nc"

# Informações de Longitude and Latitude

xvals <- unique(values(init(BHC, "x"))) # extrai os valores de longitude
yvals <- unique(values(init(BHC, "y"))) # extrai os valores de latitude
nx <- length(xvals) # verifica o tamanho dos valores de longitude
ny <- length(yvals) # verifica o tamanho dos valores de latitude
lon <- ncdim_def("longitude", "degrees_east", xvals) # define a componente longitude do netcdf
lat <- ncdim_def("latitude", "degrees_south", yvals) # define a componente latitude do netcdf


# Missing value a ser usado

mv <- -9999

# Componente Tempo

time <- ncdim_def(name = "Time", 
                  units = "months", 
                  vals = 1:360, # entrar com toda a série desejada, neste caso, 30 anos = 360 meses #
                  unlim = TRUE,
                  longname = "Month_of_year")

# Definição das variáveis da componente de balanço hídrico climático

var_bhc <- ncvar_def(name = "cwb",
                      units = "mm/month",
                      dim = list(lon, lat, time),
                      longname = "Monthly_Total_climatic_water_balance",
                      missval = mv,
                      compression = 9)


# Adicionando as variáveis ao arquivo netcdf

ncout <- nc_create(filename, var_bhc, force_v4 = TRUE)
print(paste("The file has", ncout$nvars,"variables"))
print(paste("The file has", ncout$ndim,"dimensions"))

# adicionando algumas informações globais sobre o dado

ncatt_put(ncout, 0, "CRU TS 4.05 climatic water balance", "Modified using Thales Vaz Penha (2022) code ")
ncatt_put(ncout, 0, "Source Data obtained from","British Atmospheric Data Centre, RAL, UK")
ncatt_put(ncout, 0, "References Information", "The original data is available at http://badc.nerc.ac.uk/data/cru/")
ncatt_put(ncout, 0, "Created on", date())


# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# Colocando os valores de balanço hídrico climático ao arquivo netcdf
# é preciso fazer um loop através das camadas para adicionar os valores de balanço hidrico climático e casar com o valor de tempo corretamente

for (i in 1:nlayers(BHC)) { 
  #message("Processing layer ", i, " of ", nlayers(BHC))
  ncvar_put(nc = ncout, 
            varid = "cwb", 
            vals = values(BHC[[i]]), 
            start = c(1, 1, i), 
            count = c(-1, -1, 1))
}
nc_close(ncout)

# Stop the clock | tempo decorrido
proc.time() - ptm

#####################################################################################################################
#####################################################################################################################


# Abrindo o arquivo raster (NETCDF) com as informações de balanço hídrico climático mensal

clim_water_balance <- nc_open("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.climatic_water_balance.nc")


# Extração das variáveis de precipitação para transformação de raster para matriz


lon <- length(ncvar_get(clim_water_balance, "longitude")) # verifica o número de longitudes
lat <- length(ncvar_get(clim_water_balance, "latitude")) # verifica o número de latitudes
lon1 <- ncvar_get(clim_water_balance, "longitude") # extrai os valores de longitudes
lat1 <- ncvar_get(clim_water_balance, "latitude") # extrai os valores de latitudes
ppt  <- ncvar_get(clim_water_balance, "cwb") # extrai os valores de precipitação
#ppt  <- ppt[ , ,1069:1428] # 1990-2020 # Cliamtologia dos últimos 30 anos] # long:40 lat:32 prec:360
Anual = 360/12 #anos na série temporal

############################################################################################################

#Redimensionando o dado de entrada para uma matriz

clim_water_balance <- sapply(1:dim(ppt)[3], function(x)t(ppt[,,x]))

############################################################################################################


# Cálculo do SPEI com escala de 6 meses (SPEI-6)


############################################################################################################
# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# construção do array de três dimensões

spei_6 <- array(list(),(lon*lat))

for (i in 1:(lon*lat)) {
  spei_6[[i]] <- SPEI::spei(clim_water_balance[i,], scale=6, na.rm=TRUE)
}

# Stop the clock | tempo decorrido
proc.time() - ptm

#############################################################################################################

#Retorno ao formato array 

sapply(spei_6, '[[',2 )->matriz_ppt 
ppt_6 <- array(aperm(matriz_ppt, c(2,1),c(40,32,360)));spei_c <- array(t(ppt_6), dim=c(40,32,360))

#############################################################################################################

#Salvando o índice cálculado em formato NETCDF

for(i in 1:360) { 
  nam <- paste("SPEI", i, sep = "")
  assign(nam,raster((spei_c[ , ,i]), xmn=min(lon1), xmx=max(lon1), ymn=min(lat1), ymx=max(lat1), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")) )
}

# agrega todos os raster SPEI calculados em um único raster
CRU_spei <- stack(mget(paste0("SPEI", 1:360)))

# compatibiliza a resolusão espacial do dado SPEI com o de entrada
CRU_spei <- resample(CRU_spei, clim_water_balance)

#remove os SPEI em demasia
rm(list=ls(pattern="SPEI"))

# Exporta o resultado do índice em formato NETCDF
outfile <- "spei6_cru_1990_2019.nc"
crs(CRU_spei) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(CRU_spei, outfile, overwrite=TRUE, format="CDF", varname="SPEI", varunit="units",longname="SPEI CRU", xname="lon", yname="lat")


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
