# -------------------------------------------------------------------------	#
#
#
#					SCRIPT PARA C�LCULO DOS �NDICES DE SECA 
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
# Intitui��o: University of East Aglia Cimatic Research Unit (CRU)
# Fonte: https://catalogue.ceda.ac.uk/uuid/c26a65020a5e4b80b20018f148556681 | https://crudata.uea.ac.uk/cru/data/hrg/ | https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/
# Resolu��o espacial: 0.5�
# Janela temporal do dado: 1901-2020
# Formato: NetCDF
# Refer�ncias: Harris, I., Osborn, T.J., Jones, P. et al. Version 4 of the CRU TS monthly high-resolution gridded multivariate climate dataset. Sci Data 7, 109 (2020). https://doi.org/10.1038/s41597-020-0453-3
#
#
# -------------------------------------------------------------------------	#

############################################################################################################


rm(list=ls()) # remove da mem�ria global os dados anteriormente trabalhados


memory.limit (9999999999) # expande o limite de mem�ria de armazenamento e processamento de dados no R


setwd("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/") # Define o diret�rio de �rea de trabalho


############################################################################################################


# ----- Bibliotecas exigidas ------ # caso n�o tenha as bibliotecas instaladas, utilize o comando "install.packages" no lugar de "library" e depois utilize o comando "library"

library('raster') ; library('rgdal') ; library('ncdf4') ; library('utils')
library('sp') ; library('RNetCDF') ; library('SPEI') ; library('rasterVis')


############################################################################################################
#
#
# -------------------------------------------------------------------------	#
#
#
# ----- Defini��o da �rea de estudo
# 
# Limites da regi�o Northeastern South America (NES) do IPCC-WGI
# Fonte: https://github.com/SantanderMetGroup/ATLAS/tree/v1.6/reference-regions
# Refer�ncia: https://essd.copernicus.org/articles/6/2959/2020/essd-6-2959-2020.html
# 
# Limites do Semi�rido Brasileiro
# Fonte: http://antigo.sudene.gov.br/delimitacao-do-semiarido
#
#
# -------------------------------------------------------------------------	#


# Op��o 1: Abrir shapefile com os limites da �rea de estudo

NES = readOGR("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/Shapefile/IPCC-WGI-NES_region_v4.shp")

# Op��o 2: Definir o quadrangular que engloba a �rea de estudo

coords = matrix(c(-34.0, -20.0,
                  -50.0, -20.0,
                  -50.0, 0.0,
                  -34.0, 0.0), 
                  ncol = 2, byrow = TRUE)

NES = Polygon(coords)
NES = SpatialPolygons(list(Polygons(list(NES), ID = "a")), proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
plot(NES, axes = TRUE)

############################################################################################################


# conhecendo aspectos b�sicos do dado de entrada de precipita��o e suas caracter�sticas

# Abrir o arquivo de entrada NetCDF de precipita��o do CRU

pr <- nc_open("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/NetCDF/cru_ts4.05.1901.2020.pre.dat.nc")

# Visualizar as caracter�sticas do arquivo de entrada NetCDF de precipita��o do CRU

print(pr)

# obter informa��es da dimens�o temporal do dado de entrada de precipita��o

time <- ncvar_get(pr, "time") # 1440 intervalos (meses) entre 1900 e 2020
time <- as.vector(time)

# verificando a data de in�cio da s�rie

tunits <- ncatt_get(pr,"time","units")
tunits


############################################################################################################
#
#           Recorte espa�o-temporal dos dados de entrada para os limites da �rea de estudo
#
#
#               PRECIPITA��O - dado de entrada necess�rio para o c�lculo do SPI e SPEI
#
#                                                e 
#
#           EVAPOTRANSPIRA��O POTENCIAL - dado de entrada necess�rio para o c�lculo do SPEI
#
#
############################################################################################################


# PRECIPITA��O  - CRU TS 4.05 dataset (1901-2020) - https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/cruts.210305643.v4.05/pre/


############################################################################################################


# abre e transforma o arquivo NetCDF de precipita��o em Raster para manipula��o

precipitation = raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/NetCDF/cru_ts4.05.1901.2020.pre.dat.nc")

# Verifica��o da estrutura do dado

str(precipitation)


# Identifica��o das datas de in�cio e fim da s�rie do dado

precipitation_date <- getZ(precipitation)            # extra��o da informa��o de data
grep("1990-01-16", precipitation_date)    # valor [1069] corresponde a data "01-01-1990"
grep("2019-6-16", precipitation_date)    # valor [1428] corresponde ao final da s�rie de 30 anos "01-6-2019" 

# Recorte temporal do dado de entrada. Em seguida, recupera��o da informa��o de data.
# precipitation_30y � a vari�vel de precipita��o com a s�rie temporal de 30 anos

precipitation_30y <- subset(precipitation, 1069:1428)          # recorte temporal para 1990-2019
precipitation_30y@z$Date <- precipitation@z$Date[1069:1428]    # acrescentando o campo $Date de volta ao campo @z 


# Start the clock | contagem do tempo de processamento
ptm <- proc.time()

# lopping com recorte espacial do dado de entrada
for(i in 1:360){
  message(paste("reading layer", i))
  # Preparo do raster layer para o ano i
  precipitation_30y_NES <- precipitation_30y  
  #Recorte espacial com os limites da �rea de estudo
  precipitation_30y_NES <- crop(precipitation_30y_NES, NES)
  #Recorte espacial com extra��o das informa��es da grade raster nos limites da �rea de estudo
  precipitation_30y_NES <- mask(precipitation_30y_NES, NES)
  gc()
}

# Stop the clock | tempo decorrido
proc.time() - ptm

# renomeamento da vari�vel c�lculada

prec <- precipitation_30y_NES

#######################################################################################################
#
#
# Transforma o novo dado recortado espa�ialmente e temporalmente em arquivo NETCDF para facilitar a manipula��o nos c�lculos dos �ndices de seca
#
#
########################################################################################################

# nome de sa�da do arquivo NetCDF a ser salvo

filename <- "cru_ts4.05.1990.2019.precipitation.nc"

# Informa��es de Longitude e Latitude 

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
                  vals = 1:360, # entrar com toda a s�rie desejada, neste caso, 30 anos = 360 meses #
                  unlim = TRUE,
                  longname = "Month_of_year")

# Defini��o das vari�veis da componente de precipita��o

var_prec <- ncvar_def(name = "precipitation",
                      units = "mm/month",
                      dim = list(lon, lat, time),
                      longname = "Monthly_Total_Precipitation",
                      missval = mv,
                      compression = 9)

# Adicionando as vari�veis ao arquivo netcdf

ncout <- nc_create(filename, var_prec, force_v4 = TRUE)
print(paste("The file has", ncout$nvars,"variables"))
print(paste("The file has", ncout$ndim,"dimensions"))

# adicionando algumas informa��es globais sobre o dado

ncatt_put(ncout, 0, "CRU TS 4.05 precipitation", "Modified using Thales Vaz Penha (2022) code ")
ncatt_put(ncout, 0, "Source Data obtained from","British Atmospheric Data Centre, RAL, UK")
ncatt_put(ncout, 0, "References Information", "The original data is available at http://badc.nerc.ac.uk/data/cru/")
ncatt_put(ncout, 0, "Created on", date())


# Start the clock | contagem do tempo de processamento
ptm <- proc.time()

# Colocando os valores de precipita��o ao arquivo netcdf
# � preciso fazer um loop atrav�s das camadas para adicionar os valores de precipita��o e casar com o valor de tempo corretamente

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


# EVAPOTRANSPIRA��O POTENCIAL  - CRU TS 4.05 dataset (1901-2020) - https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.05/cruts.210305643.v4.05/pet/


############################################################################################################


# abre e transforma o arquivo NetCDF de evapipita��o em Raster para manipula��o

evapotranspiration = raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_entrada/NetCDF/cru_ts4.05.1901.2020.pet.dat.nc")

# Verifica��o da estrutura do dado

str(evapotranspiration)


# Identifica��o das datas de in�cio e fim da s�rie do dado

evapotranspiration_date <- getZ(evapotranspiration)            # extra��o da informa��o de data
grep("1990-01-16", evapotranspiration_date)    # valor [1069] corresponde a data "01-01-1990"
grep("2019-6-16", evapotranspiration_date)    # valor [1428] corresponde ao final da s�rie de 30 anos "01-6-2019"  

# Recorte temporal do dado de entrada. Em seguida, recupera��o da informa��o de data.
# evapotranspiration_30y � a vari�vel de evapotranspira��o potencial com a s�rie temporal de 30 anos

evapotranspiration_30y <- subset(evapotranspiration, 1069:1428)          # recorte temporal para 1990-2019
evapotranspiration_30y@z$Date <- evapotranspiration@z$Date[1069:1428]    # acrescentando o campo $Date de volta ao campo @z 



# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# lopping com recorte espacial do dado de entrada

for(i in 1:360){
  message(paste("reading layer", i))
  # Preparo do raster layer para o ano i
  evapotranspiration_30y_NES <- evapotranspiration_30y  
  #Recorte espacial com os limites da �rea de estudo
  evapotranspiration_30y_NES <- crop(evapotranspiration_30y_NES, NES)
  #Recorte espacial com extra��o das informa��es da grade raster nos limites da �rea de estudo
  evapotranspiration_30y_NES <- mask(evapotranspiration_30y_NES, NES)
  gc()
}

# Stop the clock | tempo decorrido
proc.time() - ptm


# renomeamento da vari�vel c�lculada

evap <- evapotranspiration_30y_NES

#######################################################################################################
#
#
# Transforma o novo dado recortado espa�ialmente e temporalmente em arquivo NETCDF para facilitar a manipula��o nos c�lculos dos �ndices de seca
#
#
########################################################################################################

# nome de sa�da do arquivo NetCDF a ser salvo

filename <- "cru_ts4.05.1990.2019.evapotranspiration.nc"

# Informa��es de Longitude e Latitude 

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
                  vals = 1:360, # entrar com toda a s�rie desejada, neste caso, 30 anos = 360 meses #
                  unlim = TRUE,
                  longname = "Month_of_year")

# Defini��o das vari�veis da componente de evapotrasnpira��o potencial

var_evap <- ncvar_def(name = "evapotranspiration",
                      units = "mm/month",
                      dim = list(lon, lat, time),
                      longname = "Monthly_Total_evapotranspiration",
                      missval = mv,
                      compression = 9)

# Adicionando as vari�veis ao arquivo netcdf
ncout <- nc_create(filename, var_evap, force_v4 = TRUE)
print(paste("The file has", ncout$nvars,"variables"))
print(paste("The file has", ncout$ndim,"dimensions"))

# adicionando algumas informa��es globais sobre o dado

ncatt_put(ncout, 0, "CRU TS 4.05 evapotranspiration", "Modified using Thales Vaz Penha (2022) code ")
ncatt_put(ncout, 0, "Source Data obtained from","British Atmospheric Data Centre, RAL, UK")
ncatt_put(ncout, 0, "References Information", "The original data is available at http://badc.nerc.ac.uk/data/cru/")
ncatt_put(ncout, 0, "Created on", date())


# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# Colocando os valores de evapotranspira��o potencial ao arquivo netcdf
# � preciso fazer um loop atrav�s das camadas para adicionar os valores de evapotranspira��o potencial e casar com o valor de tempo corretamente

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
# Manipula��o do dado de precipita��o para calcular o �ndice de seca SPI
#
#
############################################################################################################


# Abrindo o arquivo raster (NETCDF) com as informa��es de precipita��o mensal

precip <- nc_open("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.precipitation.nc")


# Extra��o das vari�veis de precipita��o para transforma��o de raster para matriz array


lon <- length(ncvar_get(precip, "longitude")) # verifica o n�mero de longitudes
lat <- length(ncvar_get(precip, "latitude")) # verifica o n�mero de latitudes
lon1 <- ncvar_get(precip, "longitude") # extrai os valores de longitudes
lat1 <- ncvar_get(precip, "latitude") # extrai os valores de latitudes
ppt  <- ncvar_get(precip, "precipitation") # extrai os valores de precipita��o
#ppt  <- ppt[ , ,1069:1428] # 1990-2020 # Cliamtologia dos �ltimos 30 anos] # long:40 lat:32 prec:360
Anual = 360/12 #anos na s�rie temporal

############################################################################################################

#Redimensionando o dado de entrada para uma matriz array

precip <- sapply(1:dim(ppt)[3], function(x)t(ppt[,,x]))

############################################################################################################
#
#
# C�lculo do SPI com escala de 6 meses (SPI-6)
#
#
############################################################################################################

#Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# constru��o do array de tr�s dimens�es

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

#Salvando o �ndice c�lculado em formato NETCDF

for(i in 1:360) { 
  nam <- paste("SPI", i, sep = "")
  assign(nam,raster((spi_c[ , ,i]), xmn=min(lon1), xmx=max(lon1), ymn=min(lat1), ymx=max(lat1), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")) )
}

# agrega todos os raster SPI calculados em um �nico raster
CRU_spi <- stack(mget(paste0("SPI", 1:360)))

# compatibiliza a resolus�o espacial do dado SPI com o de entrada
CRU_spi <- resample(CRU_spi, precip)


#remove os SPEI em demasia
rm(list=ls(pattern="SPI"))

# Exporta o resultado do �ndice em formato NETCDF
outfile <- "spi6_cru_1990_2019.nc"
crs(CRU_spi) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(CRU_spi, outfile, overwrite=TRUE, format="CDF", varname="SPI", varunit="units",longname="SPI CRU", xname="lon", yname="lat")


############################################################################################################
#
#            SPEI
#
############################################################################################################


# Manipula��o do dado de precipita��o e evapotranspira��o potencial para calcular o �ndice de seca SPEI


# Abrindo o arquivo raster (NETCDF) com as informa��es de precipita��o mensal

precip <- raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.precipitation.nc")

# Abrindo o arquivo raster (NETCDF) com as informa��es de evapotranspira��o potencial mensal

evap.transp <- raster::brick("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.evapotranspiration.nc")



# Transforma��es de unidades 

precip # unidade: mm/m�s

evap.transp # unidade: mm/dia


#compatibilizando as unidades

evap.transp = (evap.transp * 30) # mm/dia para mm/m�s



# C�lculo do balan�o h�drico clim�tico com base no saldo de precipita��o menos evapotranspira��o


BHC = (precip - evap.transp)



# nome de sa�da do arquivo NetCDF a ser salvo

filename <- "cru_ts4.05.1990.2019.climatic_water_balance.nc"

# Informa��es de Longitude and Latitude

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
                  vals = 1:360, # entrar com toda a s�rie desejada, neste caso, 30 anos = 360 meses #
                  unlim = TRUE,
                  longname = "Month_of_year")

# Defini��o das vari�veis da componente de balan�o h�drico clim�tico

var_bhc <- ncvar_def(name = "cwb",
                      units = "mm/month",
                      dim = list(lon, lat, time),
                      longname = "Monthly_Total_climatic_water_balance",
                      missval = mv,
                      compression = 9)


# Adicionando as vari�veis ao arquivo netcdf

ncout <- nc_create(filename, var_bhc, force_v4 = TRUE)
print(paste("The file has", ncout$nvars,"variables"))
print(paste("The file has", ncout$ndim,"dimensions"))

# adicionando algumas informa��es globais sobre o dado

ncatt_put(ncout, 0, "CRU TS 4.05 climatic water balance", "Modified using Thales Vaz Penha (2022) code ")
ncatt_put(ncout, 0, "Source Data obtained from","British Atmospheric Data Centre, RAL, UK")
ncatt_put(ncout, 0, "References Information", "The original data is available at http://badc.nerc.ac.uk/data/cru/")
ncatt_put(ncout, 0, "Created on", date())


# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# Colocando os valores de balan�o h�drico clim�tico ao arquivo netcdf
# � preciso fazer um loop atrav�s das camadas para adicionar os valores de balan�o hidrico clim�tico e casar com o valor de tempo corretamente

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


# Abrindo o arquivo raster (NETCDF) com as informa��es de balan�o h�drico clim�tico mensal

clim_water_balance <- nc_open("SeuDiretorio(modifique aqui)/Dados_Penha_&_Cunha_2022/Dados_saida/cru_ts4.05.1990.2019.climatic_water_balance.nc")


# Extra��o das vari�veis de precipita��o para transforma��o de raster para matriz


lon <- length(ncvar_get(clim_water_balance, "longitude")) # verifica o n�mero de longitudes
lat <- length(ncvar_get(clim_water_balance, "latitude")) # verifica o n�mero de latitudes
lon1 <- ncvar_get(clim_water_balance, "longitude") # extrai os valores de longitudes
lat1 <- ncvar_get(clim_water_balance, "latitude") # extrai os valores de latitudes
ppt  <- ncvar_get(clim_water_balance, "cwb") # extrai os valores de precipita��o
#ppt  <- ppt[ , ,1069:1428] # 1990-2020 # Cliamtologia dos �ltimos 30 anos] # long:40 lat:32 prec:360
Anual = 360/12 #anos na s�rie temporal

############################################################################################################

#Redimensionando o dado de entrada para uma matriz

clim_water_balance <- sapply(1:dim(ppt)[3], function(x)t(ppt[,,x]))

############################################################################################################


# C�lculo do SPEI com escala de 6 meses (SPEI-6)


############################################################################################################
# Start the clock! | contagem do tempo de processamento
ptm <- proc.time()

# constru��o do array de tr�s dimens�es

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

#Salvando o �ndice c�lculado em formato NETCDF

for(i in 1:360) { 
  nam <- paste("SPEI", i, sep = "")
  assign(nam,raster((spei_c[ , ,i]), xmn=min(lon1), xmx=max(lon1), ymn=min(lat1), ymx=max(lat1), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")) )
}

# agrega todos os raster SPEI calculados em um �nico raster
CRU_spei <- stack(mget(paste0("SPEI", 1:360)))

# compatibiliza a resolus�o espacial do dado SPEI com o de entrada
CRU_spei <- resample(CRU_spei, clim_water_balance)

#remove os SPEI em demasia
rm(list=ls(pattern="SPEI"))

# Exporta o resultado do �ndice em formato NETCDF
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
