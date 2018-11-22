rm(list = ls())
library("raster")
rcl<- function(x,y){
  reclassify(x, y, include.lowest=FALSE, right=NA)}

corte <- function(x){substr(x,pos,pos+6)}

# lee archivos .tif de un directorio
f.1 <- function(x,y){
  x<- list.files(y,corte(x),full.names = T)
  return(x[grepl(x = x,pattern = ".tif$")])}


#########
#### funciones de performance
m <- function(x){
  return(timestamp(quiet = T))
}

f.t2<- function(x){
  hh<- substr(x,21,22);mm<- substr(x,24,25);ss<- substr(x,27,28)
  hh[2]<- substr(x[2],21,22);mm[2]<- substr(x[2],24,25);ss[2]<- substr(x[2],27,28)
  hh <- as.numeric(hh);mm <- as.numeric(mm);ss <- as.numeric(ss);
  if((ss[2]-ss[1])>=0){ss[3] <- (ss[2]-ss[1])}else{ss[3] <- 60+(ss[2]-ss[1]);mm[1] <- mm[1]+1 }
  if((mm[2]-mm[1])>=0){mm[3] <-  mm[2]-mm[1]}else{mm[3] <-  60+(mm[2]-mm[1]);hh[1] <- hh[1]+1}
  if((hh[2]-hh[1])>=0){hh[3] <-  hh[2]-hh[1]}else{hh[3] <-  24+(hh[2]-hh[1])}
  return(paste(sprintf("%02d",hh[3]),sprintf("%02d",mm[3]),sprintf("%02d",ss[3]),sep=":"))
}
pr <- function(x){
  unique(getValues(x))
}
###################

#   ##############################################################
# 0-40: Soil              # 237: inland water       
# 40-100: snow cover      # 239: ocean
# 200: missing data       # 250: cloud
# 201: no decision        # 254: detector saturated
# 211: night              # 255: fill
#   ##############################################################
dir.mod <-"~/TESIS/test/mod/" 
dir.myd <- "~/TESIS/test/myd/" 
dir.mod.c <-"~/TESIS/test/c_mod/" 
dir.myd.c <- "~/TESIS/test/c_myd/" 
dir.mod.myd.c <- "~/TESIS/test/c_mod_myd/" 
dir.mod.myd <- "~/TESIS/test/mod_myd/" 
# lmod <- list.files("~/TESIS/test/mod/")
# lmyd <- list.files("~/TESIS/test/myd/")
pos <- 10
lmod <- list.files(path=dir.mod,pattern = ".tif$")
lmyd <- list.files(path=dir.myd,pattern = ".tif$")
# lmod <- f.1(lmod,dir.mod)
# lmyd <- f.1(lmyd,dir.myd)
nodata <- c(200,201,211,237,239,254,250,255)

suelo <- c(seq(0,39))
nieve <- c(seq(40,100))
acc <- vector()
tabla <- data.frame()

### matríz para armar la imágenes de nubes
capa.nubes <- matrix(ncol=2,c(250,seq(0,100),200,201,211,237,239,254,255,1,rep(NA,108)))

snow.bare.clouds <- as.matrix(data.frame(col1=c(suelo,nieve,nodata),
                                         col2=c(rep(0,length(suelo)),
                                                rep(1,length(nieve)),
                                                rep(NA,length(nodata)))))

bool.clouds <- as.matrix(data.frame(col1=c(NA,0,1),col2=c(1,0,0)))

na20 <- as.matrix(data.frame(col1=c(NA),col2=c(0)))
cero2na <- as.matrix(data.frame(col1=c(0),col2=c(NA)))

###############################################################################################
{d1 <- m()
for(i in 1:5000){
  # Leo las imágenes
  # levanto los rasters
  mod <- raster("~/TESIS/test/MOD10A1.A2016089.h12v12.006.2016104051228.NDSI_Snow_Cover.tif")
  
  # mod <- raster(paste(dir.mod,lmod,sep="/"))
  c.mod<- rcl(mod, capa.nubes)
  
  mod<- rcl(mod, snow.bare.clouds)# NA = nubes, 0 = suelo, 1 = nieve
  
  
  # grabo mod, grabo c.mod, renombro c.mod
  ### solo para el testing
  lmod <- list.files(path=dir.mod,pattern = ".tif$")
  writeRaster(mod,paste(dir.mod,lmod,sep="/"),format="GTiff", overwrite=T,datatype='INT1U')
  writeRaster(c.mod,paste(dir.mod.c,"/MOD10A1.A",corte(lmod),".clouds.cover.tif"
                          ,sep=""),format="GTiff", overwrite=T,datatype="INT1U")
  file.rename(paste(dir.mod,lmod,sep="/"),
              paste(dir.mod,"/",substr(lmod,1,nchar(lmod)-8),"_",
                    sprintf("%03d",i),".tif",sep=""))
  
  #############################################################################################
  
  myd <- raster("~/TESIS/test/MYD10A1.A2016089.h12v12.006.2016091123302.NDSI_Snow_Cover.tif")
  
  # myd <- raster(paste(dir.myd,lmyd[j],sep="/"))
  c.myd<- rcl(myd, capa.nubes)# 0 = nubes, 1 = suelo, 2 = nieve
  myd<- rcl(myd, snow.bare.clouds)# 0 = nubes, 1 = suelo, 2 = nieve
  lmyd <- list.files(path=dir.myd,pattern = ".tif$")
  writeRaster(myd,paste(dir.myd,lmyd,sep="/"),format="GTiff", overwrite=T,datatype='INT1U')
  writeRaster(c.myd,paste(dir.myd.c,"/MYD10A1.A",corte(lmyd),".clouds.cover.tif"
                          ,sep=""),format="GTiff", overwrite=T,datatype="INT1U")
  file.rename(paste(dir.myd,lmyd,sep="/"),
              paste(dir.myd,"/",substr(lmyd,1,nchar(lmyd)-8),
                    "_",sprintf("%03d",i),".tif",sep=""))
  
  #############################################################################################
  # c.mod[which(is.na(getValues(c.mod)))] <- 0
  # c.myd[which(is.na(getValues(c.myd)))] <- 0
  writeRaster((rcl(c.mod,na20)+rcl(c.myd,na20)- rcl(c.mod,na20)*rcl(c.myd,na20)),
              paste(dir.mod.myd.c,"/MOD.MYD.A",substr(lmyd,pos,pos+6),
                    ".clouds.max.tif",sep=""),format="GTiff", overwrite=T,datatype="INT1U")
  ### este está ok así
  writeRaster((c.mod*c.myd),paste(dir.mod.myd.c,"/MOD.MYD.A",substr(lmyd,pos,pos+6),
                                  ".clouds.min.tif",sep=""),format="GTiff", overwrite=T,
              datatype="INT1U")
  ###########################################################################################
  
  
  mask.mod <-  reclassify(mod, bool.clouds, include.lowest=FALSE, right=NA) 
  
  mod.myd <- mask.mod*rcl(myd+1,na20)
  mod.myd <- rcl(mod+1,na20) + mod.myd
  mod.myd<- (rcl(mod.myd,cero2na)-1)
  
  writeRaster(mod.myd,paste(dir.mod.myd,"/MOD.MYD.A",substr(lmyd,pos,pos+6),
                            ".snow.cover.area.tif",sep=""),format="GTiff", overwrite=T,
              datatype='INT1U')
  for(h in 1:7){
    
    tab<- getValues(mod.myd)
    largo<- mod.myd@ncols*mod.myd@nrows
    clouds.cover<- length(na.omit(getValues(mod.myd)))
    snow.cover<- length(which(getValues(mod.myd)==1))
    tabla[i,1] <- corte(lmod)
    tabla[i,2] <- (clouds.cover*100)/largo
    tabla[i,3] <- snow.cover*100/(tabla[i,2]*largo/100)
    tabla[i,4]
  }          
}

d1[2] <- m()       
f.t2(d1)
write.csv(d1,paste(dir.mod,"t1.csv"))
}

### 
# Este test concluye en unos resultados para las imágenes de 
# mod/myd NA 0 1 
# cmod/cmyd NA 1
# modmyd NA 0 1 
# 
