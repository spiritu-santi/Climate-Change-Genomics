####################### auxilary functions not created by us###################
# Functions for creating site frequency spectra (SFS) from a VCF file.
# Functions for manipulating and plotting SFS.
# By Shenglin Liu, Feb 12, 2020.



## Generate a SFS (table object) from the gt object.
# It will output a SFS based on raw count without accounting for the missing values.
# pops: a character or integer vector; IDs of populations to be included in the SFS.
gt2sfs.raw<-function(gt,pops)
{
  popmap<-gt$popmap
  chrom<-gt$genotype
  
  nrow.vcf<-nrow(chrom)
  ncol.vcf<-ncol(chrom)
  
  n.pop<-length(pops)
  
  # Number of chromosomes.
  ns.chr<-sapply(pops,function(x){sum(popmap==x)})*2
  
  # SFS based on raw count.
  cnt<-matrix(0,nrow.vcf,n.pop)
  ext<-list()
  for(i in 1:n.pop)
  {
    index<-which(popmap==pops[i])
    cnt[,i]<-rowSums(chrom[,index],na.rm=T)
    ext<-c(ext,list(0:ns.chr[i]))
  }
  ext<-as.matrix(expand.grid(ext))
  cnt<-data.frame(rbind(cnt,ext))
  sfs.raw<-table(cnt)-1
  names(dimnames(sfs.raw))<-pops
  
  sfs.raw
}

## Fold a SFS.
fold.sfs<-function(sfs)
{
  sfs[]<-sfs+rev(sfs)
  dims<-dim(sfs)
  cnt.pool<-rowSums(expand.grid(lapply(dims-1,function(x)0:x)))
  index<-cnt.pool>(sum(dims-1)/2)
  sfs[index]<-0
  index<-cnt.pool==(sum(dims-1)/2)
  sfs[index]<-sfs[index]/2
  sfs
}


## Plot a SFS (barplot for 1D, image for 2D).
# cr: color range.
plot.sfs<-function(sfs,cr=c(1,20000),colScheme=function(x)rainbow(x,s=0.8),...)
{
  pops<-names(dimnames(sfs))
  dims<-dim(sfs)
  n.pop<-length(dims)
  sfs[1]<-sfs[length(sfs)]<-0
  
  # Color coding for number of SNPs.
  cr<-log10(cr)			## Color range
  nb<-101				## Number of breaks
  breaks<-10^seq(cr[1],cr[2],length.out=nb)
  
  if(n.pop==1)
  {
    barplot(sfs,names.arg=c(0,rep(NA,dims-2),dims-1),main=pops,...)
  }
  if(n.pop==2)
  {
    image(sfs,axes=F,col=colScheme(nb-1),breaks=breaks,
          xlab=pops[1],ylab=pops[2],...)
    axis(1,at=c(0,1),labels=c(0,dims[1]-1))
    axis(2,at=c(0,1),labels=c(0,dims[2]-1))
    box()
  }
  if(n.pop>2)stop("Cannot plot SFS with dimension higher than 2!")
}

## Write a 2D-SFS to a file in fastSimCoal format.
write.2D.fsc<-function(sfs,f.output)
{
  fsc<-as.matrix(sfs)
  pops<-names(dimnames(sfs))
  # add row names and column names
  rownames(fsc)<-paste("d",pops[1],"_",1:nrow(fsc)-1,sep="")
  colnames(fsc)<-paste("d",pops[2],"_",1:ncol(fsc)-1,sep="")
  # output
  cat("1 observations\n\t",file=f.output)
  oldw<-getOption("warn")
  options(warn=-1)
  write.table(fsc,file=f.output,append=T,sep="\t",col.names=T,row.names=T,quote=F)
  options(warn=oldw)
}



########### functions GF ##############
convert_env_trns <- function(path= "./"){ #set path
  file_list <- list.files(path = path,pattern = ".asc$|.bil$|.tif$")
  ras <- paste(path,file_list[1],sep = "")
  raster_final <- raster(ras)
  temp <- values(raster_final)
  tab_final <- data.frame(cell=1:length(temp))
  nb_layers <- length(file_list)
  for(file in file_list){
    nb_layers <- nb_layers-1
    ras <- raster(paste(path,file,sep = ""))
    ras <- values(ras)
    tab_final <- cbind(tab_final,ras)
    names(tab_final)[ncol(tab_final)] <- file
    print(paste("missing layers:",nb_layers))
    
  }
  
  env_trns <- tab_final
  names_layers <- sub(pattern = ".asc$|.bil$|.tif$",replacement = "",x = names(env_trns))
  names(env_trns) <-names_layers 
  env_trns <- env_trns[which(!is.na(env_trns$bio_1)),]
  return(env_trns)
}

euclidian_distance <- function(proj_fut=future,pred_pres=present){
  num <- ncol(proj_fut)
  tot <- rep(0,nrow(proj_fut))
  for(i in 1:num ){
    sum <- (proj_fut[,i]-pred_pres[,i])^2
    tot <- tot+sum
  }
  tot <- sqrt(tot)
  return(tot)
}



pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}


pcaToRaster_two_periods <- function(snps_pres=NA,snps_futuro=NA, rast=NA, mapCells=NA){
  require(raster)
  snps_pres <- data.frame(period="present",snps_pres)
  snps_futuro <- data.frame(period="future",snps_futuro)
  snpPreds <- rbind(snps_pres,snps_futuro)
  pca <- prcomp(snpPreds[,-1], center=TRUE, scale.=FALSE)
  pca <- data.frame(period=snpPreds$period,pca$x)
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  
  ##assigns color to raster
  rast1_p <- rast2_p <- rast3_p <- rast1_f <- rast2_f <- rast3_f <- rast
  
  
  
  rast1_p[mapCells] <- pca[which(pca$period=="present"),"PC1"]
  rast2_p[mapCells] <- pca[which(pca$period=="present"),"PC2"]
  rast3_p[mapCells] <- pca[which(pca$period=="present"),"PC3"]
  
  rast1_f[mapCells] <- pca[which(pca$period=="future"),"PC1"]
  rast2_f[mapCells] <- pca[which(pca$period=="future"),"PC2"]
  rast3_f[mapCells] <- pca[which(pca$period=="future"),"PC3"]
  
  #pca
  pca_present <- pca[which(pca$period=="present"),]
  names(pca_present)[1] <- "cell"
  pca_present$cell <- mapCells
  
  pca_future <- pca[which(pca$period=="future"),]
  names(pca_future)[1] <- "cell"
  pca_future$cell <- mapCells
  
  ##stacks color rasters
  outRast_p <- stack(rast1_p, rast2_p, rast3_p)
  
  names(outRast_p)[1:3]<-c("pca1","pca2","pca3")
  outRast_f <- stack(rast1_f, rast2_f, rast3_f)
  names(outRast_f)[1:3]<-c("pca1","pca2","pca3")
  
  outRast <- list(present=outRast_p,future = outRast_f,
                  pca_present=pca_present,
                  pca_future=pca_future)
  return(outRast)
}


#############################
useful <- function(ras=NA){
  if(is.na(ras)){
    stop("add a path to BIO layer")
  }
  require(maps); require(raster)
  #create a map to plot the data
  mex <- map("world","mexico",plot = F)
  # restart settings 
  opar <- par()
  # create mask for future models 
  ras <-  list.files(path = ras,pattern = ".asc$|.bil$|.tif$")[1]
  mask <- raster(ras)
  mask[mask>0]<- 0  #create a raster mask that will be used to create raster objects in different parts of the script. It is not important which raster layer you use it is only to have the cells of the raster layers.
  return(list(mask=mask,opar=opar,mex=mex))
}

u <- useful("datasets/input/present/")
mask <- u$mask ; mex <- u$mex ; opar <- u$opar


####################################################
################## SDM functions ###################
####################################################

#the next function runs the model by setting the M matrix to the eco-geografic regions where the populations grow, and selecting uncorrelated bio variables the 
run.maxent <- function(path_regions = "datasets/input/official/wwf_terr_ecos.shp",path_pres = "datasets/input/present/",path_input ="datasets/output/maxent_input.csv",genetic_group = "g1",cor.val=0.8,ext= NA,rdp=10000,
                       maxent_jar="~/Escritorio/sdm/datasets/input/maxent.jar",do.ENMeval=FALSE,visible=FALSE,patrn="*.asc$|*.bil$",N=20,plot_pdf=FALSE){
  require("biomod2");require("sp");require("raster");require("dplyr");require("rgdal");require("fuzzySim");require("rgeos");require("magrittr");require("tools");require("readr");require("dismo");require("ENMeval")
  if(!inherits(ext,"Extent")){
    if(!inherits(ext,"list")){
      stop("add extent object or list with x and y extent")
    }else{
      ext <- extent(ext)
    }
  }
  ext<-ext
  #get the ecogeographic regions
  regions <- rgdal::readOGR(path_regions) 
  input <- read.csv(path_input)
  input <- input[input$ID==genetic_group,] ; input$temp <- 1 ; input <- input[,-1]; names(input) <- c("x","y",genetic_group); print(head(input))
  #get bio variables
  pres_clim <- stack(list.files(path=path_pres, pattern = patrn, full.names=TRUE)) %>% crop(x=., y=extent(ext))
  
  outputFolder <- paste("datasets/maxent/",genetic_group,sep = "")
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder)
  }
  #next lines remove the correlated variables
  covarData <- raster::extract(pres_clim, input[,c("x","y")])
  covarData <- cbind(input, covarData[,paste("bio_",1:19,sep = "")])
  correlacion <- corSelect(
    data = covarData,
    sp.cols = grep(genetic_group,names(covarData)),
    var.cols = grep("bio",names(covarData)),
    cor.thresh = cor.val,
    use = "pairwise.complete.obs",
    method="pearson"
  )
  select_var <- correlacion$selected.vars
  
  #crop variables based on the eco-geographic regions
  pres_clim <- pres_clim[[select_var]]
  coords <- input
  crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  sp::coordinates(coords) <- c("x", "y")
  sp::proj4string(coords) <- crs.wgs84
  crs(regions) <- crs(coords)
  
  ecoregionsOfInterest <- sp::over(coords, regions) %>% filter(!is.na(ECO_ID))
  idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
  polygonsOfInterest <- regions[regions$ECO_ID %in% idsEcoRegions, ]
  plot(polygonsOfInterest)
  points(coords,col="red",pch="+")
  pts_b <- gBuffer(coords, width=3)
  pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')
  plot(pts_b, add =T, border = "red")
  polygonsOfInterest <- gIntersection(polygonsOfInterest, pts_b, drop_lower_td = T)
  polygonsOfInterest <- as(polygonsOfInterest, 'SpatialPolygonsDataFrame')
  plot(polygonsOfInterest, add = T, border = "blue")
  writeOGR(polygonsOfInterest, layer = 'ecoregionsOI_M', outputFolder, driver="ESRI Shapefile", overwrite_layer = T)
  
  selectedVariablesCrop <- raster::crop(pres_clim, polygonsOfInterest)
  plot(selectedVariablesCrop,1)
  lines(polygonsOfInterest)
  myExpl <- raster::mask(selectedVariablesCrop,polygonsOfInterest)%>% stack() #Species variables delimited by M
  plot(myExpl)
  plot(myExpl,1); points(coords,pch="+")
  #background
  ext<-extent(myExpl[[1]])
  background <- dismo::randomPoints(myExpl[[1]], n = rdp, excludep = T, coords, ext) %>% data.frame()
  background$temp <-"0" ; names(background)[ncol(background)] <- genetic_group
  plot(myExpl,1,col="lightgrey",legend=F)
  points(background,pch="-",col="blue")
  points(coords,col="red",pch="+")
  
  # Data Formating ####
  presencias<-data.frame(coords)[,1:3] 
  ## BIOMOD####
  DataSpecies<-rbind(presencias,background)
  myResp <- as.numeric(DataSpecies[,genetic_group])
  myRespCoord = DataSpecies[c("x", "y")]
  myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                       expl.var = myExpl,
                                       resp.xy = myRespCoord,
                                       resp.name = genetic_group)
  
  #if do.ENMeval ==T, performs test to estimate characteristics that best fit the models; set the characteristics of the model
  if(do.ENMeval){
    myENMeval <- ENMevaluate(occ=coords@coords,
                             env=myExpl,
                             bg.coords = background[,1:2],
                             RMvalues = c(0.5,1,1.5),
                             clamp=FALSE,
                             algorithm = 'maxnet',
                             method = 'jackknife',
                             parallel = TRUE,
                             numCores = 8)
    write.csv(myENMeval@results[order(myENMeval@results$delta.AICc),],file = paste(outputFolder,"/enmeval_params.csv",sep = ""),row.names = F)
    features <- myENMeval@results[which (myENMeval@results$delta.AICc == 0),]
    linear <- grep("L",features$features); if(length(linear)==1){linear=TRUE}else{linear=FALSE}
    quadratic <- grep("Q",features$features); if(length(quadratic)==1){quadratic=TRUE}else{quadratic=FALSE}
    product <- grep("P",features$features); if(length(product)==1){product=TRUE}else{product=FALSE}
    threshold <- grep("T",features$features); if(length(threshold)==1){threshold=TRUE}else{threshold=FALSE}
    hinge <- grep("H",features$features); if(length(hinge)==1){hinge=TRUE}else{hinge=FALSE}
    
    myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list( path_to_maxent.jar = maxent_jar,
                                                                     memory_allocated = 1024,
                                                                     visible = visible,
                                                                     linear = linear,
                                                                     quadratic = quadratic,
                                                                     product = product,
                                                                     threshold = threshold,
                                                                     hinge = hinge,
                                                                     betamultiplier = features$rm)) 
    
    
  }else{
    myBiomodOption <- BIOMOD_ModelingOptions(MAXENT.Phillips = list( path_to_maxent.jar = maxent_jar,
                                                                     memory_allocated = 1024,
                                                                     visible = visible,
                                                                     linear = TRUE,
                                                                     quadratic = TRUE,
                                                                     product = TRUE,
                                                                     threshold = TRUE,
                                                                     hinge = TRUE,
                                                                     betamultiplier = 1)) 
    
    
  }
  
  
  setwd(outputFolder)
  # run model
  myBiomodModelOut <- BIOMOD_Modeling(
    myBiomodData,
    models = c('MAXENT.Phillips'),
    models.options = myBiomodOption, 
    NbRunEval=N,
    DataSplit =75,
    models.eval.meth = c('TSS','ROC'),
    SaveObj = TRUE,
    rescal.all.models = FALSE,
    do.full.models = FALSE,
    modeling.id = genetic_group)
  
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  write.csv(myBiomodModelEval, file = "myBiomodModelEval.csv",
            row.names = FALSE)
  
  # project model at the present
  myBiomodProj <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = myExpl,
    proj.name = genetic_group,
    selected.models = 'all',
    binary.meth = "TSS",
    compress = 'gzip',
    clamping.mask = FALSE,
    output.format = '.RData')
  
  ### ensemble_modeling#####
  myBiomodEM <- BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('ROC'),
    eval.metric.quality.threshold = c(0.7),
    prob.mean.weight = T)
  
  myBiomodEMEval<-get_evaluations(myBiomodEM)
  write.csv(myBiomodEMEval, file = "myBiomodEMEval.csv",
            row.names = FALSE)
  
  # Creating ensembles projections 
  myBiomodEM_proj <-BIOMOD_EnsembleForecasting(EM.output  = myBiomodEM,
                                               projection.output = myBiomodProj,
                                               selected.models = 'all',
                                               proj.name = genetic_group,
                                               binary.meth = "TSS")
  file_t <- list.files(path = paste(genetic_group,"/proj_",genetic_group,sep=""),pattern = "ensemble_TSSbin.grd",full.names = T)
  currentPred <- stack(file_t)
  
  writeRaster(currentPred, 
              file.path("TSSbin.tif"), 
              suffix='names',
              bylayer=TRUE, 
              overwrite= TRUE)
  plot(currentPred,2)
  if(plot_pdf){
    pdf(paste("projections_",genetic_group,".pdf",sep=""),useDingbats = F)
  }
  plot(myBiomodEM_proj@proj@val@layers[[2]],col=colorRampPalette(colors=c("ivory","olivedrab2","darkolivegreen4","darkgreen","gray9"))(50))
  maps::map("world",add=T)
  points(coords,pch=21,cex=0.8,bg="white")
  title("Present")
  if(plot_pdf){
    dev.off()
  }  
  
  output <- list(coords=coords,
                 env_pres=myExpl,
                 background=background,
                 genetic_group=genetic_group,
                 path_pres=path_pres,
                 myBiomodData=myBiomodData,
                 select_var=select_var,
                 crs=crs.wgs84,
                 polygonsOfInterest=polygonsOfInterest,
                 myBiomodOption=myBiomodOption,
                 myBiomodModelOut=myBiomodModelOut,
                 myBiomodModelEval=myBiomodModelEval,
                 myBiomodProj=myBiomodProj,
                 myBiomodEM=myBiomodEM,
                 myBiomodEMEval=myBiomodEMEval,
                 myBiomodEM_proj=myBiomodEM_proj,
                 currentPred=currentPred)
  return(output)
}


projections.maxent <- function(path_proj,ext,select_var,myBiomodEM,
                               year="year_2050",genetic_group, polygonsOfInterest,
                               coords,patrn="*.asc$|*.bil$|.tif$",pts_b=5,extended=FALSE,clim_stack=FALSE,plot_pdf=FALSE,nam="extended"){
  require("biomod2");require("sp");require("raster");require("dplyr");require("rgdal");require("fuzzySim");require("rgeos");require("magrittr");require("tools");require("readr");require("dismo");require("ENMeval")
  if(clim_stack){
    clim <- stack(list.files(path=path_proj, pattern = patrn, full.names=TRUE))
    names(clim) <- strsplit(names(clim),split = ".",fixed = T) %>% sapply(.,FUN = "[",length(.[[1]])) %>% paste("bio_",.,sep = "")
    clim <- crop(x=clim, y=ext) %>% .[[select_var]]
  }else{
    clim <- stack(list.files(path=path_proj, pattern = patrn, full.names=TRUE)) %>% crop(x=., y=ext) %>% .[[select_var]]
  }
  if(extended){
    pts_b <- gBuffer(coords, width=pts_b)
    pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')
    plot(polygonsOfInterest)
    plot(pts_b, add =T, border = "red")
    selectedVariablesCrop <- raster::crop(clim, pts_b)
    myExpl <- raster::mask(selectedVariablesCrop ,  pts_b) %>% stack()
    plot(myExpl,1)
    plot(polygonsOfInterest,add=T)
  }else{
    pts_b <- gBuffer(coords, width=pts_b)
    pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')
    plot(polygonsOfInterest)
    plot(pts_b, add =T, border = "red")
    polygonsOfInterest <- gIntersection(polygonsOfInterest, pts_b, drop_lower_td = T)
    polygonsOfInterest <- as(polygonsOfInterest, 'SpatialPolygonsDataFrame')
    plot(polygonsOfInterest, add = T, border = "blue")
    selectedVariablesCrop <- raster::crop(clim, polygonsOfInterest)
    myExpl <- raster::mask(selectedVariablesCrop ,  polygonsOfInterest) %>% stack()
    plot(myExpl,1)
    plot(polygonsOfInterest,add=T)
    
  }
  
  myBiomodEM_env <- BIOMOD_EnsembleForecasting(EM.output  = myBiomodEM,
                                               new.env = myExpl,
                                               selected.models = 'all',
                                               proj.name = paste(year,sep = "_"),
                                               binary.meth = "TSS")
  file_t <- list.files(path = paste(genetic_group,"/proj_",year,sep=""),pattern = "ensemble_TSSbin.grd",full.names = T)
  Proj <- stack(file_t)
  plot(Proj,1)
  writeRaster(Proj, 
              file.path(paste(nam,"_proj_",year,"_TSSbin.tif",sep = "")), 
              suffix='names',
              bylayer=TRUE, 
              overwrite= TRUE)
  if(plot_pdf){
    pdf(paste(nam,"_projections_",year,"_",genetic_group,".pdf",sep=""),useDingbats = F)
  }
  plot(myBiomodEM_env@proj@val@layers[[2]],col=colorRampPalette(colors=c("ivory","olivedrab2","darkolivegreen4","darkgreen","gray9"))(50))
  maps::map("world",add=T)
  points(coords,pch=21,cex=0.8,bg="white")
  title(paste(names(coords),nam,year,"Projection"))
  if(plot_pdf){
    dev.off()
  }
  size <- Proj[[2]] %>% rasterToPoints()
  size <- nrow(size[which(size[,3]>0),])
  output <- list(myBiomodEM_env=myBiomodEM_env,Proj=Proj,clim=myExpl,year=year,genetic_group=genetic_group,size=size)
  return(output)
}


create.ras <- function(cluster){
  # List of all sdms in datasets/maxent
  model_list <-  list.files(path="datasets/maxent",pattern=".tif$",full.names = T,recursive = T) %>% .[grep(cluster,.)] %>% .[grep("extended",.)]
  model_list <- model_list[grep("wmeanByROC",model_list)] # we are taking the weighted mean models
  st_o <- lapply(model_list,raster::raster)
  values(st_o[[1]])[values(st_o[[1]]) == 0] <- NA
  values(st_o[[2]])[values(st_o[[2]]) == 0] <- NA
  values(st_o[[3]])[values(st_o[[3]]) == 0] <- NA
  names(st_o) <- sub(".*_proj_","",model_list) %>% strsplit(.,"_TSS") %>% lapply(.,"[",1) %>% unlist()
  return(st_o)
}

##### migration

run.migration <- function(model_gf=NULL,cords_pop=NULL,present=NULL,future=NULL,qt=NULL,qt_cost=0.9,km_dist=10,plot_all=FALSE){
  load("datasets/output/GO_objects.R")
  pops <- rownames(GO_objects$gfData)
  matrix_list <- vector("list",length(pops))
  names(matrix_list) <- pops
  # create a table that will have summary data
  matrix_table <- data.frame(pop=pops,min=NA,median=NA,num_sites=NA,row.names = "pop")
  for(i in 1:nrow(cords_pop)){
    print(rownames(cords_pop)[i])
  # get coordinates of population i
  coord_temp <- coords[pops[i],]
  sp::coordinates(coord_temp) <- c("longitude", "latitude")
  # create a buffer around the population considering the potential areas of migration (M in BAM model)
  pts_b <- rgeos::gBuffer(coord_temp, width=3) %>% as('SpatialPolygonsDataFrame')
  # get the present and future rasters and crop based on the buffer polygon
  matrix_future <- raster::crop(future, pts_b) %>% raster::mask(.,pts_b)%>% stack() %>% rasterToPoints() %>% data.frame()
  coords_furure <- matrix_future[,c("x","y")]
  # predict allelic frequencies in present and future
  matrix_future <- predict(model_gf,matrix_future[,grep("bio",names(matrix_future))])
  matrix_present <- data.frame(raster::extract(present,coord_temp@coords))
  matrix_present<- predict(model_gf,matrix_present[,grep("bio",names(matrix_present))])
  # create a matrix that contains the present bio values of population i repeated as many times as the potential future migration areas (this allows resting matrices for calculating the Euclidian distance.  
  matrix_present <- rbind(matrix_present, matrix_present[rep(1, nrow(matrix_future)-1), ])
  
  #obtain genetic offset between the current location of a population and the future potential areas of migration to obtain the future genetic offset (as in Gougherty et al. 2020)
  gen_off <- (matrix_future-matrix_present)^2
  gen_off <- sqrt(rowSums(gen_off))
  gen_off <- data.frame(coords_furure,gen_off)
  gen_off <- gen_off[which(gen_off$gen_off<=qt),]
  #if there are no migration areas, set that the population will become extinct 
  if(nrow(gen_off)==0){
    matrix_list[[i]] <- "Extinction"
    next
  }
  
  # first we obtain the turnvover function of the adaptive alleles across the candidate bio
  costs <- GO_objects$temp_cand_overall
  costs <- data.frame(bio=costs$x,turnover=costs$y)
  bio_cand <- GO_objects$bio_cand
  # get turnover functions of populations with function predict and relativize by the max value
  # obtain the current candidate raster layer. this will be used as the migration layer
  temp_pres <- raster::crop(present, pts_b) %>% raster::mask(.,pts_b)%>% stack() 
  # get the environmental values of all cells
  matrix_pres <- rasterToPoints(temp_pres) %>% data.frame()
  coords_pres <- matrix_pres[,c("x","y")]
  # predict the allelic turnover across the landscape
  turn_pres <- predict(model_gf,matrix_pres[,grep("bio",names(matrix_pres))])
  # get the value for the candidate value, and create a table that has all the bio values and the GF predicted values
  turn_costs <- data.frame(bio_cand=matrix_pres[,bio_cand],turn=turn_pres[,bio_cand])
  turn_costs <- turn_costs[!duplicated(turn_costs$bio_cand),]
  
  # get the allelic value for the population
  pop_turn<- predict(model_gf,GO_objects$gf_candidate$X[i,])
  pop_turn <- pop_turn[,bio_cand]
  pop <- pops[i]
  coords_temp <- cords_pop[i,]
  
  # for the candidate variable, get the raster layer (migration matrix)
  mig <- temp_pres[[bio_cand]]
  
  # the next loop takes each value of the migration matrix (turn cost biocand values) and transforms it to a migration cost based on the differential between the observed value of the pop and the predicted change across the turnover function
  for(j in 1:nrow(turn_costs)){
    print(j)
    #pop_turn is the value of the pop at its bio_cand distribution; turn_costs[j,"turn"] is the value of the turnover at bio_cand[j]; the differential is the cost of migration by moving from 1 to 2; we take absolute value because there is no order
    set_cost <- abs(pop_turn- turn_costs[j,"turn"])
    # the migration matrix at the bio_cand value j is transformed by the cost created
    mig[mig==turn_costs[j,"bio_cand"]]<- set_cost #for each grid it sets a cost
  }
  
  # remove areas below 0
  mig[mig<0]<-NA
  # transform migration cost from 0 to 1
  mig <- mig/max(values(mig),na.rm = T)
  # the least cost path model uses the reciprocal cost, so we transform it to 1-mig and multiply by 1000, cost will go from 0 to 1000, with 1000 being lower migration costs and plot  
  cost_mig <- 1-mig
  cost_mig <- cost_mig*1000
  
  # create transition object, based on the costs of migration
  tr <- transition(cost_mig,transitionFunction = mean,directions = 8)
  # getall the possible coordinates where the population could settle in the future (remove all 0 areas, not exist)
  future_matrix <- gen_off[,c("x","y")]
  
  # estimate the cost of migration between the focal population and all the potential areas where it could migrate to; and transform it to a vector
  cost <- costDistance(tr,fromCoords = as.matrix(coord_temp@coords),#the coord of the populations
                       toCoords = as.matrix(future_matrix[,c("x","y")]))# all the coord values of settlement
  cost <- as.vector(cost)
  # create table cotaining the coordinates of potential migration and costs
  future_matrix <- data.frame(future_matrix,cost)
  future_matrix <- na.omit(future_matrix)
  
  # if migration is not possible (because there are no future matrix values) the matrix will have no rows, indicate the population will go extinct
  if(nrow(future_matrix)==0){
    matrix_list[[i]] <- "Extinction"
    next
  }
  #get the distance in (km) to each potential area of migration and add it to the future_matrix migration table
  distances <- distGeo(coord_temp@coords[,c("longitude","latitude")],future_matrix[,c("x","y")])/1000
  future_matrix <- data.frame(future_matrix,distances=distances)
  # order the table by increasing distance and add a color with warmer colors indicating lower migration potential (it will need to migrate more).
  temp <- future_matrix[order(future_matrix$distances),]
  par(mfrow=c(2,1))
  mig %>% rasterToPoints() %>% data.frame() -> mig_p
  mig_p$col <- NA
  # order by costs to add colors depending of the genetic offset
  mig_p <- mig_p[order(mig_p$bio_9),]
  colores = c("grey90","grey70","grey45","grey20","grey10")
  mig_p$col <- as.character(cut(1:nrow(mig_p),length(colores),labels=colores))
  colores = maize_pal("MaizMorado")
  temp$col <- rev(as.character(cut(1:nrow(temp),length(colores),labels=colores)))
  
  # plot migration matrix transformed from bio cand to costs of migration
  if(plot_all){
    plot(extent(pts_b),type="n",main=rownames(cords_pop)[i],xlab="lon",ylab="lat")
    maps::map("world","mexico",add=T)
    points(mig_p[,c("x","y")],col=alpha(colour = mig_p$col,alpha = 0.4),pch=15,cex=0.5)
    points(temp[,1:2],pch=15,col=temp$col,cex=0.5)
    points(coord_temp,pch=21,col=maize_pal("HighlandMAGIC")[6],bg="white",lwd=2,cex=1.5)         
    
  }
  # get a matrix ordered by costs of migration and add colors
  temp <- future_matrix[order(future_matrix$cost),]
  colores = maize_pal("MaizMorado")
  temp$col <- rev(as.character(cut(1:nrow(temp),length(colores),labels=colores)))
  if(length(which(cost==Inf))==length(cost)){
    temp <- data.frame(distances=1,cost=1)
  }
 #add the information to the list
  # get the top cost and lower migration distance of populations
  qt_cst <- quantile(future_matrix$cost,qt_cost)
  qt_dist <-quantile(future_matrix$distances,1-qt_cost)
  # if any migration is at a low migration distance but a high cost, remove them, since it is because migration resistance is high
  condition <- which(future_matrix$cost>qt_cst & future_matrix$distances<qt_dist)
  if(length(condition)>0){
    future_matrix <- future_matrix[-condition,]
  }
  matrix_list[[i]]<- future_matrix
  #add summary statistics to the table
  matrix_table[i,"min"]<-min(future_matrix$distances,na.rm=T)
  matrix_table[i,"median"]<-median(future_matrix$distances,na.rm=T)
  matrix_table[i,"num_sites"]<-length(which(future_matrix$distances<km_dist)) #number of areas that are suitable below km_dist Kms
  # } # commented out as part of the loop
  # create an output and return it
  
  }
  output <- list(matrix_list=matrix_list,matrix_table=matrix_table)
  return(output)
}

