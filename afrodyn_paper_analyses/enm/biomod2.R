#!/usr/bin/env Rscript

rm(list = ls())

args = commandArgs(trailingOnly = TRUE)

###
# Biomod2
###
setwd("~/Dropbox/projects/AJH_AFRODYN/enm/")
library(biomod2)
library(ConR)
library(maps)

###
# Environmental data
###
#fili<-list.files("~/Téléchargements",pattern="*tif",full.names = T)
#fili_tif<-list.files("~/Téléchargements",pattern="*tif")
#for(i in 1:length(fili)){
#  foo<-raster(fili[i])
#  e<-extent(predictors[[1]])
#  foo<-crop(foo,e)
#  foo<-resample(foo, predictors[[1]])
#  writeRaster(foo,fili_tif[i],format='GTiff')
#}

soil_tif <- list.files(path = "./soil/",
                       pattern = "*tif",
                       full.names = T)
soil_files = lapply(soil_tif, raster)

#remove categorical var TAXNWRB (soil type)
soil_files <- soil_files[1:(length(soil_files) - 1)]
soil_files

#present-day climate data
temp = list.files(path = "cropped/wc2/",
                  pattern = "*.tif",
                  full.names = T)
myfiles = lapply(temp, raster)


#get elevation data
#library(elevatr)
#locs_list = list.files(path="sp_locs",full.names = T)
#
#all_loc<-lapply(locs_list,read.csv)
#
#all_loc[[1]]$collector<-rep("NA",length(all_loc[[1]]$species))
#colnames(all_loc[[1]])<-colnames(all_loc[[2]])
#
#
#all_loc_bind<-rbind(all_loc[[1]],
#      all_loc[[2]],
#      all_loc[[3]],
#      all_loc[[4]],
#      all_loc[[5]],
#      all_loc[[6]],
#      all_loc[[7]])

#ele<-get_elev_raster(locations=myfiles[[1]],prj=sp::proj4string(myfiles[[1]]),z=5)

#saveRDS(ele,"ele.Rdata")

#load in elevation data
ele <- readRDS("ele.Rdata")

#make elevation same resolution as other rasters
ele <- resample(ele, myfiles[[1]])

#crop to extent of clim rasters
e <- extent(myfiles[[1]])
ele <- crop(ele, e)
names(ele) <- "alt"

#crop soil to extent of clim rasters
soil_files <- crop(stack(soil_files), e)

#combine bioclim, altitude and soil vars 
all_pres_vars <- c(myfiles, ele, soil_files)

#make present stack
pres <- stack(all_pres_vars)

########
#Checked correlation 24/4/20
########

#checked correlation with ENMTools pearson > 0.9 removed
#same variables removed in past climate layers

library(ENMTools)
#rcm <- raster.cor.matrix(pres)
#which(rcm > 0.9, arr.ind = T)

#check no corrs after dropping layers
#rcm <- raster.cor.matrix(predictors)
#which(rcm > 0.9, arr.ind = T)

library(biomod2)

#drop correlates
predictors <-
  dropLayer(pres,
            c(
              "bio_7",
              "bio_10",
              "bio_11",
              "bio_16",
              "bio_17",
              "PHIHOX_M_sl2_250m",
              "TAXNWRB_250m"
            ))

predictors <- stack(predictors)

###
#past climate data
###

#ccl
temp = list.files(path = "cropped/ccl/",
                  pattern = "*.tif",
                  full.names = T)
cclfiles = lapply(temp, raster)
srtm <- c(cclfiles, ele, soil_files)
srtm <- stack(srtm)

#remove correlates
ccl_predictors <-
  dropLayer(
    srtm,
    c(
      "bio_7",
      "bio_10",
      "bio_11",
      "bio_16",
      "bio_17",
      "PHIHOX_M_sl2_250m",
      "TAXNWRB_250m"
    )
  )

#past climate data
#mel
temp = list.files(path = "cropped/mel/",
                  pattern = "*.tif",
                  full.names = T)
melfiles = lapply(temp, raster)
srtm <- c(melfiles, ele, soil_files)
srtm <- stack(srtm)

#remove correlates
mel_predictors <-
  dropLayer(
    srtm,
    c(
      "bio_7",
      "bio_10",
      "bio_11",
      "bio_16",
      "bio_17",
      "PHIHOX_M_sl2_250m",
      "TAXNWRB_250m"
    )
  )

#past climate data
#mrl
temp = list.files(path = "cropped/mrl/",
                  pattern = "*.tif",
                  full.names = T)
mrlfiles = lapply(temp, raster)
srtm <- c(mrlfiles, ele, soil_files)
srtm <- stack(srtm)

#remove correlates
mrl_predictors <-
  dropLayer(
    srtm,
    c(
      "bio_7",
      "bio_10",
      "bio_11",
      "bio_16",
      "bio_17",
      "PHIHOX_M_sl2_250m",
      "TAXNWRB_250m"
    )
  )

locs_list = list.files(path = "sp_locs", full.names = T)
sp <- c("anni",
        "anon",
        "mona",
        "green",
        "podoa",
        "podob",
        "sclero")



#rename
names(predictors)
names(ccl_predictors)
names(mrl_predictors)
names(mel_predictors)

ccl_predictors <- stack(ccl_predictors)
mrl_predictors <- stack(mrl_predictors)
mel_predictors <- stack(mel_predictors)

library(biomod2)

#names
sp <- c("anni", "anon", "mona", "green", "podoa", "podob", "sclero")

####
# loop starts
####

#interactive loop to be used for batch runs

for (i in as.numeric(args[1])) {
  print("Running species:")
  print(sp[as.numeric(args[1])])
    
#non-interactive

#for (i in 1) {
#  print("Running species:")
#  print(sp[i])

  locs <- read.csv(locs_list[i])
  locs <- locs[, c(2:3)]
  head(locs)
  colnames(locs)[1:2] <- c("Longitude", "Latitude")
  
  
  ###
  # Convex hull
  ###
  
  #need to reverse lat and lon
  conv.hull <- EOO.computing(locs[, 2:1], export_shp = T)
  
  #####
  
  #5 degree barrier
  conv.hull_1 <-
    raster::buffer(x = conv.hull[[2]], width = 5) ### Buffering the hull convex
  crs(conv.hull_1) <- crs(predictors[[1]])
  
  #### Masking
  masked_clim.var <- raster::mask(predictors[[1]], conv.hull_1)
  
  
  map('world', col = 'gray90', fill = TRUE,xlim=c(-30,40),ylim=c(-10,30))
  plot(conv.hull_1, add = T)
  
  #mask layers to convex hull extent
  predictors_mask <- mask(predictors, conv.hull_1)
  predictors_mask <- stack(predictors_mask)
  mrl_predictors_mask <- mask(mrl_predictors, conv.hull_1)
  mrl_predictors_mask <- stack(mrl_predictors_mask)
  ccl_predictors_mask <- mask(ccl_predictors, conv.hull_1)
  ccl_predictors_mask <- stack(ccl_predictors_mask)
  mel_predictors_mask <- mask(mel_predictors, conv.hull_1)
  mel_predictors_mask <- stack(mel_predictors_mask)
  
  ##' 1. Formatting Data
  myBiomodData <-
    BIOMOD_FormatingData(
      
      #number of presences
      resp.var = rep(1, length(locs$Longitude)),
      
      #coordinates
      resp.xy = locs,
      
      #predictors
      expl.var = predictors_mask,
      
      #pseudo absences rep
      PA.nb.rep = 100,
      
      #number of pseudo-absence selected for each repetition (when PA.nb.rep > 0) of the selection
      PA.nb.absences = 1000,
      
      #strategy for selecting the Pseudo Absences
      PA.strategy = "sre",
      
      #quantile used for ‘sre’ Pseudo Absences selection
      PA.sre.quant = 0.05,
      
      #response variable name (character). The species name.
      resp.name = sp[i]
    )
  
  ###
  # Defining Models Options using default options.
  ###
  
  myBiomodOption <- BIOMOD_ModelingOptions()
  
  ###
  # Modelling
  ###
  
  #This function allows to calibrate and evaluate a range of species distribution models techniques run over a given species.
  #Calibrations are made on the whole sample or a random subpart.
  #The predictive power of the different models is estimated using a range of evaluation metrics.
  
  myBiomodModelOut <- BIOMOD_Modeling(
    
    #data
    myBiomodData,
    
    #models to be used
    models = c('GLM',
               'GBM',
               'ANN',
               'RF',
               'MAXENT.Phillips'),
    
    #default options
    models.options = myBiomodOption,
    
    #Number of Evaluations to run
    NbRunEval = 10,
    
    #% of data used to calibrate the models, the remaining part will be used for testing
    DataSplit = 70,
    
    #response points weights
    Yweights = NULL,
    
    #Number of permutation to estimate variable importance
    VarImport = 10,
    
    #vector of names of evaluation metric
    models.eval.meth = c('KAPPA', 'TSS'),
    SaveObj = TRUE,
    
    #if true, all model prediction will be scaled with a binomial GLM
    rescal.all.models = FALSE,
    
    #if true, models calibrated and evaluated with the whole dataset are done
    do.full.models = FALSE
  )
  
  ##' print a summary of modeling stuff
  get_evaluations(myBiomodModelOut)
  
  ####
  #Ensemble Modelling
  ####
  
  #BIOMOD_EnsembleModeling combines models and make ensemble predictions built with BIOMOD_Modeling.
  #The ensemble predictions can also be evaluated against the original data given to BIOMOD_Modeling. 
  #Biomod2 proposes a range of options to build ensemble models and predictions and to assess the modeling uncertainty.
  #The created ensemble models can then be used to project distributions over space and time as classical biomod2 models.

  myBiomodEM <-
    BIOMOD_EnsembleModeling(
      
      #a "BIOMOD.models.out" returned by BIOMOD_Modeling
      modeling.output = myBiomodModelOut,
      
      #a character vector (either 'all' or a sub-selection of model names) that defines the models
      #kept for building the ensemble model
      chosen.models = 'all',
      
      #Character. Flag defining the way the models will be combined to build the ensemble models.
      em.by = 'all',
      
      #vector of names of evaluation metric used to build ensemble models.
      eval.metric = c('TSS'),
      
      #the minimum scores below which models will be excluded of the ensemble-models building.
      eval.metric.quality.threshold = c(0.6),
      
      #the evaluation methods used to evaluate ensemble models
      models.eval.meth = c('KAPPA', 'TSS'),
      
      #Logical. Estimate the mean probabilities across predictions
      prob.mean = TRUE,
      
      #Logical. Estimate the coefficient of variation across predictions
      prob.cv = FALSE,
      
      #Logical . Estimate the confidence interval around the prob.mean
      prob.ci = FALSE,
      
      #Numeric. Significance level for estimating the confidence interval.
      prob.ci.alpha = 0.05,
      
      #Logical. Estimate the mediane of probabilities
      prob.median = FALSE,
      
      #Logical. Estimate the committee averaging across predictions
      committee.averaging = FALSE,
      
      #Logical. Estimate the weighted sum of probabilities
      prob.mean.weight = TRUE,
      
      #Define the relative importance of the weights. 
      #A high value will strongly discriminate the 'good' models from the 'bad' ones (see the details section). 
      #If the value of this parameter is set to 'proportional' (default), then the attributed weights
      #are proportional to the evaluation scores given by 'weight.method'(eval.metric)
      prob.mean.weight.decay = 'proportional'
    )
  
  # print summary
  myBiomodEM
  
  # get evaluation scores
  get_evaluations(myBiomodEM)
  
  ###
  # PROJECTION and ensemble FORECASTING
  ###
  
  #biomod2 is able to project potential distributions of species in other areas, other resolutions or other time scales.

  ##present projection
  presBiomodProjectionFuture <-
    BIOMOD_Projection(
      
      #"BIOMOD.models.out" object produced by a BIOMOD_Modeling run
      modeling.output = myBiomodModelOut,
      
      #A set of explanatory variables onto which models will be projected
      new.env = predictors_mask,
      
      proj.name = 'pres',
      
      #'all' when all models have to be used to render projections
      selected.models = 'all',
      
      binary.meth = NULL,
      compress = FALSE,
      
      #if TRUE, a clamping mask will be saved on hard drive 
      build.clamping.mask = TRUE
    )
  
  #This function use projections of ‘individual models’ and ensemble models from BIOMOD_EnsembleModeling
  #to build an ensemble of species' projections over space and time.
  
  ##present ensemble
  presBinBIOMOD_EnsembleForecasting <-
    BIOMOD_EnsembleForecasting(
      
      #a "BIOMOD.EnsembleModeling.out" returned by BIOMOD_EnsembleModeling
      myBiomodEM,
      
      #a "BIOMOD.projection.out" returned by BIOMOD_Projection
      projection.output = presBiomodProjectionFuture,
      
      #if not 'all', a character vector containing a subset of ensemble-models you want make projacion
      selected.models = 'all',
      
      proj.name = "ensemble_pres_bin",
      
      #vector specifying the names of evaluation metrics and associated thresholds to transform the probabilities of
      #presence into presence and absence (binary transformation).
      binary.meth = c('KAPPA', 'TSS'),
      
      filtered.meth = NULL,
      compress = TRUE
    )
  
  ##mel
  melBiomodProjectionFuture <-
    BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = mel_predictors_mask,
      proj.name = 'mel',
      selected.models = 'all',
      binary.meth = NULL,
      compress = FALSE,
      build.clamping.mask = TRUE
    )
  
  melBinBIOMOD_EnsembleForecasting <-
    BIOMOD_EnsembleForecasting(
      myBiomodEM,
      projection.output = melBiomodProjectionFuture,
      selected.models = 'all',
      proj.name = "ensemble_mel_bin",
      binary.meth = c('KAPPA', 'TSS'),
      filtered.meth = NULL,
      compress = TRUE
    )
  
  
  ##ccl
  cclBiomodProjectionFuture <-
    BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = ccl_predictors_mask,
      proj.name = 'ccl',
      selected.models = 'all',
      binary.meth = NULL,
      compress = FALSE,
      build.clamping.mask = TRUE
    )
  
  
  cclBinBIOMOD_EnsembleForecasting <-
    BIOMOD_EnsembleForecasting(
      myBiomodEM,
      projection.output = cclBiomodProjectionFuture,
      selected.models = 'all',
      proj.name = "ensemble_ccl_bin",
      binary.meth = c('KAPPA', 'TSS'),
      filtered.meth = NULL,
      compress = TRUE
    )
  
  ##mrl
  mrlBiomodProjectionFuture <-
    BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = mrl_predictors_mask,
      proj.name = 'mrl',
      selected.models = 'all',
      binary.meth = NULL,
      compress = FALSE,
      build.clamping.mask = TRUE
    )
  
  
  mrlBinBIOMOD_EnsembleForecasting <-
    BIOMOD_EnsembleForecasting(
      myBiomodEM,
      projection.output = mrlBiomodProjectionFuture,
      selected.models = 'all',
      proj.name = "ensemble_mrl_bin",
      binary.meth = c('KAPPA', 'TSS'),
      filtered.meth = NULL,
      compress = TRUE
    )

}
