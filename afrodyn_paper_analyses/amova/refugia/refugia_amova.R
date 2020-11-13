rm(list = ls())
setwd("~/Dropbox/afrodyn_data/amova/refugia/")
library(raster)
library(maps)
library(GSIF)
library(assigner)
library(radiator)
library("adegenet")
library(codep)
library(adespatial)
library(car)
library(vcfR)
library(elevatr)
library(distancetocoast)
library("rnaturalearth")
library("rnaturalearthdata")
library("hierfstat")
library("pegas")
library(poppr)
library(scales)
library(MASS)

######
# 3rd quantile = df)[ 5]
######

######
# mean = df)[ 4]
######

######
# median = df)[ 3]
######

######
##1st quantile = df)[ 2]
######

sp <-
  c("anni", "anon",  "green", "mona", "podoa", "podob", "sclero")
locs_list <- list()
bin_list <- list()

####
# clim-based refugia
####

load("../../stable_areas/clim_stab.Rdata")

####
# maley refugia
####

#create empty raster
r <- raster(ncol = 5000, nrow = 5000)
crs(r) = '+proj=utm +zone=12 +datum=WGS84'
e <- extent(0, 30, -10, 10)
r <- crop(r, e)

#approximate midpoints of refugia
xy <- c(10, 5)
xy <- rbind(xy, c(11, 3))
xy <- rbind(xy, c(10.5, 1))
xy <- rbind(xy, c(10.5, -2))
xy <- rbind(xy, c(11.5, -1.5))
xy <- rbind(xy, c(12.5, -4.5))

#calculate distance from points for each cell
d2 <- distanceFromPoints(r, xy)

lf <- list.files("../../stable_areas/locs")
lff <-
  list.files("../../stable_areas/locs", full.names = T)
rd <-
  list.files("../../stable_areas/data", full.names = T)

for (i in 1:length(lff)) {
  locs <- read.csv(lff[i])
  locs <- locs[, c("index", "long", "lat")]
  locs <- locs[order(locs$index), ]
  locs$sp <- rep(sp[i], length(locs$long))
  locs_list[[i]] <- locs
}

###
# Read in data
###

load("../vcfs.Rdata")

#genetic data
vcfs <-
  list(anni_vcf,
       anon_vcf,
       green_vcf,
       fries_vcf,
       podo_a_vcf,
       podo_b_vcf,
       sclero_vcf)
names(vcfs) <-
  c("anni", "anon",  "green", "mona", "podo_a", "podo_b", "sclero")


#location data
list_locs <-
  list.files("../../stable_areas/locs/", full.names = T)
list_locs

###
# clim-based refugia
###

#empty lists
list_amova <- list()
list_amova_signif <- list()
list_rand_amova <- list()
list_rand_amova_signif <- list()

par(mfrow = c(1, 3))

for (i in 1:length(vcfs)) {
  pdf(paste(names(vcfs)[i], "_refugia_clim.pdf", sep = ""))
  
  ####
  # enm-based refugia
  ####
  enm_stab <-
    raster(paste("../../enm/stable_areas/", sp[i], "_ccl_tss.tif", sep = ""))
  
  #extract distance from refugia info
  clim_df <- extract(comb_stab, locs_list[[i]][, 2:3])
  enm_df <- extract(enm_stab, locs_list[[i]][, 2:3])
  maley_df <- extract(d2, locs_list[[i]][, 2:3])
  
  ###
  # CALCULATE BINARY THRESHOLDS (CLIM/MALEY)
  ###
  
  bin_df <-
    data.frame(
      locs_list[[i]]$index,
      
      #difference between past/present clim (absolute)
      clim_df < summary(clim_df)[2],
      
      #3 = median, #4 = mean, #5 = 3rd quant
      enm_df,
      
      # distance from refugia centroid
      maley_df < summary(maley_df)[2] #2 = 1st quant, #3 = median, #4 = mean, #5 = 3rd quant
    )
  
  #add binary thresholds to list
  bin_list[[i]] <- bin_df
  
  #convert bin_df to factors
  bin_df <- data.frame(bin_list[[i]])
  colnames(bin_df) <- c("index", "clim", "enm", "maley")
  
  #plot classifications on rasters
  
  plot(comb_stab, xlim = c(0, 20), ylim = c(-10, 10))
  points(locs_list[[i]]$long, locs_list[[i]]$lat, col = as.numeric(as.logical(clim_df < summary(clim_df)[2])) +
           1)
  legend(
    "topright",
    legend = c("non-stable clim", "stable clim"),
    col = c(1, 2),
    pch = 1
  )
  
  plot(enm_stab, xlim = c(0, 20), ylim = c(-10, 10))
  points(locs_list[[i]]$long, locs_list[[i]]$lat, col = (enm_df +
                                                           1))
  legend(
    "topright",
    legend = c("non-stable enm", "stable enm"),
    col = c(1, 2),
    pch = 1
  )
  
  plot(d2)
  points(locs_list[[i]]$long, locs_list[[i]]$lat, col = as.numeric(as.logical(maley_df < summary(maley_df)[2])) +
           1)
  legend(
    "topright",
    legend = c("non-stable maley", "stable maley"),
    col = c(1, 2),
    pch = 1
  )
  
  
  ### convert to genind
  aa.genind <-
    vcfR2genind(vcfs[[i]], return.alleles = TRUE)
  
  if (i == 5) {
    #podo_a has 3 inds with missing coords
    aa.genind <- aa.genind[1:34]
  }
  
  #calculate genetic distances
  distgenEUCL <-
    dist(
      aa.genind,
      method = "euclidean",
      diag = FALSE,
      upper = FALSE,
      p = 2
    )
  
  #add row 2 (clim) as strata in genind
  strata(aa.genind) <- data.frame(bin_df[, 2])
  names(strata(aa.genind)) <- "clim"
  
  #run amova against clim
  #no within-ind calculations
  res_amova <-
    poppr.amova(
      aa.genind,
      ~ clim,
      threads = 1,
      quiet = TRUE,
      within = F
    )
  
  #add amova to list
  list_amova[[i]] <- res_amova
  
  #perform random permutations to get significance
  res_amova_signif   <-
    randtest(res_amova, nrepet = 9999)
  
  #plot hist of random permutations and observed val
  plot(res_amova_signif, main = "list_amova_signif")
  list_amova_signif[[i]] <- res_amova_signif
  
  ###
  #randomized population structure
  ###
  
  #randomly shuffle population assignments and add to strata
  aa.genind.new <- aa.genind
  new_df <-
    data.frame(strata(aa.genind)[sample(nInd(aa.genind)), 1])
  names(new_df) <- "clim"
  strata(aa.genind.new) <- new_df
  
  #run amova
  aa.genind.amova <-
    poppr.amova(
      aa.genind.new,
      ~ clim,
      threads = 1,
      quiet = TRUE,
      within = F
    )
  list_rand_amova[[i]] <- aa.genind.amova
  
  #permute based on this list
  aa.genind.amova.test <-
    randtest(aa.genind.amova, nrepet = 9999)
  list_rand_amova_signif[[i]] <- aa.genind.amova.test
  plot(aa.genind.amova.test, main = "list_rand_amova_signif")
  
  dev.off()
  
}

#extract components of covariance from each AMOVA
for (p in 1:length(list_amova)) {
  foo <-
    cbind(list_amova[[p]]$componentsofcovariance[1:3, ],
          rep(sp[p], 3),
          list_amova_signif[[p]]$pvalue)
  if (p == 1) {
    res_cov <- foo
  } else {
    res_cov <- rbind(res_cov, foo)
  }
}

res_cov_clim <- res_cov
colnames(res_cov_clim) <-
  c("sigma", "percentage", "sp", "p")
write.csv(res_cov_clim, "res_cov_clim.csv")

save.image("clim_refugia.Rdata")


#######################################
# enm-based refugia
#######################################

list_amova <- list()
list_amova_signif <- list()
list_rand_amova <- list()
list_rand_amova_signif <- list()
res_cov <- NULL

#check whether names of genind are in right order compared to location tables

par(mfrow = c(1, 3))

for (i in 1:length(vcfs)) {
  pdf(paste(names(vcfs)[i], "_refugia_enm.pdf", sep = ""))
  
  ### convert to genind
  aa.genind <-
    vcfR2genind(vcfs[[i]], return.alleles = TRUE)
  
  ####
  # enm-based refugia
  ####
  
  enm_stab <-
    raster(paste(
      "~/Dropbox/afrodyn_data/enm/stable_areas/",
      sp[i],
      "_ccl_tss.tif",
      sep = ""
    ))
  
  clim_df <- extract(comb_stab, locs_list[[i]][, 2:3])
  enm_df <- extract(enm_stab, locs_list[[i]][, 2:3])
  maley_df <- extract(d2, locs_list[[i]][, 2:3])
  
  #if all inds are present in 1 category
  #no comparison can be made
  if (length(table(enm_df)) == 1) {
    list_amova[[i]] <- NA
    list_amova_signif[[i]] <- NA
    list_rand_amova[[i]] <- NA
    list_rand_amova_signif[[i]] <- NA
  } else {
    #if (i == 2) {
    #    #anon,podo_a has no inds in ENM refugia, so skipped
    #    #podo_a has all inds in refugia, so skipped
    #    list_amova[[i]] <- NA
    #    list_amova_signif[[i]] <- NA
    #    list_rand_amova[[i]] <- NA
    #    list_rand_amova_signif[[i]] <- NA
    #  } else if (i == 5) {
    #    list_amova[[i]] <- NA
    #    list_amova_signif[[i]] <- NA
    #    list_rand_amova[[i]] <- NA
    #    list_rand_amova_signif[[i]] <- NA
    #  } else {
    #    clim_df <- extract(comb_stab, locs_list[[i]][, 2:3])
    #    enm_df <- extract(enm_stab, locs_list[[i]][, 2:3])
    #    maley_df <- extract(d2, locs_list[[i]][, 2:3])
    
    
    ###
    # CALCULATE BINARY THRESHOLDS (CLIM/MALEY)
    ###
    
    bin_df <-
      data.frame(
        locs_list[[i]]$index,
        
        #difference between past/present clim (absolute)
        clim_df < summary(clim_df)[2],
        
        #3 = median, #4 = mean, #5 = 3rd quant
        enm_df,
        
        # distance from refugia centroid
        maley_df < summary(maley_df)[2] #2 = 1st quant, #3 = median, #4 = mean, #5 = 3rd quant
      )
    
    #add binary thresholds to list
    bin_list[[i]] <- bin_df
    
    #convert bin_df  to factors
    bin_df <- data.frame(bin_list[[i]])
    colnames(bin_df) <-
      c("index", "clim", "enm", "maley")
    
    distgenEUCL <-
      dist(
        aa.genind,
        method = "euclidean",
        diag = FALSE,
        upper = FALSE,
        p = 2
      )
    
    bin_df <- data.frame(bin_list[[i]])
    colnames(bin_df) <-
      c("index", "clim", "enm", "maley")
    
    if (i == 5) {
      aa.genind <- aa.genind[1:34]
    }
    strata(aa.genind) <- data.frame(bin_df[, 3])
    names(strata(aa.genind)) <- "enm"
    
    res_amova <-
      poppr.amova(
        aa.genind,
        ~ enm,
        threads = 1,
        quiet = TRUE,
        within = F
      )
    list_amova[[i]] <- res_amova
    res_amova_signif   <-
      randtest(res_amova, nrepet = 9999)
    plot(res_amova_signif, main = "list_amova_signif")
    list_amova_signif[[i]] <- res_amova_signif
    
    ###
    #randomized population structure
    ###
    
    aa.genind.new <- aa.genind
    new_df <-
      data.frame(strata(aa.genind)[sample(nInd(aa.genind)), 1])
    names(new_df) <- "enm"
    strata(aa.genind.new) <- new_df
    aa.genind.amova <-
      res_amova <-
      poppr.amova(
        aa.genind.new,
        ~ enm,
        threads = 1,
        quiet = TRUE,
        within = F
      )
    list_rand_amova[[i]] <- aa.genind.amova
    
    start_time <- Sys.time()
    aa.genind.amova.test <-
      randtest(aa.genind.amova, nrepet = 9999)
    list_rand_amova_signif[[i]] <- aa.genind.amova.test
    plot(aa.genind.amova.test, main = "list_rand_amova_signif")
    end_time <- Sys.time()
    print(paste("time elapsed:", end_time - start_time))
    
    
  }
  
  dev.off()
  
}

num <- c(1:7)

succ <- num[!is.na(list_amova)]

for (p in num) {
  if (p %in% succ == T) {
    foo <-
      cbind(list_amova[[p]]$componentsofcovariance[1:3, ],
            rep(sp[p], 3),
            list_amova_signif[[p]]$pvalue)
  if (p == 1) {
    res_cov <- foo
  } else {
    res_cov <- rbind(res_cov, foo)
  }
  }
}

res_cov_enm <- res_cov
colnames(res_cov_enm) <-
  c("sigma", "percentage", "sp", "p")
write.csv(res_cov_enm, "res_cov_enm.csv")

save.image("enm_refugia.Rdata")
}


###
# maley-based refugia
###

list_amova <- list()
list_amova_signif <- list()
list_rand_amova <- list()
list_rand_amova_signif <- list()

#check whether names of genind are in right order compared to location tables

par(mfrow = c(1, 3))

for (i in 1:length(vcfs)) {
  pdf(paste(names(vcfs)[i], "_refugia_maley.pdf", sep = ""))
  
  clim_df <- extract(comb_stab, locs_list[[i]][, 2:3])
  enm_df <- extract(enm_stab, locs_list[[i]][, 2:3])
  maley_df <- extract(d2, locs_list[[i]][, 2:3])
  
  ###
  # CALCULATE BINARY THRESHOLDS (CLIM/MALEY)
  ###
  
  bin_df <-
    data.frame(
      locs_list[[i]]$index,
      
      #difference between past/present clim (absolute)
      clim_df < summary(clim_df)[2],
      
      #3 = median, #4 = mean, #5 = 3rd quant
      enm_df,
      
      # distance from refugia centroid
      maley_df < summary(maley_df)[2] #2 = 1st quant, #3 = median, #4 = mean, #5 = 3rd quant
    )
  
  #add binary thresholds to list
  bin_list[[i]] <- bin_df
  
  #convert bin_df to factors
  bin_df <- data.frame(bin_list[[i]])
  colnames(bin_df) <- c("index", "clim", "enm", "maley")
  
  ### convert to genind
  aa.genind <-
    vcfR2genind(vcfs[[i]], return.alleles = TRUE)
  
  if (i == 5) {
    #podo_a has 3 inds with missing coords
    aa.genind <- aa.genind[1:34]
  }
  
  distgenEUCL <-
    dist(
      aa.genind,
      method = "euclidean",
      diag = FALSE,
      upper = FALSE,
      p = 2
    )
  
  bin_df <- data.frame(bin_list[[i]])
  colnames(bin_df) <- c("index", "clim", "enm", "maley")
  
  strata(aa.genind) <- data.frame(bin_df[, 4])
  names(strata(aa.genind)) <- "maley"
  
  res_amova <-
    poppr.amova(
      aa.genind,
      ~ maley,
      threads = 1,
      quiet = TRUE,
      within = F
    )
  list_amova[[i]] <- res_amova
  
  res_amova_signif   <-
    randtest(res_amova, nrepet = 9999)
  plot(res_amova_signif, main = "list_amova_signif")
  list_amova_signif[[i]] <- res_amova_signif
  
  ###
  #randomized population structure
  ###
  
  aa.genind.new <- aa.genind
  new_df <-
    data.frame(strata(aa.genind)[sample(nInd(aa.genind)), 1])
  names(new_df) <- "maley"
  strata(aa.genind.new) <- new_df
  aa.genind.amova <-
    poppr.amova(
      aa.genind.new,
      ~ maley,
      threads = 1,
      quiet = TRUE,
      within = F
    )
  list_rand_amova[[i]] <- aa.genind.amova
  
  aa.genind.amova.test <-
    randtest(aa.genind.amova, nrepet = 9999)
  list_rand_amova_signif[[i]] <- aa.genind.amova.test
  plot(aa.genind.amova.test, main = "list_rand_amova_signif")
  
  
  dev.off()
}



for (p in 1:length(list_amova)) {
  foo <-
    cbind(list_amova[[p]]$componentsofcovariance[1:3, ],
          rep(sp[p], 3),
          list_amova_signif[[p]]$pvalue)
  if (p == 1) {
    res_cov <- foo
  } else {
    res_cov <- rbind(res_cov, foo)
  }
}

res_cov_maley <- res_cov
colnames(res_cov_maley) <-
  c("sigma", "percentage", "sp", "p")
write.csv(res_cov_maley, "res_cov_maley.csv")

save.image("maley_refugia.Rdata")


###############################
# GROUPED BARPLOT
###############################


bet_maley <-
  res_cov_maley[grep("Between", rownames(res_cov_maley)),]
bet <- rbind(res_cov_clim, res_cov_enm, res_cov_maley)

ref_type <- c(rep("clim", length(res_cov_clim$sigma)),
              rep("enm", length(res_cov_enm$sigma)),
              rep("maley", length(res_cov_maley$sigma)))

ref_type_bet <- cbind(ref_type, bet)

ref_type_bet <-
  ref_type_bet[grep("Between", rownames(ref_type_bet)),]

# library
library(ggplot2)

# Grouped
ggplot(ref_type_bet, aes(fill = sp, y = percentage, x = ref_type)) +
  geom_bar(position = "dodge", stat = "identity") + xlab("Model type")  + ylab("% variance explained by refugia")
