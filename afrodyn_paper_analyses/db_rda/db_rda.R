        #' ---
        #' title: "db-RDA and db-MEM on AFRODYN species"
        #' author: "Andrew J. Helmstetter"
        #' date: "April 9th, 2020"
        #' output:
        #'    html_document:
        #'      toc: true
        #'      highlight: haddock
        #' ---
        #'
        #+ messages = "hide"
        #+ message = FALSE
        #+ messages = FALSE
        #+ results = "hide"
        #+ include = F
        rm(list = ls())
        setwd("~/Dropbox/afrodyn_data/db_rda/")
        
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
        library(raster)
        library(maps)
        library(GSIF)
        library("rnaturalearth")
        library("rnaturalearthdata")
        library("hierfstat")
        library("pegas")
        library(poppr)
        library(vcfR)
        library(scales)
        par(mar = c(3, 3, 3, 3))
        
        #+ fig.height = 10, fig.width = 10
        
        #https://github.com/laurabenestan/db-RDA-and-db-MEM/blob/master/README.md
        #https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html
        par(mar = c(2, 2, 2, 2))
        
        #create empty raster
        r <- raster(ncol = 5000, nrow = 5000)
        crs(r) = '+proj=utm +zone=12 +datum=WGS84'
        r <- crop(r, e)
        
        #approximate midpoints of refugia
        xy <- c(10, 5)
        xy <- rbind(xy, c(11, 3))
        xy <- rbind(xy, c(10.5, 1))
        xy <- rbind(xy, c(10.5, -2))
        xy <- rbind(xy, c(11.5, -1.5))
        xy <- rbind(xy, c(12.5, -4.5))
        xy
        
        #world shape for mask
        world <- ne_countries(scale = "medium", returnclass = "sf")
        
        #calculate distance from points for each cell
        d2 <- distanceFromPoints(r, xy)
        
        lf <- list.files("../stable_areas/locs")
        lff <- list.files("../stable_areas/locs", full.names = T)
        rd <- list.files("../stable_areas/data", full.names = T)
        
        temp = list.files(path = "~/Dropbox/afrodyn_data/enm/cropped/wc2/",
                          pattern = "*.tif",
                          full.names = T)
        currclim = lapply(temp, raster)
        currclim <- stack(currclim)
        
        #checked correlation with ENMTools pearson > 0.9 removed
        currclim <-
          dropLayer(currclim,
                    c(
                      "bio_7",
#                      "bio_8",
#                      "bio_9",
                      "bio_10",
                      "bio_11",
                      "bio_16",
                      "bio_17"
                    ))
        
        temp = list.files(path = "~/Dropbox/afrodyn_data/enm/cropped/ccl/",
                          pattern = "*.tif",
                          full.names = T)
        cclclim = lapply(temp, raster)
        cclclim <- stack(cclclim)
        
        #checked correlation with ENMTools pearson > 0.9 removed
        cclclim <-
          dropLayer(cclclim,
                    c(
                      "bio_7",
#                      "bio_8",
#                      "bio_9",
                      "bio_10",
                      "bio_11",
                      "bio_16",
                      "bio_17"
                    ))
        
        ###
        #soil order
        ###
        
        library(GSIF)
        
        #extract soil suborder
        data("landmask20km")
        soil_so_raw <- raster(landmask20km['suborder'])
        soil_so <-
          projectRaster(soil_so_raw, crs = "+proj=longlat +datum=WGS84", method = "ngb")
        levels(soil_so_raw)
        
        # convert suborder to order
        old <- 1:64
        new <- c(
          1:4,
          rep(5, 6),
          rep(11, 5),
          rep(16, 7),
          rep(23, 5),
          rep(28, 6),
          rep(34, 6),
          rep(40, 5),
          rep(45, 8),
          rep(53, 5),
          rep(58, 6)
        )
        
        values(soil_so)[values(soil_so) %in% old] <-
          new[match(values(soil_so), old, nomatch = 0)]
        values(soil_so)[values(soil_so) %in% old] <-
          new[match(values(soil_so), old, nomatch = 0)]
        table(values(soil_so))
        
        ##
        # soil characteristics
        ##
        
        temp = list.files(path = "~/Dropbox/afrodyn_data/enm/soil/",
                          pattern = "*.tif",
                          full.names = T)
        currsoil = lapply(temp, raster)
        currsoil <- stack(currsoil)
        currsoil <- dropLayer(currsoil, c("PHIHOX_M_sl2_250m", "TAXNWRB_250m"))
        
        ###
        #stability files
        ###
        
        load("../stable_areas/clim_stab.Rdata")
        
        ###
        # Geographic locations
        ###
        
        #species names placeholders
        spe <- c("anni", "anon", "green", "mona", "podo_a", "podo_b", "sclero")
        sp <- c("anni", "anon", "green", "mona", "podoa", "podob", "sclero")
        
        library(vegan)
        
        #list placeholders
        gen_list <- list()
        distgenEUCL_list <- list()
        locs_list <- list()
        clim_list <- list()
        soil_list <- list()
        
        #read in location data
        lf <- list.files("../stable_areas/locs")
        lff <- list.files("../stable_areas/locs", full.names = T)
        
        #read in genetic data
        rd <- list.files("../stable_areas/data", full.names = T)
        
        for (i in 1:length(lff)) {
          
          #combine location data
          locs <- read.csv(lff[i])
          locs <- locs[, c("index", "long", "lat")]
          locs <- locs[order(locs$index),]
          locs$sp <- rep(spe[i], length(locs$long))
          if (i == 1) {
            all_locs <- locs
            locs_list[[i]] <- locs
          } else {
            all_locs <- rbind(all_locs, locs)
            locs_list[[i]] <- locs
          }
          
          #read in table with index/diversity stat
          
          if (i == 5) {
            #podo has 3 inds with no coords
            load(rd[i])
            gen <- vcfR2genind(vcf,return.alleles=T)
            gen <- gen[1:(length(indNames(gen)) - 3)]
            pop(gen) <- indNames(gen)
            gp <- genind2genpop(gen)
            gen_list[[i]] <- gp
            
          } else {
            load(rd[i])
            gen <- vcfR2genind(vcf,return.alleles=T)
            pop(gen) <- indNames(gen)
            gp <- genind2genpop(gen)
            gen_list[[i]] <- gp
            
          }
          
          #First, estimating individual genetic distances is a crucial step before performing the db-RDA.
          #The individual genetic distances will be considered as the response variables.
          distgenEUCL_list[[i]] <-
            dist(
              gen_list[[i]],
              method = "euclidean",
              diag = FALSE,
              upper = FALSE,
              p = 2
            )
          
          #Second, import the environmental dataset, which will contain the explanatory variables
          #extract data from climate dataset based on coordinates
          #put in list
          
          for (j in 1:length(names(currclim))) {
            if (j == 1) {
              
        
              clim_list[[i]] <- raster::extract(currclim[[j]], locs_list[[i]][, 2:3])
              
            } else {
              clim_list[[i]] <-
                cbind(clim_list[[i]], raster::extract(currclim[[j]], locs_list[[i]][, 2:3]))
              
            }
            
          }
          
          #extract data from soil dataset based on coordinates
          #put in list
          
          for (j in 1:length(names(currsoil))) {
            if (j == 1) {
              soil_list[[i]] <- raster::extract(currsoil[[j]], locs_list[[i]][, 2:3])
              
            } else {
              soil_list[[i]] <-
                cbind(soil_list[[i]], raster::extract(currsoil[[j]], locs_list[[i]][, 2:3]))
              
            }
            
          }
          
          #get elevation data
          locs_for_ele <- locs_list[[i]][, 2:3]
          colnames(locs_for_ele) <- c("x", "y")
          ele <-
            get_elev_raster(
              locations = locs_list[[i]][, 2:3],
              prj = sp::proj4string(currclim[[1]]),
              z = 5
            )
          
          #raster::extract and bind to climate list
          clim_list[[i]] <-
            cbind(clim_list[[i]], raster::extract(ele, locs_list[[i]][, 2:3]))
          
          ###
          # STABILITY ENM
          ###
          
          
          #read in ENM data
          enm_stab <-
            raster(paste("../enm/stable_areas/", sp[i], "_mel_tss.tif", sep = ""))
        
          #raster::extract ENM stability data
          
          #skip ENM stab if all coordinates 0 (not in stable area)
          if (length(unique(raster::extract(enm_stab, locs_list[[i]][, 2:3]))) > 1) {
            clim_list[[i]] <-
              cbind(clim_list[[i]],
                    raster::extract(ann_temp, locs_list[[i]][, 2:3]),
                    raster::extract(ann_prec, locs_list[[i]][, 2:3]),
                    raster::extract(enm_stab, locs_list[[i]][, 2:3], ))
            
          } else {
            clim_list[[i]] <-
              cbind(clim_list[[i]],
                    raster::extract(ann_temp, locs_list[[i]][, 2:3]),
                    raster::extract(ann_prec, locs_list[[i]][, 2:3]))
            
          }
          
          #get soil data (suborder)
          #old = suborder numbers
          #new = order numbers
          
          #swap in order numbers for suborder numbers in raster
          values(soil_so)[values(soil_so) %in% old] <-
            new[match(values(soil_so), old, nomatch = 0)]
          
          #raster::extract order from raster
          so <- data.frame(raster::extract(soil_so, locs_list[[i]][, 2:3]))
          colnames(so) <- "number"
          
          #get unique orders
          so_uniq <- na.omit(unique(so$number))
          
          #reformat data frame to one column 0/1 per soil order
          for (g in 1:length(so_uniq)) {
            so[, g + 1] <- so$number == so_uniq[g]
            so[, g + 1] <- as.integer(as.logical(so[, g + 1]))
            
            soil_names <- vector()
            for (h in 1:length(so_uniq)) {
              soil_names[h] <-
                as.character(levels(soil_so_raw)[[1]]$levels[so_uniq[h]][1])
              
            }
          }
          
          so <- so[, c(2:length(so[1, ]))]
          so[is.na(so)] <- 0
          
          #get distance from refugia data
          ref_dist <- raster::extract(d2, locs_list[[i]][, 2:3])
          
          #get distance from coast data
          #coast
          coast_dist <-
            raster::extract(distance_to_coastline_lowres, locs_list[[i]][, 2:3])
          #in km
          coast_dist <- coast_dist / 1000
          
          
          #combine all data into 1 list
          clim_list[[i]] <- cbind(clim_list[[i]], so, ref_dist, coast_dist)
          clim_list[[i]] <- cbind(clim_list[[i]], soil_list[[i]])
          
          
          #give names to each variable
          #contingent on whether enm_stab is present for species
          if (length(unique(extract(enm_stab, locs_list[[i]][, 2:3]))) == 1) {
            colnames(clim_list[[i]]) <-
              c(
                names(currclim),
                "alt",
                "clim_stab_temp",
                "clim_stab_prec",
                soil_names,
                "refugia_dist",
                "coast_dist",
                names(currsoil)
              )
          } else {
            colnames(clim_list[[i]]) <-
              c(
                names(currclim),
                "alt",
                "clim_stab_temp",
                "clim_stab_prec",
                "enm_stab",
                soil_names,
                "refugia_dist",
                "coast_dist",
                names(currsoil)
              )
            
          }
          
        }
      
        #check in same order
      cbind(row.names(as.matrix(distgenEUCL_list[[7]], labels = TRUE)), as.character(locs_list[[7]]$index))
      
      #Look at sites in space by keeping only latitude and longitude and saving this in the object called Coor.
      #Keep latitude and longitude in this order as for the function gcd.hf, latitude needs to be first.
      
      Coor_list <- list()
      Coorxy_list <- list()
      DistSpatial_list <- list()
      dbmem_list <- list()
      
      #'
      #' Locations
      #'
      
      par(mar = c(4, 4, 4, 4))
      for (i in 1:length(locs_list)) {
        Coor_list[[i]] = locs_list[[i]][, 3:2]
        Coorxy_list[[i]] = locs_list[[i]][, 2:3]
      
        #Compute spatial distances among sites accounting for the earth curvature.
        #dist also calculates euclidean distances
        DistSpatial_list[[i]] <- gcd.hf(Coor_list[[i]])
        
        #Compute MEM by keeping default setting for truncation
        #(length of the longest edge of the minimum spanning tree will be used as the threshold) and just positive MEM.
        
        ##
        # warnings due to gcd.hf
        ##
        
        dbmem_list[[i]] <- dbmem(DistSpatial_list[[i]])
        
        ##look at output
        summary(dbmem_list[[i]])
      }
      
      par(mfrow = c(1, 1))
      library(adegraphics)
      
      
      #Visualising the 19 MEMs by using the Coorxy object.
      #The 1rst dbmems are large spatial scales, the last dbmems are small spatial scales.
      #dbmem can be used in stats like any other env variables.
      for (i in 1:length(locs_list)) {
        adegraphics::s.label(Coor_list[[i]],
                             nb = attr(dbmem_list[[i]], "listw"),
                             main = spe[i])
        
      
        ade4::s.value(Coorxy_list[[i]], dbmem_list[[i]][, 1])
      }
      
      
      #PCOA
      
      #Perform a PCoA on genetic distance matrix.
      #The advantage of the PCoA is that it does not take into account the missing observations
      #(corresponding to 0 in the genomic data frame) in the response variable comparatively to the PCA.
      
      Pcoa_list <- list()
      Y_list <- list()
      X_list <- list()
      
      for (i in 1:length(locs_list)) {
        
        #pcoa on genetic distances
        Pcoa_list[[i]] <- ape::pcoa(distgenEUCL_list[[i]])
        Pcoa_list[[i]]$correction
        Pcoa_list[[i]]$note
        head(Pcoa_list[[i]]$values)
        
        #Genetic distances are then in adequate multivariate space format for running the db-RDA.
        #Extract Pcoa principal components, which will be the response variable in the db-RDA.
        X_list[[i]] = Pcoa_list[[i]]$vectors
        
        #Look at genotypes distribution in relation to the first two PCoA axes.
        plot(X_list[[i]][, 1], X_list[[i]][, 2], main = spe[i])
        
        #Create a matrix with all expanatory variables, 19 MEMs
        Y_list[[i]] <- cbind(dbmem_list[[i]], clim_list[[i]])
        
      }
      
      
      #look at correlation
      #MEMs are orthogonal to each other so correlation unlikely
      #will need to remove correllated vars
      cor(Y_list[[1]])>0.9
      
      write.csv(cor(Y_list[[1]])>0.9,"anni_cor.csv")
      write.csv(cor(Y_list[[2]])>0.9,"anon_cor.csv")
      write.csv(cor(Y_list[[3]])>0.9,"green_cor.csv")
      write.csv(cor(Y_list[[4]])>0.9,"mona_cor.csv")  
      write.csv(cor(Y_list[[5]])>0.9,"podo_a_cor.csv")
      write.csv(cor(Y_list[[6]])>0.9,"podo_b_cor.csv")  
      write.csv(cor(Y_list[[7]])>0.9,"sclero_cor.csv")
      
      ## anni corr
      
      which(cor(Y_list[[1]])>0.9,arr.ind=T)
      
      #bio_13             18  17
      #bio_12             17  18
      
      #clim_stab_temp     30  23
      #bio_1              23  30
    
      #ORCDRC_M_sl2_250m  47  46
      #OCSTHA_M_sd1_250m  46  47
      
      Y_list[[1]] <-
        Y_list[[1]][,-which(colnames(Y_list[[1]]) %in% c("bio_13","clim_stab_temp","OCSTHA_M_sd1_250m"))]
      
      table(which(cor(Y_list[[1]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[1]])>0.9,arr.ind=T)[,2])
      
      which(cor(Y_list[[2]])>0.9,arr.ind=T)
  
      #ORCDRC_M_sl2_250m  40  39
      #OCSTHA_M_sd1_250m  39  40
      
      #coast_dist         32  31
      #refugia_dist       31  32
      
      #bio_13             14  13
      #bio_12             13  14
      
      Y_list[[2]] <-
        Y_list[[2]][,-which(colnames(Y_list[[2]]) %in% c("bio_13", "coast_dist","OCSTHA_M_sd1_250m"))]
      
      table(which(cor(Y_list[[2]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[2]])>0.9,arr.ind=T)[,2])
      
      #tail(cor(Y_list[[3]])>0.9)
      #bio_13             11  10
      #bio_12             10  11
      
      #bio_6              21  16
      #bio_1              16  21
      
      #coast_dist         31  30
      #refugia_dist       30  31
      
      #ORCDRC_M_sl2_250m  39  38
      #OCSTHA_M_sd1_250m  38  39
      
      which(cor(Y_list[[3]])>0.9,arr.ind=T)
      
      Y_list[[3]] <-
        Y_list[[3]][,-which(colnames(Y_list[[3]]) %in% c("bio_13", "bio_6", "coast_dist","OCSTHA_M_sd1_250m"))]
      
      table(which(cor(Y_list[[3]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[3]])>0.9,arr.ind=T)[,2])
      
      which(cor(Y_list[[4]])>0.9,arr.ind=T)
    
      #bio_2              19   1
      #MEM1                1  19
      
      #bio_13             13  12
      #bio_12             12  13
      
      #bio_2              19  32
      #coast_dist         32  19
      
      #coast_dist         32  31
      #refugia_dist       31  32
      
      #ORCDRC_M_sl2_250m  40  39
      #OCSTHA_M_sd1_250m  39  40
      
      #bio12+bio13
      Y_list[[4]] <-
        Y_list[[4]][,-which(colnames(Y_list[[4]]) %in% c("bio_2", "bio_13", "coast_dist","OCSTHA_M_sd1_250m"))]
      
      table(which(cor(Y_list[[4]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[4]])>0.9,arr.ind=T)[,2])
      
      which(cor(Y_list[[5]])>0.9,arr.ind=T)
    
      
      #bio_19              6   3
      #bio_14              3   6
      
      #bio_4               7   4
      #bio_15              4   7
      
      #clim_stab_temp     10   8
      #bio_5               8  10
      
      #alt                 9  14
      #Aquox              14   9
      
      #alt                 9  20
      #CLYPPT_M_sl2_250m  20   9
      
      Y_list[[5]] <-
        Y_list[[5]][,-which(
          colnames(Y_list[[5]]) %in% c(
            "bio_1",
            "bio_2",
            "bio_3",
            "bio_12",
            "bio_6",
            "refugia_dist",
            "SLTPPT_M_sl2_250m",
            "alt",
            "clim_stab_temp",
            "bio_15",
            "bio_19"
          )
        )]
      
      table(which(cor(Y_list[[5]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[5]])>0.9,arr.ind=T)[,2])
      
      which(cor(Y_list[[6]])>0.9,arr.ind=T)
      
      #bio_15              7   1
      #MEM1                1   7
      #
      #bio_4              13   1
      #bio_4              13   7
      #MEM1                1  13
      #bio_15              7  13
      #
      #bio_13              5   4
      #bio_12              4   5
      #
      #bio_1              10  14
      #bio_5              14  10
      #
      #bio_1              10  17
      #clim_stab_temp     17  10
      #
      #bio_5              14  17
      #clim_stab_temp     17  14
      
      #ORCDRC_M_sl2_250m  27  26
      #OCSTHA_M_sd1_250m  26  27
      
      Y_list[[6]] <-
        Y_list[[6]][,-which(colnames(Y_list[[6]]) %in% c("bio_4", "bio_15", "bio_13", "bio_5","clim_stab_temp","OCSTHA_M_sd1_250m"))]
      
      table(which(cor(Y_list[[6]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[6]])>0.9,arr.ind=T)[,2])
      
      which(cor(Y_list[[7]])>0.9,arr.ind=T)
      
      #bio_13             18  17
      #bio_12             17  18
      
      #bio_4              26  20
      #bio_15             20  26
      
      #bio_5              27  23
      #bio_1              23  27
      
      #CLYPPT_M_sl2_250m  42  29
      #alt                29  42
      
      #ORCDRC_M_sl2_250m  45  44
      #OCSTHA_M_sd1_250m  44  45
      
      Y_list[[7]] <-
        Y_list[[7]][,-which(colnames(Y_list[[7]]) %in% c("CLYPPT_M_sl2_250m", "bio_13", "bio_15", "bio_5","OCSTHA_M_sd1_250m"))]
      
      table(which(cor(Y_list[[7]])>0.9,arr.ind=T)[,1]-which(cor(Y_list[[7]])>0.9,arr.ind=T)[,2])
      
      # perform a db-RDA global model including all explanatory variables.
      rda1_list <- list()
      rda0_list <- list()
      rdaG_list <- list()
      sel_list <- list()
      anova_list <- list()
      library(car)
        for (i in 1:length(spe)) {
          rda1_list[[i]] = rda(X_list[[i]], Y_list[[i]])
          
          #Looking at the Variance Inflation Factor (VIF) for multicolinearity within the model.
          #When VIF is greater than around 10, this is problematic.
          
          ###
          #remove var with vif>10
          ###
          
          #loop through and gradually remove based on vif
          
          #while the biggest vif value is greater than 10
          while (max(na.omit(vif(rda1_list[[i]]))) > 10) {
            
            #sort based on vif (decreasing)
            large_vif <-
              names(sort(vif(rda1_list[[i]])[vif(rda1_list[[i]]) > 10], decreasing = T))
            
            #remove this individual from the list
            Y_list[[i]] <-
              Y_list[[i]][,-which(names(Y_list[[i]]) %in% large_vif[1])]
            
            #rerun rda
            rda1_list[[i]] = rda(X_list[[i]], Y_list[[i]])
            
            #print species and vif
            print(spe[i])
            print(vif(rda1_list[[i]]))
            
          }
          
          #final rda model
          print(spe[i])
          print(vif(rda1_list[[i]]))
          
          #Then we look at the explained variance by the global db-RDA model.
          #Adjusted R2 accounts for the number of variables and then is a measure of the unbiased amount of explained variation.
          RsquareAdj(rda1_list[[i]])
          #RsquareAdj(rda1) was equal 0.45 so that means this model explains 44% of the entire genomic variation.
          
          #Test for the significancy/probability of the model by using a ANOVA (Analysis of variance) test on the db-RDA global model.
          #Permutation tests are preferently used here as it does not required data normality distribution.
          #Permutation tests or ramdomization test compute the distribution when the null hypothesis
          #is true using a large number of times (nperm) and recomputing the statistic each time.
          #The test compare the true value of the statistic to the reference null distribution
          anova_list[[i]]<-anova(rda1_list[[i]], perm = 999)
          
          #Then, use the function OrdiR2Stepfor selecting the relevant variables to use in the db-RDA.
          #This function allows to add and remove variables in order to maximise the explained variance.
          #To avoid overfitting, selected variables should not explained more than the global model (e.g. 45%).
          rda0_list[[i]] <- rda(X_list[[i]] ~ 1, Y_list[[i]])
          
          #OrdiR2step will move towards the global model with all explanatory variables
          rdaG_list[[i]] <- rda(X_list[[i]] ~ ., Y_list[[i]])
          
          #Selection of variables until the rda global model is reached
          sel_list[[i]] <-
            ordiR2step(
              rda0_list[[i]],
              scope = formula(rdaG_list[[i]]),
              direction = "both",
              trace = F
            )
          
          #selection of variables
          print(spe[i])
          print(sel_list[[i]]$anova)
          
          #' Take a look at the points.
          #'
          #+ echo = FALSE
          
          #look at spatial scale of most important var
          ade4::s.value(Coorxy_list[[i]], dbmem_list[[i]][, 1])
          
          #+ echo = TRUE
          
          
        }
        
        
        ##
        # Build a model with the selected variables and visualize the results
        ##
        
        Ysel_list <- list()
        rdaS_list <- list()
        site_list <- list()
        form_split_list <- list()
        
        for (i in 1:length(spe)) {
          
          #extract variables from formula
          form <- as.character(sel_list[[i]]$call)[2]
          form <- gsub("X_list\\[\\[i\\]\\] \\~ ", "", form)
          form <- gsub(" \\+ ", ",", form)
          
          form_split <- strsplit(form, split = ",")[[1]]
          form_split_list[[i]] <- strsplit(form, split = ",")[[1]]
          
          Y_list[[i]][, form_split]
          
          Ysel_list[[i]] <- Y_list[[i]][, form_split]
          
          #Make sure that you have the right number of variables to include in your model.
          rdaS_list[[i]] <- rda(X_list[[i]] , Ysel_list[[i]])
          #summary(rdaS, scaling=1)
          
          RsquareAdj(rdaS_list[[i]])
          
          #Check the RDA summary. Scaling 1 allows the interpretation to focus on the ordination of objects because
          #the distances among objects approximate their Euclidean distances in the space of response variables.
          
          #Do a db-RDA biplot. First save the results of the db-RD in an object named site.
          site_list[[i]] = cbind(scores(
            rdaS_list[[i]],
            display = "sites",
            choices = c(1, 2),
            scaling = 1
          ),
          locs_list[[i]][, "sp"])
          
          #Assess how much percent of variation is explained by each axis in order to add this information to the biplot.
          summary(eigenvals(rdaS_list[[i]], model = "constrained"))
          
        }
        
        
        #Create the frame of the RDA. We selected a scaling equal to 2,
        #that means a correlation biplot (response variable focused):
        #distances between objects are not approximate Euclidean distances.
        #Angles between all vectors reflect linear correlation.
        
        #' Take a look at the points.
        #'
        #+ echo = FALSE
        
        library(RColorBrewer)
        cols <- brewer.pal(length(spe), "Set2")
        
        sca <- c(100, 100, 100, 175, 125, 150, 175)
        
        for (i in 1:length(spe)) {
          plot(
            rdaS_list[[i]],
            scaling = 2,
            main = spe[i],
            type = "none",
            xlab = c("db-RDA-1"),
            ylab = c("db-RDA-2"),
            xlim = c(-35, 35),
            ylim = c(-35, 35)
          )
          
          points(
            as.numeric(site_list[[i]][, 1]) * (sca[i] / 5),
            as.numeric(site_list[[i]][, 2]) * (sca[i] / 5),
            col = "black",
            bg = cols[i],
            pch = 21,
            cex = 1.2
          )
          
          #Add arrows showing selected variables pointing to their db-RDA scores.
          #Again you can rescale the biplot in order to minimise overlapping of arrows and points.
          
          arrows(
            0,
            0,
            scores(
              rdaS_list[[i]],
              display = "bp",
              choices = 1,
              scaling = 1
            ) * sca[i],
            scores(
              rdaS_list[[i]],
              display = "bp",
              choices = 2,
              scaling = 1
            ) * sca[i],
            col = "black",
            length = 0.1
          )
          
          text(
            scores(
              rdaS_list[[i]],
              display = "bp",
              choices = 1,
              scaling = 1
            ) * sca[i],
            scores(
              rdaS_list[[i]],
              display = "bp",
              choices = 2,
              scaling = 1
            ) * sca[i],
            labels = form_split_list[[i]],
            col = "black",
            cex = 0.8,
            pos = 3
          )
          
        }
        
        #+ echo = T
        
        ##save image
        save.image("dba.Rdata")

        