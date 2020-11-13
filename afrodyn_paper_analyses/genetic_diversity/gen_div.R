  rm(list = ls())
  setwd("~/Dropbox/afrodyn_data/genetic_diversity/")
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
  library(vegan)
  library(distancetocoast)
  
  ###
  # read in data
  ###
  
  #location data
  lf <- list.files("~/Dropbox/afrodyn_data/genetic_diversity/locs/")
  lff <-
    list.files("~/Dropbox/afrodyn_data/genetic_diversity/locs/", full.names = T)
  
  #genetic data
  rd <-
    list.files("~/Dropbox/afrodyn_data/genetic_diversity/data/", full.names = T)
  
  #species names
  spe <- c("anni", "anon", "green", "mona", "podo_a", "podo_b", "sclero")
  
  #list placeholders
  ref_mods <- list()
  coast_mods <- list()
  dat_list <- list()
  
  ###
  # create distance-from-refugia raster
  ###
  
  par(mar = c(3, 3, 3, 3))
  
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
  xy
  
  #world shape for mask
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  #calculate distance from points for each cell
  d2 <- distanceFromPoints(r, xy)
  
  #removed as introduces NAs
  #d2 <- mask(d2, world)
  par(mfrow=c(1,1))
  plot(d2)
  
  ##
  # gen div loop
  ##
  
  pdf("gen_div_plots.pdf")
  
  for (i in 1:length(lf)) {
      
      #format locs
      locs <- read.csv(lff[i])
      locs <- locs[, c("index", "long", "lat")]
      locs <- locs[order(locs$index), ]
      
      #extract distances
      dist <- extract(d2, locs[, c(2:3)])
      
      #read in table with index/diversity stat
      load(rd[i])
      
      #format data
      gen <- vcfR2genind(vcf)
      pop(gen) <- c(1:length(indNames(gen)))
      pop(gen)
      
      #calculate expected heterozygosity per individual
      Hs(gen)
      
      #calculate observed heterozygosity
      
      #make list of gen files per individual
      n.pop <- seppop(gen)
      
      #placeholder lists
      mean.hexp <- vector()
      mean.hobs <- vector()
      inds <- vector()
      tot_alleles <- vector()
      
      #
      for (j in 1:length(n.pop)) {
        inds[j] <- length(indNames(n.pop[[j]]))
        sum_pop <- summary(n.pop[[j]])
        
        #mean He (removed NA)
        mean.hexp[j] <- mean(sum_pop$Hobs[!(is.na(sum_pop$Hobs))])
        
        #mean Ho (removed NA)
        mean.hobs[j] <- mean(sum_pop$Hexp[!(is.na(sum_pop$Hexp))])
        
        #total number of alleles
        tot_alleles[j] <- sum_pop$pop.n.all
          
      }
      
      #calculate allelic richness
      
      ar <- allelic.richness(gen)

    ardf <- data.frame(ar)
    
    #remove NAs
    ardf <- na.omit(ardf)
    
    #remove min.all column
    ar_vect <- colMeans(ardf)[2:length(ardf[1, ])]
    
    
    #make df of all the info
    gendf <-
      data.frame(indNames(gen),
                 Hs(gen),
                 mean.hexp,
                 mean.hobs,
                 ar_vect,
                 tot_alleles)
    colnames(gendf) <-
      c("index",
        "Hs",
        "mean.hexp",
        "mean.hobs",
        "allelic_richness",
        "tot_alleles")
    
    #order by index
    gendf <- gendf[order(gendf$index), ]
    
    #make sure that gendf contains only inds with location data
    gendf <- gendf[gendf$index %in% locs$index, ]
    
    #cbind location and diversity data
    ref_he <- cbind(gendf, dist)
    colnames(ref_he) <-
      c("index",
        "Hs",
        "mean.hexp",
        "mean.hobs",
        "allelic_richness",
        "tot_alleles",
        "dist")
    
    #remove He outliers?
    #ref_he<-ref_he[!(ref_he$Hs>0.3),]
    
    #expected het vs distance from refugia central point
    
    par(mfrow = c(2, 2))
    
    ###
    # REFUGIA
    ###
    
    #loop through Hs, mean hexp, mean hobs and allelic richness
    #build linear model of distance against diversity metric
    
    for (k in 2:5) {
      mod <- lm(ref_he$dist ~ ref_he[, k])
      summary(mod)
      #ref_mods[[paste(spe[i], colnames(ref_he[, k]), sep = "")]] <-
        mod
      
        #plot distance against gen div
      plot(
        ref_he$dist ~ ref_he[, k],
        xlab = colnames(ref_he)[k],
        ylab = "Distance from nearest refugium",
        col = alpha(1, 0.5),
        pch = 16
      )
      
      #add legend with r2 and p value
      legend(
        "topright",
        legend = c(spe[i],
                   paste(
                     "adj r2 =", round(summary(mod)$adj.r.squared, 5)
                   ),
                   paste(
                     "p =", round(summary(mod)$coefficients[, 4][2], 7)
                   )),
        bty = "n",
        cex = 0.5
      )
      
      #add abline
      abline(mod,
             col = alpha(2, 0.5),
             lty = 2,
             lwd = 3)
    }
    
    ###
    # COAST
    ###
    
    #calculate distance to coast
    #All are RasterLayer layer with the distance (metres) to the Natural Earth coastline
    
    #No NAs here compared to distance from refugia
    coast_dist <-
      raster::extract(distance_to_coastline_lowres, locs[, c(2:3)])
    
    #in km
    coast_dist <- coast_dist / 1000
    
    for (k in 2:5) {
      mod <- lm(coast_dist ~ ref_he[, k])
      summary(mod)
      ref_mods[[paste(spe[i], colnames(ref_he[, k]), sep = "")]] <-
        mod
      
      plot(
        coast_dist ~ ref_he[, k],
        xlab = colnames(ref_he)[k],
        ylab = "Distance from coast",
        col = alpha(1, 0.5),
        pch = 16
      )
      
      legend(
        "topright",
        legend = c(spe[i],
                   paste(
                     "adj r2 =", round(summary(mod)$adj.r.squared, 5)
                   ),
                   paste(
                     "p =", round(summary(mod)$coefficients[, 4][2], 7)
                   )),
        bty = "n",
        cex = 0.5
      )
      abline(mod,
             col = alpha(2, 0.5),
             lty = 2,
             lwd = 3)
    }
    
    
    #add data to list for downstream plotting
    dat_list[[spe[i]]] <- cbind(ref_he, coast_dist)
    
    par(mfrow = c(1, 1))
    
  }
  
  
  dev.off()
  
  #save data
  save.image("gen_div.Rdata")