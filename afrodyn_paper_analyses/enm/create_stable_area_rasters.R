    rm(list = ls())
    setwd("~/enm_10001/")
    library(biomod2)
    library(ConR)
    
    par(mfrow=c(1,2))
    par(mar=c(2,2,2,2))
    
    locs_list = list.files(path="sp_locs",full.names = T)
    sp<-c("anni",
          "anon",
          "mona",
          "green",
          "podoa",
          "podob",
          "sclero")
    
    pdf("stableareas.pdf")
    
    for(i in 1:length(sp)){
    
      #read locations
      locs<-read.csv(locs_list[i])
      locs<-locs[,c(2:3)]
      head(locs)
      colnames(locs)[1:2]<-c("Longitude","Latitude")
    
      pres_kp<-raster(paste(sp[i],"/proj_ensemble_pres_bin/proj_ensemble_pres_bin_",sp[i],"_ensemble_KAPPAbin.grd",sep=""))
    
    
      ###
      # Convex hull
      ###
      conv.hull <- EOO.computing(locs[,2:1],export_shp=T)
    
      # 5 degrees around initial convex hull
      conv.hull_1 <- raster::buffer(x=conv.hull[[2]], width=10) ### Buffering the hull convex
      crs(conv.hull_1) <- crs(pres_kp)
    
      ##
      # read in ensemble projection rasters for past/present and KAPPA/TSS metric
      ##
    
      pres_kp<-raster(paste(sp[i],"/proj_ensemble_pres_bin/proj_ensemble_pres_bin_",sp[i],"_ensemble_KAPPAbin.grd",sep=""))
      pres_tss<-raster(paste(sp[i],"/proj_ensemble_pres_bin/proj_ensemble_pres_bin_",sp[i],"_ensemble_TSSbin.grd",sep=""))
  
      
      ccl_kp<-raster(paste(sp[i],"/proj_ensemble_ccl_bin/proj_ensemble_ccl_bin_",sp[i],"_ensemble_KAPPAbin.grd",sep=""))
    
      #reclassify values of presence (1) to 0.5
      ccl_kp<-reclassify(ccl_kp, cbind(1, 0.5))
    
      ccl_tss<-raster(paste(sp[i],"/proj_ensemble_ccl_bin/proj_ensemble_ccl_bin_",sp[i],"_ensemble_TSSbin.grd",sep=""))
      ccl_tss<-reclassify(ccl_tss, cbind(1, 0.5))
  
      mrl_kp<-raster(paste(sp[i],"/proj_ensemble_mrl_bin/proj_ensemble_mrl_bin_",sp[i],"_ensemble_KAPPAbin.grd",sep=""))
      mrl_kp<-reclassify(mrl_kp, cbind(1, 0.5))
    
      mrl_tss<-raster(paste(sp[i],"/proj_ensemble_mrl_bin/proj_ensemble_mrl_bin_",sp[i],"_ensemble_TSSbin.grd",sep=""))
      mrl_tss<-reclassify(mrl_tss, cbind(1, 0.5))
    
      mel_kp<-raster(paste(sp[i],"/proj_ensemble_mel_bin/proj_ensemble_mel_bin_",sp[i],"_ensemble_KAPPAbin.grd",sep=""))
      mel_kp<-reclassify(mel_kp, cbind(1, 0.5))
    
      mel_tss<-raster(paste(sp[i],"/proj_ensemble_mel_bin/proj_ensemble_mel_bin_",sp[i],"_ensemble_TSSbin.grd",sep=""))
      mel_tss<-reclassify(mel_tss, cbind(1, 0.5))
    
      #check plot
      plot(pres_kp)
      plot(pres_tss)
      plot(ccl_kp)
      plot(ccl_tss)
      plot(mrl_kp)
      plot(mrl_tss)
      plot(mel_kp)
      plot(mel_tss)
    
      ##
      # Plot rasters
      ##
    
      par(mfrow=c(1,2))
      par(mar=c(5,5,5,5))
    
      #add togather present and past binary rasters
    
      #this produces three categories:
      #0.5 = past presence only
      #1 = present-day presence
      #1.5 = present during both time periods
      plot(pres_kp+ccl_kp,main=paste(sp[i],"Stable areas CCL KAPPA"),legend=T)
      legend("bottomleft",legend = c("0.5 = past only","1  = present only","1.5 = present & past"),pch=16,cex=0.5,bty="n")
      points(locs,cex=0.5)
    
      plot(pres_tss+ccl_tss,main=paste(sp[i],"Stable areas CCL TSS"),legend=T)
      legend("bottomleft",legend = c("0.5 = past only","1  = present only","1.5 = present & past"),pch=16,cex=0.5,bty="n")
      points(locs,cex=0.5)
    
      plot(pres_kp+mel_kp,main=paste(sp[i],"Stable areas MEL KAPPA"),legend=T)
      legend("bottomleft",legend = c("0.5 = past only","1  = present only","1.5 = present & past"),pch=16,cex=0.5,bty="n")
      points(locs,cex=0.5)
    
      plot(pres_tss+mel_tss,main=paste(sp[i],"Stable areas MEL TSS"),legend=T)
      legend("bottomleft",legend = c("0.5 = past only","1  = present only","1.5 = present & past"),pch=16,cex=0.5,bty="n")
      points(locs,cex=0.5)
    
      plot(pres_kp+mrl_kp,main=paste(sp[i],"Stable areas MRL KAPPA"),legend=T)
      legend("bottomleft",legend = c("0.5 = past only","1  = present only","1.5 = present & past"),pch=16,cex=0.5,bty="n")
      points(locs,cex=0.5)
    
      plot(pres_tss+mrl_tss,main=paste(sp[i],"Stable areas MRL TSS"),legend=T)
      legend("bottomleft",legend = c("0.5 = past only","1  = present only","1.5 = present & past"),pch=16,cex=0.5,bty="n")
      points(locs,cex=0.5)
    
      ##
      # stable rasters
      ##
      
      tmp<-reclassify(pres_kp+ccl_kp,cbind(rbind(0.5,1,1.5),rbind(0,0,1)))
      
      plot(tmp)
      
      #reclassify to stability binary rasters
      #1.5 (suitable in both past and present) -> 1
      writeRaster(reclassify(pres_kp+ccl_kp,cbind(rbind(0.5,1,1.5),rbind(0,0,1))),file=paste(sp[i],"ccl_kp.tif",sep = "_"),format="GTiff", overwrite=T)
      writeRaster(reclassify(pres_tss+ccl_tss,cbind(rbind(0.5,1,1.5),rbind(0,0,1))),file=paste(sp[i],"ccl_tss.tif",sep = "_"),format="GTiff", overwrite=T)
      writeRaster(reclassify(pres_kp+mel_kp,cbind(rbind(0.5,1,1.5),rbind(0,0,1))),file=paste(sp[i],"mel_kp.tif",sep = "_"),format="GTiff", overwrite=T)
      writeRaster(reclassify(pres_tss+mel_tss,cbind(rbind(0.5,1,1.5),rbind(0,0,1))),file=paste(sp[i],"mel_tss.tif",sep = "_"),format="GTiff", overwrite=T)
      writeRaster(reclassify(pres_kp+mrl_kp,cbind(rbind(0.5,1,1.5),rbind(0,0,1))),file=paste(sp[i],"mrl_kp.tif",sep = "_"),format="GTiff", overwrite=T)
      writeRaster(reclassify(pres_tss+mrl_tss,cbind(rbind(0.5,1,1.5),rbind(0,0,1))),file=paste(sp[i],"mrl_tss.tif",sep = "_"),format="GTiff", overwrite=T)
    
    }
    
    #mv into folder
    system("mv *kp.tif ~/Dropbox/afrodyn_data/enm/stable_areas/")
    system("mv *tss.tif ~/Dropbox/afrodyn_data/enm/stable_areas/")
    
    dev.off()
    
