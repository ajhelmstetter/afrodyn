setwd("~/Dropbox/afrodyn_data/clustering/sensitivity/")
library(scales)
library(RColorBrewer)
library(maps)
par(mfrow=c(1,1))

palette(alpha(brewer.pal(8,"Set2"),0.5))

list_rdat<-list.files("~/Dropbox/afrodyn_data/stable_areas/data",full.names = T)
list_locs<-list.files("~/Dropbox/afrodyn_data/stable_areas/locs",full.names = T)
list_locs

j<-2:8
i<-1

sp<-c("anni","anon","green","mona","podoa","podob","sclero")

clust_list<-list()

for(i in 1:length(list_rdat)){

  par(mfrow=c(2,2))
  load(list_rdat[i])

  p<-1

  for(k in 1:length(j)){

    ### find clusters
    grp <- find.clusters(aa.genlight, max.n.clust=20,
                         n.clust=j[k], n.pca=300, glPca = pca.1, n.iter=1e6, n.start=1000)

    dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1, n.pca=300, n.da = (j[k]-1))

    opt_pc <- optim.a.score(dapc1,plot=F)

    dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1, n.pca=opt_pc$best, n.da = (j[k]-1))

    locs<-read.csv(list_locs[i])
    locs<-locs[order(locs$index),]

    if(i == 5){
      locs$cluster<-dapc1$grp[1:(length(dapc1$grp)-3)]
    } else {
      locs$cluster<-dapc1$grp
    }



    par(mar=c(0.5,0.5,0.5,0.5))
    map('world', xlim=c(5,16), ylim=c(-5,6), col='gray90', fill=TRUE)
    legend("bottomleft",legend=paste(sp[i],"k =",j[k]))
    #map.axes() #Add axes
    palette(brewer.pal(8,"Set2"))
    points(locs$long,locs$lat,col=alpha(locs$cluster,0.5),pch=16,cex=2)

    if(k==1){
      loc_df<-data.frame(locs$cluster)
    } else {
      loc_df<-data.frame(loc_df,locs$cluster)
    }


  }

  colnames(loc_df)<-paste("k",seq(2,8,1),sep="")
  clust_list[[i]]<-loc_df
}

names(clust_list)<-c("anni","anon","green","mona","podoa","podob","sclero")

save.image("clust_sense.Rdata")

###
# different methods to infer significance
###

load("~/Dropbox/projects/AJH_AFRODYN/genetic_diversity/gen_div.Rdata")
library(distancetocoast)
library(MASS)
sp<-c("A. affinis","A. mannii","G. suaveolens","M. enghiana","P. acaulis","P. barteri","S. mannii")
p<-1
spe<-c("aa","am","gs","me","pa","pb","sm")

best_mod_list<-list()

for(i in 1:length(lf)){

  png(paste("~/Dropbox/projects/AJH_AFRODYN/figures/gen_div/",spe[i],".png",sep=""),height=2500,width=2000,res=300)


  par(mfrow=c(3,2))

  ref_he<-dat_list[[i]]

  for(k in 3:5){
    mod<-glm(ref_he[,k]~ref_he$dist+ref_he$coast_dist)
    mod2<-glmulti(mod)
    mod<-mod2@objects[[1]]
    summary(mod)
    ref_mods[[paste(spe[i],colnames(ref_he[,k]),sep="")]]<-mod

    plot(mod2@objects[[1]]$model[,1]~mod2@objects[[1]]$model[,2],
         xlab=colnames(ref_he)[k],
         ylab="Distance from nearest refugium",
         col=alpha(1,0.5),
         pch=16)

    legend("topright",
           legend=c(sp[i],
                    paste("p =",round(summary(mod)$coefficients[,4][2],7))),
           bty="n",
           cex=0.5
    )
    abline(mod,col=alpha(2,0.5),lty=2,lwd=3)
    mtext(letters[p],side=2,line=2,adj=0,las=1,padj=-7)
    p<-p+1
  }

  dev.off()

  #best_mod_list[[i]]<-mod
  p<-1

}



  par(mfrow=c(3,2))

  ref_he<-dat_list[[i]]

  for(k in 3:5){
    mod<-glm(ref_he$dist~ref_he[,k])
    summary(mod)
    ref_mods[[paste(spe[i],colnames(ref_he[,k]),sep="")]]<-mod

    plot(ref_he$dist~ref_he[,k],
         xlab=colnames(ref_he)[k],
         ylab="Distance from nearest refugium",
         col=alpha(1,0.5),
         pch=16)

    legend("topright",
           legend=c(sp[i],
                    paste("adj r2 =",round(summary(mod)$adj.r.squared,5)),
                    paste("p =",round(summary(mod)$coefficients[,4][2],7))),
                    bty="n",
                    cex=0.5
                    )
    abline(mod,col=alpha(2,0.5),lty=2,lwd=3)
    mtext(letters[p],side=2,line=2,adj=0,las=1,padj=-7)
    p<-p+1
  }


  for(k in 3:5){
    mod<-lm(ref_he$coast_dist~ref_he[,k])
    summary(mod)
    ref_mods[[paste(spe[i],colnames(ref_he[,k]),sep="")]]<-mod

    plot(ref_he$coast_dist~ref_he[,k],
         xlab=colnames(ref_he)[k],
         ylab="Distance from coast",
         col=alpha(1,0.5),
         pch=16)

    legend("topright",
           legend=c(spe[i],
                    paste("adj r2 =",round(summary(mod)$adj.r.squared,5)),
                    paste("p =",round(summary(mod)$coefficients[,4][2],7))),
           bty="n",
           cex=0.5
    )
    abline(mod,col=alpha(2,0.5),lty=2,lwd=3)
      mtext(letters[p],side=2,line=2,adj=0,las=1,padj=-9)
      p<-p+1
  }

  dev.off()

  p<-1

}
