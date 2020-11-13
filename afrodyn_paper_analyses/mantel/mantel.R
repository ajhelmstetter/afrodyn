rm(list=ls())
setwd("~/Dropbox/afrodyn_data/mantel/")
library(distancetocoast)
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

par(mar=c(3,3,3,3))

##
# distance from refugia raster
##

#create empty raster
r <- raster(ncol=5000,nrow=5000)
crs(r) = '+proj=utm +zone=12 +datum=WGS84'
e<-extent(0,30,-10,10)
r<-crop(r,e)

#approximate midpoints of refugia
xy <- c(10,5)
xy<-rbind(xy,c(11,3))
xy<-rbind(xy,c(10.5,1))
xy<-rbind(xy,c(10.5,-2))
xy<-rbind(xy,c(11.5,-1.5))
xy<-rbind(xy,c(12.5,-4.5))
xy

#world shape for mask
world <- ne_countries(scale = "medium", returnclass = "sf")

#calculate distance from points for each cell
d2 <- distanceFromPoints(r, xy)
#d2<-mask(d2,world)

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

## check closest refugia  
foo<-read.csv(lff[5])
plot(d2,xlim=c(10,15),ylim=c(-5,5))
points(foo$long,foo$lat)

###
# Geographic locations
### 

#species names
spe <- c("anni", "anon", "green", "mona", "podo_a", "podo_b", "sclero")


library(vegan)

#empty lists for storage
gen_list<-list()
distgenEUCL_list<-list()
locs_list<-list()
clim_list<-list()
ref_mant_list<-list()
ibd_list<-list()
dist_ref_list<-list()

pdf("mantel.pdf")

for(i in 1:length(lff)){
  
  par(mfrow=c(1,2))
  
  #create location data df and list
  locs<-read.csv(lff[i])
  locs<-locs[,c("index","long","lat")]
  locs<-locs[order(locs$index),]
  locs$sp<-rep(spe[i],length(locs$long))
  if(i == 1){
    all_locs<-locs
    locs_list[[i]]<-locs
  } else {
    all_locs<-rbind(all_locs,locs)
    locs_list[[i]]<-locs
  }
  
  #read in table with index/diversity stat
  
  if(i == 5){
    
    load(rd[i])
    
    #make genind
    gen<-vcfR2genind(vcf)
    
    #podo a has 3 inds with no coords
    gen<-gen[1:(length(indNames(gen))-3)]
    
    #1 pop per individual
    pop(gen)<-indNames(gen)
    
    #create genpop object
    gp<-genind2genpop(gen)

    #store in list
    gen_list[[i]]<-gp
    
  } else {
    
    load(rd[i])
    
    #make genind
    gen<-vcfR2genind(vcf)
    
    #1 pop per individual
    pop(gen)<-indNames(gen)
    
    #create genpop object
    gp<-genind2genpop(gen)
    
    #store in list
    gen_list[[i]]<-gp
    
  }
  
  #calculate genetic distance (euclidean)
  distgenEUCL_list[[i]] <- dist(gen_list[[i]],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
  
  par(mar=c(4,4,4,4))
  
  #mantel test
  ibd_list[[i]]<-mantel(distgenEUCL_list[[i]], dist(locs_list[[i]][,c(2:3)]), method="pearson", permutations=9999)
  
  #plot histogram
  hist(ibd_list[[i]]$perm, 
       main=paste(spe[i],"genetic distance vs geographic distance"),
       xlim=c((min(ibd_list[[i]]$statistic,ibd_list[[i]]$perm)),(max(ibd_list[[i]]$statistic,ibd_list[[i]]$perm)+0.1)),
       col="lightgrey",
       xlab="mantel statistic",ylab="frequency",
       breaks=30,
       border=F,
       xaxt = "n",
       cex.main=0.5,
       cex.axis=0.5)
  
  axis(1, at = seq(round((min(ibd_list[[i]]$statistic,ibd_list[[i]]$perm))-0.1,1),(max(ibd_list[[i]]$statistic,ibd_list[[i]]$perm)+0.1),0.1)) 
  
  #line showing empirical value
  abline(v=ibd_list[[i]]$statistic,lty=2,col=alpha(2,0.5),lwd=2)
  
  #plot pairwise genetic distance vs pairwise spatial distance
  plot(distgenEUCL_list[[i]],dist(locs_list[[i]][,c(2:3)]),xlab="genetic distance",ylab="spatial distance",
       main=paste(spe[i],"genetic distance vs geographic distance"),cex.main=0.5)

}

dev.off()

save.image("mantel.Rdata")

###
# avg dist between individuals
###

par(mfrow=c(3,3))
avg<-vector()
stdev<-vector()
for(i in 1:length(locs_list)){
  hist(dist(locs_list[[i]][,c(2:3)]),main=spe[i],xlab="pairwise geographic distance")
  avg[i]<-mean(dist(locs_list[[i]][,c(2:3)]))
  stdev[i]<-sd(dist(locs_list[[i]][,c(2:3)]))
}

data.frame(avg,stdev,row.names = spe)



