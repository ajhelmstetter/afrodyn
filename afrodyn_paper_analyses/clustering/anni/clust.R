rm(list=ls())
library(MASS)
library(pegas)
library(vcfR)
setwd("~/Dropbox/afrodyn_data/clustering/anni/")

##
# Sort IDs and names
##

id_ind<-read.csv("ID_index.csv")
id_vou<-read.csv("ID_voucher.csv")
only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
sb<-data.frame(as.character(foo$index),as.character(foo$voucher))
colnames(sb)<-c('index','voucher')
sb$index<-as.character(sb$index)
sb$voucher<-as.character(sb$voucher)
sb<-sb[match(sort(sb$index),sb$index),]
sb<-sb[order(sb$index),]

###
# read in VCF
###
vcf <- read.vcfR("anni_filtered_maf.vcf",checkFile = T, convertNA = T) #read in all data
head(vcf)
vcf

### convert to genlight
aa.genlight <- vcfR2genlight(vcf, n.cores=1)

### change index to voucher
indNames(aa.genlight) <-sb$voucher

### PCA
pca.1 <- glPca(aa.genlight, nf=300, n.cores=2) # retain first 300 axes (for later use in find.clusters); slow function
#save.image("~/Dropbox/afrodyn_data/stable_areas/data/anni.Rdata")

plot(pca.1$scores[,1],pca.1$scores[,2])
text(pca.1$scores[,1],pca.1$scores[,2], labels=names(pca.1$scores[,1]), cex= 0.7)

### find clusters
grp <- find.clusters(aa.genlight, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000)
300
4
dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1)
80
3

#find number of axes to use
opt_pc <- optim.a.score(dapc1)
opt_pc$best

#redo for optim a score
dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1)
1
1

#cluster scatterplot
scatter(dapc1)

scatter(dapc1, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

#cluster membership
table(dapc1$grp)

#transform vcf to genind for downstream
aa.genind <- vcfR2genind(vcf)

#set populations based on DAPC resuts
dapcgrp<-dapc1$grp
dapcgrp<-dapcgrp[order(match(names(dapcgrp),indNames(aa.genind)))]
pop(aa.genind) <- dapcgrp

###
#CROSS VALIDATION # DO NOT RUN UNLESS NEEDED
###

set.seed(999)
#pdf("cross_val.pdf")
#pramx <- xvalDapc(tab(aa.genind, NA.method = "mean"), pop(aa.genind))
#system.time(xval <- xvalDapc(tab(aa.genind, NA.method = "mean"), pop(aa.genind),n.pca.max=300,
#                             result = "groupMean", center = TRUE, scale = FALSE,
#                             n.pca = NULL, n.rep = 30, xval.plot = TRUE,
#                            parallel = "multicore", ncpus = 3))
#dev.off()

###
# TESS3
###

rm(list=ls())

require(tess3r)
require(maps)
require(LEA)
library(vcfR)
library(scales)
library(RColorBrewer)

#convert vcf to lfmm
#had to remove a SNP ln 5721 that couldn't be imported
vcf2lfmm("anni_filtered_tess.vcf", force = TRUE)

#format names and coords
vcf <- read.vcfR("anni_filtered_maf.vcf",checkFile = T, convertNA = T)
aa.genlight <- vcfR2genlight(vcf, n.cores=1)
id_ind<-read.csv("ID_index.csv")
id_vou<-read.csv("ID_voucher.csv")
only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
foo<-foo[order(match(foo$index,indNames(aa.genlight))),]
coordinates<-cbind(foo$ddlong,foo$ddlat)

#read in LFMM file
genotype<-read.lfmm("anni_filtered_tess.lfmm")
genotype[genotype == 9] <- NA
genotype[1,]


#plot map of coords
library(maps)
plot(coordinates, pch = 19, cex = .5,
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

#run TESS3
tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:10,
                   method = "projected.ls", ploidy = 2, openMP.core.num = 2)

#cross validation
#png("~/Dropbox/projects/AJH_AFRODYN/AFRODYN_shared/afrodyn_paper/review2_files/supp_figures/anni_xval.png",height=750,width=750)
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
#dev.off()

# retrieve tess3 Q matrix for K = 4 clusters
q.matrix <- qmatrix(tess3.obj, K = 4)

#set palette
my.palette <- CreatePalette(color.vector = c("#00bfc4",
                                             "#c77cff",
                                             "#7cae00",
                                             "#f8766d"
))

#make rownames readable
rownames(q.matrix)<-foo$voucher

# STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix, border = 1, space = 0,
        ylab = "Ancestry proportions", las = 3, cex.names = .3,
        col.palette=my.palette) -> bp


#TESS map 
library(raster)
library(rworldmap)

plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),
     main = NULL,
     xlab = "Longitude", ylab = "Latitude",
     resolution = c(300,300), cex = 1, cex.lab=1.5, cex.axis = 1.5,
     col.palette = my.palette)

save.image("clustering.Rdata")