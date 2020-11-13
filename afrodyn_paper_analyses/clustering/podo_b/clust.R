rm(list=ls())
library(MASS)
library(pegas)
library(vcfR)
setwd("~/Dropbox/afrodyn_data/clustering/podo_b/")

##
# Sort IDs and names
##

library(ape)
library(MASS)
id_ind<-read.csv("ID_index.csv")
id_vou<-read.csv("ID_voucher.csv")


#use raxml phylo to get index names of p. barteri individuals
phy<-read.tree('~/Dropbox/afrodyn_data/raxml/Podococcus_b/RAxML_bipartitions.podo_b_out')
phy<-ladderize(phy)
plot(phy)
nodelabels()
phy<-root(phy,node = 160)
plot(phy)
nodelabels()
edgelabels()
podo_b<-extract.clade(phy,node=160,root.edge = 1)
plot(podo_b)
btip<-data.frame(podo_b$tip.label)

id_ind<-id_ind[id_ind$index%in%btip$podo_b.tip.label,]

only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
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
vcf <- read.vcfR("podo_b_filtered_maf.vcf",checkFile = T, convertNA = T) #read in all data
vcf 

### convert to genlight
aa.genlight <- vcfR2genlight(vcf, n.cores=1)
#indNames(aa.genlight) <- sb$voucher

### PCA
pca.1 <- glPca(aa.genlight, nf=300, n.cores=2) # retain first 300 axes (for later use in find.clusters); slow function
#save.image("~/Dropbox/afrodyn_data/stable_areas/data/podo_b.Rdata")

#pdf("scatterplot.pdf")
plot(pca.1$scores[,1],pca.1$scores[,2])
text(pca.1$scores[,1],pca.1$scores[,2], labels=names(pca.1$scores[,1]), cex= 0.7)
#dev.off()

### find clusters
grp <- find.clusters(aa.genlight, max.n.clust=20, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) 
300
5
dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1)
40
4

opt_pc <- optim.a.score(dapc1)
opt_pc

dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1)
5
4

#cluster scatterplot
#pdf("DAPC_scatter.pdf")
scatter(dapc1) 
#dev.off() 

#pdf("DAPC_scatter2.pdf")
scatter(dapc1, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
#dev.off() 

table(dapc1$grp)

###
#CROSS VALIDATION # DO NOT RUN UNLESS NEEDED
###

#transform vcf to genind
aa.genind <- vcfR2genind(vcf)

#set populations based on DAPC resuts
dapcgrp<-dapc1$grp
dapcgrp<-dapcgrp[order(match(names(dapcgrp),indNames(aa.genind)))]
pop(aa.genind) <- dapcgrp

library(poppr)
set.seed(999)
#pdf("cross_val.pdf")
#pramx <- xvalDapc(tab(aa.genind, NA.method = "mean"), pop(aa.genind))
#system.time(xval <- xvalDapc(tab(aa.genind, NA.method = "mean"), pop(aa.genind),n.pca.max=300,
#                             result = "groupMean", center = TRUE, scale = FALSE,
#                             n.pca = NULL, n.rep = 30, xval.plot = TRUE,
#                             parallel = "multicore", ncpus = 3))
#dev.off()

###
# TESS3
###

library(vcfR)
require(tess3r)
require(maps)
require(LEA)
library(scales)
library(RColorBrewer)

#had to remove a ln 4072 that couldn't be imported
#vcf2lfmm("podo_b_filtered_maf.vcf", force = TRUE)

vcf <- read.vcfR("podo_b_filtered_maf.vcf",checkFile = T, convertNA = T) 
aa.genlight <- vcfR2genlight(vcf, n.cores=1)

id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
foo<-foo[order(match(foo$index,indNames(aa.genlight))),]
coordinates<-cbind(as.numeric(as.character(foo$ddlong)),as.numeric(as.character(foo$ddlat)))

genotype<-read.lfmm("podo_b_filtered_maf.lfmm")
genotype[genotype == 9] <- NA


library(maps)
coordinates[1:3,]

plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:10, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 2) 

#png("~/Dropbox/projects/AJH_AFRODYN/AFRODYN_shared/afrodyn_paper/review2_files/supp_figures/podo_b_xval.png",height=750,width=750)
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")
#dev.off()

# retrieve tess3 Q matrix for K = 5 clusters 
q.matrix <- qmatrix(tess3.obj, K = 5)
my.palette <- CreatePalette(color.vector = c(brewer.pal(5,"Set2")))

# STRUCTURE-like barplot for the Q-matrix 
#pdf("tess_barplot.pdf",height=4,width=10)
barplot(q.matrix, border = NA, space = 0, 
        ylab = "Ancestry proportions", 
       col.palette=my.palette) -> bp
## Use CreatePalette() to define color palettes.
axis(1, at = 1:nrow(q.matrix), labels = foo$voucher, las = 3, cex.axis = .4) 
#dev.off()

library(raster)
library(rworldmap)
#png("tess.png",height=750,width=750)
plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  
     main = NULL,
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = 1, cex.lab=1.5, cex.axis = 1.5,
     col.palette = my.palette)
#dev.off()

save.image("clustering.Rdata")
