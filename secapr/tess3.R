rm(list=ls())
require(tess3r)
require(maps)
require(LEA)

par(mfrow=c(1,1))

setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")

#had to gunzip sub2.vcf before import
#had to remove a SNP that couldn't be imported
vcf2lfmm("anni_filtered.sub2_written.vcf", force = TRUE)



vcf <- read.vcfR("anni_filtered.sub2.vcf",checkFile = T, convertNA = T) 
aa.genlight <- vcfR2genlight(vcf, n.cores=1)


id_ind<-read.csv("ID_index.csv")
id_vou<-read.csv("ID_voucher.csv")
only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]

foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
foo<-foo[order(match(foo$index,indNames(aa.genlight))),]
coordinates<-cbind(foo$ddlong,foo$ddlat)

genotype<-read.lfmm("anni_filtered.sub2_written.lfmm")
genotype[genotype == 9] <- NA

library(maps)
coordinates[1:3,]

plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:20, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 2) 

plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# retrieve tess3 Q matrix for K = 5 clusters 
q.matrix <- qmatrix(tess3.obj, K = 4)

hue_pal()(4)

hue_pal()(4)

c("#7CAE00","#C77CFF","#F8766D","#00BFC4")

my.palette <- CreatePalette(color.vector = c("#F8766D","#C77CFF","#7CAE00","#00BFC4"))


rownames(q.matrix)<-sb$voucher

library(scales)
# STRUCTURE-like barplot for the Q-matrix 

pdf(width = 15,height=5)
barplot(q.matrix, border = NA, space = 0, 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        main = "Ancestry matrix",col.palette=my.palette) -> bp
## Use CreatePalette() to define color palettes.

sb<-data.frame(as.character(foo$index),as.character(foo$voucher))
colnames(sb)<-c('index','voucher')
sb$index<-as.character(sb$index)
sb$voucher<-as.character(sb$voucher)
sb<-sb[match(sort(sb$index),sb$index),]
sb<-sb[order(sb$index),]

sb$voucher[bp$order]

axis(1, at = 1:nrow(q.matrix), labels = sb$voucher[bp$order], las = 3, cex.axis = .7) 
dev.off()


q.matrix

rownames(q.matrix)<-c(indNames(aa.genlight))
library(raster)
library(rworldmap)
plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude",
     resolution = c(300,300), cex = .4,
     col.palette = my.palette)

plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", xlim=c(8.5,13.5),ylim=c(-5,5),
     resolution = c(300,300), cex = .4,
     col.palette = my.palette)

