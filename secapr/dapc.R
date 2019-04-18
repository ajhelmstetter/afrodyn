rm(list=ls())
library(MASS)
setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")

id_ind<-read.csv("ID_index.csv")
id_vou<-read.csv("ID_voucher.csv")

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

library(pegas)
library(vcfR)
vcf <- read.vcfR("anni_filtered.sub2.vcf",checkFile = T, convertNA = T) #read in all data
head(vcf) 
vcf
#write.vcf(vcf,"anni_filtered.sub2_written.vcf.gz")
### convert to genlight
aa.genlight <- vcfR2genlight(vcf, n.cores=1)
indNames(aa.genlight) <-sb$voucher

ploidy(aa.genlight)
# look at the total data matrix (0,1,2; white = missing data)
glPlot (aa.genlight) # takes some time

# N missing SNPs per sample
x <- summary(t(as.matrix(aa.genlight)))
write.table(x[7,], file = "missing.persample.txt", sep = "\t") # NAs, if present, are in seventh row of summary 

###plot total AFS of the dataset
mySum <- glSum(aa.genlight, alleleAsUnit = TRUE)
barplot(table(mySum), col="blue", space=0, xlab="Allele counts",
        main="Distribution of ALT allele counts in total dataset") 

pca.1 <- glPca(aa.genlight, nf=300, n.cores=2) # retain first 300 axes (for later use in find.clusters); slow function

pdf("scatterplot.pdf")
plot(pca.1$scores[,1],pca.1$scores[,2])
text(pca.1$scores[,1],pca.1$scores[,2], labels=names(pca.1$scores[,1]), cex= 0.7)
dev.off()

grp <- find.clusters(aa.genlight, max.n.clust=20, glPca = pca.1, perc.pca =
                       100, n.iter=1e6, n.start=1000) 
300
4
dapc1<-dapc(aa.genlight, grp$grp, glPca = pca.1)
40
3
### Calculate Nei's distances between individuals/pops
#aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE) # Nei's 1972 distancebetween indivs
#stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.dst") # exportmatrix - for SplitsTree

pdf("DAPC_scatter.pdf")
scatter(dapc1) 
#screeplot(dapc1)
#biplot(dapc1)
dev.off() 
col<- funky(6) 

pdf("DAPC_scatter2.pdf")
scatter(dapc1, cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
dev.off() 

pdf("DAPC_all.pdf",width=20,height=5)
compoplot(dapc1, cex.names = 0.4, col=col,show.lab=T,legend=F)
dev.off() 

table(dapc1$grp)

###
#CROSS VALIDATION
###
#transform vcf to genind
aa.genind <- vcfR2genind(vcf)
#set populations based on DAPC results
dapcgrp<-dapc1$grp
dapcgrp<-dapcgrp[order(match(names(dapcgrp),indNames(aa.genind)))]
pop(aa.genind) <- dapcgrp
#gc <- as.genclone(aa.genind)

library(poppr)
set.seed(999)
pdf("cross_val1.pdf")
pramx <- xvalDapc(tab(aa.genind, NA.method = "mean"), pop(aa.genind))
system.time(xval <- xvalDapc(tab(aa.genind, NA.method = "mean"), pop(aa.genind),n.pca.max=300,
                             result = "groupMean", center = TRUE, scale = FALSE,
                             n.pca = NULL, n.rep = 30, xval.plot = TRUE,
                             parallel = "multicore", ncpus = 3))
dev.off()

xval[2:6]

###
###DAPC map
###
coords<-only_seq[only_seq$voucher%in%row.names(pca.1$scores),]
coords<-coords[match(row.names(pca.1$scores),as.character(coords$voucher)),]

# ugly quick hack to prepare the grouping of populations, not individuals
x <- data.frame(keyName=names(grp$grp), value=grp$grp, row.names=NULL) #

bar<-cbind(x[order(x$keyName),],sb[order(sb$voucher),])

dapc_grps<-data.frame(bar$index,bar$value)
dapc_grps<-dapc_grps[order(dapc_grps$bar.value),]

write.csv(dapc_grps,"dapc_grps.csv")

### select 5 random from each DAPC grp

rand5<-do.call( rbind, lapply( split(dapc_grps, dapc_grps$bar.value) ,
                               function(df) df[sample(nrow(df), 5) , ] )
)

prunetree<-setdiff(sb$index,rand5$bar.index)


write.csv(prunetree,"dapc_grps_not_rand5.csv")
write.csv(rand5,"dapc_grps_rand5.csv")

library(ggmap)
options(digits=9)
map <- get_map(location = c(lon = 10, lat = 0), zoom = 6)
mapPoints <- ggmap(map) + geom_point(aes(x = as.numeric(as.character(coords$ddlong)) , y = as.numeric(as.character(coords$ddlat)) ,
                                         colour = factor(x$value)),data = x,size=6)
mapPoints + theme(legend.position="none") + scale_color_manual(values=hue_pal()(8)[1:4])


###
###read in RAxML tree
###
library(phylotools)
library(scales)
setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")

phy<-read.tree("~/Dropbox/projects/AJH_AFRODYN/annickia/raxml_anni_out_1563514/RAxML_bipartitions.anni")
phy$tip.label<-gsub("-","",phy$tip.label)
phy<-drop.tip(phy,setdiff(phy$tip.label,sb$index))
phy<-drop.tip(phy,setdiff(sb$index,phy$tip.label))

phy<-sub.taxa.label(phy,sb)
#phy<-root(phy,"A_affinis_Dimonika_1")
plot(phy,cex=0.5)

palette(hue_pal()(4))
phy<-ladderize(phy)
phygrp<-dapc1$grp
phygrp<-phygrp[order(match(names(phygrp),phy$tip.label))]
pdf("dapc_phy.pdf")
plot(phy,show.tip.label = F)
p <- character(length(phy$node.label))

phy$node.label[1]<-0
phy$node.label<-as.numeric(phy$node.label)
p[phy$node.label >= 85] <- "black"
p[phy$node.label < 85 & phy$node.label >= 70] <- "gray"
p[phy$node.label < 70] <- "white"
nodelabels(pch=21, cex=1, bg=p)
tiplabels(pch=21, cex=2, bg=phygrp,offset=0.0001)
dev.off()

phygrp

phy_coord<-foo[order(match(foo$voucher,names(phygrp))),]
X<-phy_coord[,c(3:4)]
rownames(X)<-phy_coord[,2]
library(phytools)

lat<-phy_coord$ddlat
long<-phy_coord$ddlong
names(lat)<-rownames(X)
names(long)<-rownames(X)
names(lat)<-gsub(" ","_",names(lat))
names(long)<-gsub(" ","_",names(long))


plot(phy)
nodelabels(cex=0.5)

anni<-extract.clade(phy,node=184)

anni$tip.label<-gsub(" ","_",anni$tip.label)

lata<-lat[names(lat)%in%anni$tip.label]

longa<-long[names(long)%in%anni$tip.label]

phylo.to.map(anni,cbind(lat,long),xlim=c(5,15),ylim=c(-5,2),cex=0.5,rotate = T,fsize=0.00001)

crds<-cbind(lata,longa)

anni2$edge.length<-anni2$edge.length*1000

write.csv(lata,"lat2.csv")
write.csv(longa,"long2.csv")
#raxml concat tree for all samples (all loci)
#phy<-read.tree('~/Dropbox/projects/AJH_AFRODYN/annickia/RAxML_bipartitions.annickia')

#Loops through and drops all a1 alleles
#tips_to_drop<-vector()
#for (i in 1:length(phy$tip.label)) {
#  if(grepl("a1", phy$tip.label[i]) == T) {
#    tips_to_drop[i]<-phy$tip.label[i]
#  } else {
#    
#  }
#}
#phy<-drop.tip(phy,na.omit(tips_to_drop))

###
###read in faststructure results
###
library(pophelper)
setwd("~/programs/fastStructure/anni/")
par(mfrow=c(1,1))
#read in data
dat4<-read.table("anni_filtered_final.4.meanQ")
colnames(dat4) <- c("Cluster1","Cluster2","Cluster3","Cluster4")
#form into list and check
q1 <- list("4"=dat4)
str(q1)
is.qlist(q1)
#rename rows
rownames(q1[[1]]) <- sb$voucher

#tree is missing  A_mannii_Cristal_8
#remove from DAPC/fS
q1$`4`<-q1$`4`[order(match(rownames(q1$`4`),phy$tip.label)),]

plotQ(q1[1],sortind="all",showindlab = T, useindlab = T,dpi = 1080,indlabsize = 0.75)

x<-x[order(match(x$keyName,phy$tip.label)),]

###
###plot tree with results from DAPC and fastSTRUCTURE
###
setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")
pdf("phy_dapc.pdf")
plot(phy,align.tip.label=T,cex=0.25)
tiplabels(tip=1:length(phy$tip.label),
          pie=cbind(q1[[1]]$Cluster1,q1[[1]]$Cluster2,q1[[1]]$Cluster3,q1[[1]]$Cluster4),
          piecol=c("blue","lightblue","red","green"),cex=0.5
          #,offset=0.0007
)
tiplabels(pch=21, cex=0.75, bg=x$value)

dev.off()

###
###ladderized tree with results from DAPC and fastSTRUCTURE
###
phy<-ladderize(phy)
setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")
pdf("phy_dapc_ladd.pdf")
plot(phy,align.tip.label=T,cex=0.25)
tiplabels(tip=1:length(phy$tip.label),
          pie=cbind(q1[[1]]$Cluster1,q1[[1]]$Cluster2,q1[[1]]$Cluster3,q1[[1]]$Cluster4),
          piecol=c("blue","lightblue","red","green"),cex=0.5
          #,offset=0.0007
)
tiplabels(pch=21, cex=0.8, bg=x$value)
x$col<-factor(x$value, levels = c(1:length(unique(x$value))), 
              labels = hue_pal()(length(unique(x$value))))
tiplabels(pch=21, cex=0.75, bg=as.character(x$col))
dev.off()


###
###ggtree with results from DAPC and fastSTRUCTURE
###
library(scales)

#source("https://bioconductor.org/biocLite.R")
#biocLite("ggtree")

x <- data.frame(keyName=names(grp$grp), value=grp$grp, row.names=NULL) #
x<-x[order(match(x$keyName,phy$tip.label)),]
x<-cbind(x,((colnames(q1$`4`)[apply(q1$`4`,1,which.max)])))
colnames(x)<-c("name","dapc","fastStructure")
x$fastStructure<-gsub("Cluster","",x$fastStructure)
x$fastStructure<-as.factor(x$fastStructure)
x2<-data.frame(x$dapc,x$fastStructure)
row.names(x2)<-x$name
colnames(x2)<-c("DAPC","fastSTRUCTURE")

library(ggtree)
library(ggstance)
phy$tip.label<-gsub("A_affinis_","",phy$tip.label)
phy$tip.label<-gsub("A_pilosa_","",phy$tip.label)

row.names(x2)<-gsub("A_affinis_","",row.names(x2))
row.names(x2)<-gsub("A_pilosa_","",row.names(x2))

p <- ggtree(phy) + geom_tiplab(size=2, align=TRUE, linesize=.5,offset=0.001) + theme_tree2() +xlab("Substitutions per site")
pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>% gheatmap(x2, offset=0.005, width=0.5, colnames=F) %>% scale_x_ggtree()
pp + theme(legend.position="none")

admix<-q1$`4`
#format fastSTRUCTURE dataframe for panel
c1<-data.frame(rownames(admix),admix$Cluster1,rep("Cluster1",length(admix$Cluster1)))
colnames(c1)<-c("id","value","category")
c2<-data.frame(rownames(admix),admix$Cluster2,rep("Cluster2",length(admix$Cluster1)))
colnames(c2)<-c("id","value","category")
c3<-data.frame(rownames(admix),admix$Cluster3,rep("Cluster3",length(admix$Cluster1)))
colnames(c3)<-c("id","value","category")
c4<-data.frame(rownames(admix),admix$Cluster4,rep("Cluster4",length(admix$Cluster1)))
colnames(c4)<-c("id","value","category")

clust<-rbind(c1,c2,c3,c4)

clust$id<-gsub("A_affinis_","",clust$id)
clust$id<-gsub("A_pilosa_","",clust$id)


xdapc<-data.frame(rownames(x2))
xdapc$category<-x2$DAPC
xdapc$value<-rep(1, length(xdapc$rownames.x2.))

colnames(xdapc)<-c("id","category","value")

xdapc$category<-as.character(xdapc$category)

foo<-xdapc

for (i in 1:length(unique(xdapc$category))){ 
a1<-xdapc[xdapc$category!=i,]
a1$value<-rep(0,length(a1$id))
a1$category<-rep(i,length(a1$id))
foo<-rbind(foo,a1)
}

library(ggstance)
##Maybe used for structure plot?
p <- ggtree(phy) + geom_tiplab(size=2, align=TRUE, linesize=.5,offset=0.001)
p2 <- facet_plot(p + xlim_tree(0.045), panel = 'DAPC', data = foo, 
                 geom = geom_barh, 
                 mapping = aes(x = value, fill = as.factor(category)), 
                 stat='identity')
p3 <- facet_plot(p2, panel = 'fastSTRUCTURE', data = clust, 
                 geom = geom_barh, 
                 mapping = aes(x = value, fill = as.factor(category)), 
                 stat='identity')
p3



###
###Faststructure individual map
###

coords<-only_seq[only_seq$voucher%in%row.names(q1$`4`),]
coords<-coords[match(row.names(q1$`4`),as.character(coords$voucher)),]

library(maps)
library(plotrix)
#plot map of world in relevant region

###
###Faststructure map with populations
###
map('world', xlim=c(5,16), ylim=c(-5,6), col='gray90', fill=TRUE)
map('world', xlim=c(5,16), ylim=c(-5,6), col='gray90', fill=TRUE)
map.axes() #Add axes
only_seq<-id_vou[id_vou$ID%in%id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo<-cbind(only_seq[order(only_seq$ID),],id_ind[order(id_ind$ID),])
sb<-data.frame(as.character(foo$index),as.character(foo$voucher),as.character(foo$pop))
colnames(sb)<-c('index','voucher','pop')
sb$index<-as.character(sb$index)
sb$voucher<-as.character(sb$voucher)
sb$pop<-as.character(sb$pop)

coords<-coords[order(match(coords$voucher,rownames(admix))),]
sb<-sb[order(match(sb$voucher,rownames(admix))),]

#make data frame with all info
allDat<-cbind(coords,admix,sb)

#aggregate data frame, calculate average coord/cluster membership
aggPop<-aggregate(allDat,by=list(allDat$pop),FUN = 'mean')
#calculate number of individuals per site and add to data frame
aPN <-cbind(aggPop,table(allDat$pop))

g_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#plot pie charts showing sampling locations scaled by number of individuals per location
for (i in 1:nrow(aPN)){floating.pie(aPN$ddlong[i],aPN$ddlat[i],c(aPN$Cluster1[i],aPN$Cluster2[i],aPN$Cluster3[i],aPN$Cluster4[i]),radius=0.5/aPN$Freq[i],col=hue_pal()(8)[5:8]) }


###
# DESCRIPTIVE STATISTICS
###
setwd("~/Dropbox/projects/AJH_AFRODYN/annickia/")
library("hierfstat")
library("pegas")

#filterd and unlinked
vcf <- read.vcfR("anni_filtered.sub2.vcf") #read in all data

#transform vcf to genind
aa.genind <- vcfR2genind(vcf)

#transform vcf to genlight
aa.genlight <- vcfR2genlight(vcf, n.cores=1)

#set individual names
indNames(aa.genind) <-sb$voucher
indNames(aa.genlight) <-sb$voucher

#set populations based on DAPC results
dapcgrp<-dapc1$grp
names(dapcgrp)
dapcgrp<-dapcgrp[order(match(names(dapcgrp),indNames(aa.genind)))]
pop(aa.genind) <- dapcgrp
pop(aa.genlight) <- dapcgrp

#pop cannot be named with numbers so rename with characters
#change for number of pops
popNames(aa.genind)<-c("a","b","c","d")
popNames(aa.genind)

#calculate allelic richness per population
ar<-allelic.richness(aa.genind[,-1])
ardf<-data.frame(ar)
ardf<-na.omit(ardf)
ar_vect<-colMeans(ardf)[2:8]


#Add geographic information to individuals in genind object
coords<-only_seq[only_seq$voucher%in%indNames(aa.genind),]
coords<-coords[order(match(coords$voucher,indNames(aa.genind))),]
xy<-data.frame(coords$ddlong,coords$ddlat)
other(aa.genind)<-xy
other(aa.genind)


#summary of genind object
summary(aa.genind)
div<-summary(aa.genind)

par(mar=c(5,5,5,5))

#Number of alleles vs number of individuals in population
plot(div$n.by.pop, div$pop.n.all, xlab="Colonies sample size",
     ylab="Number of alleles",main="Alleles numbers and sample sizes",
     type="n")
text(div$n.by.pop,div$pop.n.all,lab=names(div$n.by.pop))

#Expected - Observed heterozygosity
# >0 = inbreeding
barplot(sort(div$Hexp-div$Hobs), main="Heterozygosity: expected-observed",ylab="Hexp - Hobs",xlab="Locus",axisnames = F)

#Obs heterozygosity
plot(sort(div$Hobs), xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

#Expected heterozygosity as a function of observed heterozygosity per locus
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

#Are variances of Hexp and Hobs the same?
bartlett.test(list(div$Hexp, div$Hobs)) # a test : H0: Hexp = Hobs

#Hardy-Weinberg
#hwresults<-hw.test(aa.genind, B = 1000)
head(hwresults)

#Transform to object to be used with hierfstat package
Mydata2 <- genind2hierfstat(aa.genind)

#Summary stats for dataset
basic.stats(Mydata2)

#Calculate fstats for dataset
#fstat(aa.genind)

#pairwise fstatistics among pops
#matFst <- pairwise.fst(aa.genind)
#matFst


library(StAMPP) 
# Nei's 1972 distance among individuals
aa.D.ind <- stamppNeisD(aa.genlight, FALSE) 
#stamppPhylip(aa.D.ind, file="aa.indiv_Neis_distance.phy.dst") # export for SplitsTree
# Nei's 1972 distance between pops
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE) 
#stamppPhylip(aa.D.pop, file="aa.pops_Neis_distance.phy.dst") # export matrix - for SplitsTree 

#Plot heatmap with dendrogram to pdf
library(gplots)
colnames(aa.D.ind) <- rownames(aa.D.ind)
pdf(file="Neis_dist_heatmap.pdf", width=10, height=10)
heatmap.2(aa.D.ind, trace="none", cexRow=0.4, cexCol=0.4)
dev.off() 


####
# isolation by distance
####

data(nancycats)
#keep genind for individual level analyses
toto <- aa.genind
#use genpop for population level analyses
#toto <- genind2genpop(aa.genind)

#calculate distance matrices
Dgen <- dist(toto$tab)
Dgeo <- dist(data.frame(toto$other))
ibd <- mantel.randtest(Dgen,Dgeo)
ibd

#detect IBD? line shows real value relative to sims
plot(ibd)

#Plot distances, do they cluster into more than one group? (fragmented pattern)
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red",lty=2)

#plot heatmap of distances to pdf
pdf(file="ibd_heatmap.pdf", width=10, height=10)

###some problem with loading MASS package here, so load at top of script
#library(MASS)

dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")
dev.off()

####
#populations treated separately
####

aa.genind.sep <- seppop(aa.genind, drop=TRUE) #
#library(MASS)

ibd_obs<-vector()
ibd_p<-vector()

pdf(file="ibd_heatmap_pops.pdf", width=10, height=10)

par(mfrow=c(2,1))
for (i in 1:length(aa.genind.sep)){
  Dgen <- dist(aa.genind.sep[i][[1]]@tab)
  Dgeo <- dist(data.frame(aa.genind.sep[i][[1]]@other))
  ibd <- mantel.randtest(Dgen,Dgeo,nrepet = 1000)
  
  ibd_obs[i]<-ibd$obs
  ibd_p[i]<-ibd$pvalue
  
  plot(ibd)
  
  dens <- kde2d(Dgeo,Dgen, n=300)
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(Dgeo, Dgen, pch=20,cex=.5)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  abline(lm(Dgen~Dgeo))
  title("Isolation by distance plot")
  
}

dev.off()

cbind(ibd_obs,ibd_p)

#calculate distance matrices
Dgen <- dist(aa.genind.sep$a@tab)
Dgeo <- dist(data.frame(aa.genind.sep$a@other))
ibd <- mantel.randtest(Dgen,Dgeo)
ibd

#detect IBD? line shows real value relative to sims
plot(ibd)

#Plot distances, do they cluster into more than one group? (fragmented pattern)
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red",lty=2)

pdf(file="ibd_heatmap.pdf", width=10, height=10)
#library(MASS)
dens <- kde2d(Dgeo,Dgen, n=300)
myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(Dgeo, Dgen, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
abline(lm(Dgen~Dgeo))
title("Isolation by distance plot")
dev.off()

###
#monmonier
###

D <- dist(toto$tab)

xy<-data.frame(toto$other)

#some individuals from exactly the same coords so have to jitter so chooseCN works
xy$coords.ddlat<-jitter(xy$coords.ddlat)
xy$coords.ddlong<-jitter(xy$coords.ddlong)

#not sure how to set distances
gab <- chooseCN(data.frame(xy),ask=FALSE,type=2)
mon1 <- monmonier(xy,D,gab)
d

#pairwise.fst(toto)

plot(mon1)

plot(mon1,add.arrows=FALSE,bwd=10,col="black")
points(xy, cex=2, pch=20,
       col=fac2col(pop(toto), col.pal=spectral))
legend("topright",leg=c("Pop A", "Pop B"), pch=c(20),
       col=spectral(2), pt.cex=2)


  ###
  #Allele frequency spectra
  ###
  
  ###plot total AFS of the dataset
  mySum <- glSum(aa.genlight, alleleAsUnit = TRUE)
  barplot(table(mySum), col="blue", space=0, xlab="Allele counts",
          main="Distribution of ALT allele counts in total dataset")
  ###plot AFS per one pop
  aa.genlight.sep <- seppop(aa.genlight, drop=TRUE) #
  aa.genlight.sep$`1`
  # after seppop you must remove the nonvariant positions within thepopulation
  n.alleles.1 <-colSums(as.matrix(aa.genlight.sep$`1`)) # how manyalternative alleles are in each locus?
  summary(as.factor(n.alleles.1)) # how many particular categories of alternative allele counts are in my pop?
  aa.genlight.1 <- new("genlight", (as.matrix(aa.genlight.sep$`1`))
                       [,(colSums(as.matrix(aa.genlight.sep$`1`)) > 0) &
                           (colSums(is.na(as.matrix(aa.genlight.sep$`1`))) ==
                              0)]) # remove the reference-only positions AND remove columns with NA
  aa.genlight.1
  summary(colSums(as.matrix(aa.genlight.1))) # check if there are no zeros
  # plot unfolded AFS - for one pop.
  mySum <- glSum(aa.genlight.1, alleleAsUnit = TRUE)
  barplot(table(mySum), col="blue", space=0, xlab="Allele counts",
          main="Distribution of ALT allele counts in BEL") # plot the original counts of each category 
  
  #### plot AFS for all pops in a batch
  aa.genlight.sep <- seppop(aa.genlight, drop=TRUE) # separate genlight per population
  # remove the nonvariant positions AND columns with NA within that pop.
  aa.genlight.sep.2 <- lapply (aa.genlight.sep, function (pop)
  {new("genlight", (as.matrix(pop))[,(colSums(as.matrix(pop)) > 0)
                                    & (colSums(is.na(as.matrix(pop))) == 0)])})
  ##add pop identity to list elements
  listnames<-names(aa.genlight.sep.2)
  for (i in seq(listnames)) {pop(aa.genlight.sep.2[[i]])<-
    substr(indNames(aa.genlight.sep.2[[i]]),1,3)}
  # loop over each population in a list of populations and draw AFS into one fig
  pdf("AFS_all_barplot.pdf", width=5, height=5)
  par(mfrow=c(2,3),mar=c(2,2,2,0))
  mySum <- lapply (aa.genlight.sep.2, function (pop) {
    barplot(table(glSum(pop, alleleAsUnit=T)), col="blue", space=0,
            xlab="Allele counts",
            main=paste(levels(pop(pop)),sum(table(glSum(pop, alleleAsUnit=T))),"SNPs",
                       sep=" "))
  })
  dev.off()
  par(mfrow=c(1,1)) 
  
  
  #calculated mean hobs per pop
  n.pop <- seppop(aa.genind) 
  mean.hobs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 
  mean.hobs[is.nan(mean.hobs)] <- NA 
  barplot(mean.hobs) 
  
  mean.hexp <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hexp))) 
  mean.hexp[is.nan(mean.hexp)] <- NA 
  barplot(mean.hexp) 
  
  inds<- do.call("c", lapply(n.pop, function(x) summary(x)$n))
  
  tot_alleles<- do.call("c", lapply(n.pop, function(x) summary(x)$pop.n.all))
  
  ###
  #private alleles
  ###
  
  library(poppr)
  library("ggplot2")
  Pinfpriv <- private_alleles(toto, report = "data.frame")
  priv_mean<-aggregate(Pinfpriv$count,by=list(Pinfpriv$population),FUN="mean")
  priv_sum<-aggregate(Pinfpriv$count,by=list(Pinfpriv$population),FUN="sum")
  #plot doesnt correspond to data? eg population names?
  ggplot(Pinfpriv ) + geom_tile(aes(x = population, y = allele, fill = count)) + theme(axis.title.y=element_blank(),
                                                                                       axis.text.y=element_blank(),
                                                                                       axis.ticks.y=element_blank())
  
  
  sum_stats<-cbind(unique(x$dapc),inds,tot_alleles,priv_sum$x,ar_vect,mean.hobs,mean.hexp,mean.hexp-mean.hobs,ibd_obs,ibd_p)
  colnames(sum_stats)<-c("grp","n","total alleles","private alleles","allelic richness","hobs","hexp","hexp-hobs","IBD obs","IBD p")
  sum_stats_rnd<-round(sum_stats, digits = 3)
  write.csv(sum_stats_rnd,'sum_stats.csv')
  
  sum_stats_rnd<-round(sum_stats, digits = 3)
  sum_stats_rnd<-sum_stats_rnd[c(1:length(inds)),]
  write.csv(sum_stats_rnd,'sum_stats.csv')
  
  map('world', xlim=c(2,18), ylim=c(-8,8), col='gray90', fill=TRUE)
  map('world', xlim=c(2,18), ylim=c(-8,8), col='gray90', fill=TRUE)
  map.axes() #Add axes
  palette(hue_pal()(11))
  points(as.numeric(as.character(coords$ddlong)),as.numeric(as.character(coords$ddlat)),col=factor(x$dapc),pch=16,cex=2)
  
  plot(sum_stats_rnd)
  
  library(gridExtra)
  
  #Plot your table with table Grob in the library(gridExtra)
  sum_stats_rnd<-data.frame(sum_stats_rnd)
  sum_stats_rnd<-sum_stats_rnd[order(sum_stats_rnd$grp),]
  rownames(sum_stats_rnd)<-c()
  ss <- tableGrob(sum_stats_rnd)
  
  world <- map_data("world2")
  eb2<-ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group),col='black',fill='gray')
  eb3<-eb2 + xlim(4, 18) + ylim(-7, 7) + geom_point(aes(x = as.numeric(as.character(coords$ddlong)) , y = as.numeric(as.character(coords$ddlat)) ,
                                                        colour = x$dapc),data = x,size=6)
  
  ggsave("dapc_map.pdf",eb3)

  ggsave("dapc_stats.pdf", ss)
  
  
  
