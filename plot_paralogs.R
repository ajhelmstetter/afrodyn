
#set to directory you want PDF output
setwd("/home/helmstet/Dropbox/projects/AJH_AFRODYN/annonaceae_family/paralogs/")

library(ape)

#write pdf
pdf("./paralog_trees.pdf")

#set paths to the trees folder you downloaded
#should contain all trees output from paralog_tree.sh
files <- list.files(path="/home/helmstet/Dropbox/projects/AJH_AFRODYN/annonaceae_family/paralogs/trees/", pattern="*", full.names=T, recursive=FALSE)
filenames <- list.files(path="/home/helmstet/Dropbox/projects/AJH_AFRODYN/annonaceae_family/paralogs/trees/", pattern="*", full.names=F, recursive=FALSE)
par(mfrow=c(2,1))
par(mar=c(1,1,1,1))
for (i in 1:length(files))
{
  t<-read.tree(files[i])
  plot(t,type='un',cex=0.5)
  add.scale.bar()
  title(main=filenames[i],cex=0.5)
  plot(t,cex=0.5)
  add.scale.bar()
}
dev.off()

###
# Palms only
###

#reference for frequencies of exons per gene
ex<-read.table("exons.txt")
ref_tab<-table(ex)

#list of your gene names and frequencies
#remove exon indicator e.g. "_1" so only the gene name e.g. "1877" remains for each exon
scl<-read.table("inspected_exons_tab.txt")
scl_tab<-table(scl)

ref_inex_tab<-ref_tab[names(ref_tab)%in%names(scl_tab)]

#make sure names are in same order in ref_inex_tab and scl_tab
asd<-cbind(scl_tab,ref_inex_tab,scl_tab/ref_inex_tab)
asd<-data.frame(asd)

asd<-asd[asd$V3>=0.5,]

asd<-asd[asd$scl_tab>1,]

asd

data.frame(rownames(asd))
