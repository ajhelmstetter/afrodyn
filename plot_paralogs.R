
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
