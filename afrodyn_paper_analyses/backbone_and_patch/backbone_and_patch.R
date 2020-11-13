rm(list=ls())
setwd("~/Dropbox/afrodyn_data/")
library(ape)
par(mfrow=c(1,1))

back<-read.nexus("~/Dropbox/afrodyn_data/backbone_annonaceae/back_ann_aug2020.tree")

library(phylotools)
asd<-read.csv("~/Dropbox/projects/AJH_AFRODYN/backbone_tree_annonaceae/sub_labels.csv")
back<-sub.taxa.label(back,asd)

plot(back)
bt_back<-branching.times(back)
nodelabels()
nodelabels(round(bt_back,3),cex=0.7)

#anni divergence time from outgroup
anni_out<-bt_back[names(bt_back) %in% "29"]
anni_out

anni<-read.nexus("starbeast/anni/combined_anni_star.tree")
plot(anni)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_anni<-branching.times(anni)
anni_root<-bt_anni[names(bt_anni) %in% "6"]
anni_root

anni_crown<-bt_anni[names(bt_anni) %in% "7"]
anni_crown

anni_rel<-(anni_crown/anni_root)
anni_rel

anni<-drop.tip(anni,"o")

abs_scale_stem<-anni_out*anni_rel

#make root height 1
anni$edge.length<-anni$edge.length/anni_crown

#check node heights
plot(anni)
nodelabels(branching.times(anni))

#rescale patch clade node heights
anni$edge.length<-anni$edge.length*abs_scale_stem

#rename tips
anni$tip.label<-c("A. affinis_a","A. affinis_b","A. affinis_c","A. affinis_d")

#check branching times
plot(anni)
nodelabels(round(branching.times(anni),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"Annickia_affinis")

#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem

#attach patch clade
back<-bind.tree(back,anni,where=1)

#check attachment
plot(back)

###
# Anon
###

#anon divergence time from outgroup
bt_back<-branching.times(back)
anon_out<-bt_back[names(bt_back) %in% "47"]
anon_out

anon<-read.nexus("starbeast/anon/combined_anon_star.tree")
plot(anon)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_anon<-branching.times(anon)
anon_root<-bt_anon[names(bt_anon) %in% "10"]
anon_root

anon_crown<-bt_anon[names(bt_anon) %in% "11"]
anon_crown

anon_rel<-(anon_crown/anon_root)
anon_rel

anon<-drop.tip(anon,"o")

abs_scale_stem<-anon_out*anon_rel

#make root height 1
anon$edge.length<-anon$edge.length/anon_crown

#check node heights
plot(anon)
nodelabels(round(branching.times(anon),3))

#rescale patch clade node heights
anon$edge.length<-anon$edge.length*abs_scale_stem

#rename tips
anon$tip.label<-c("A. mannii_a","A. mannii_b","A. mannii_c","A. mannii_d",
                  "A. mannii_e","A. mannii_f","A. mannii_g","A. mannii_h")

#check branching times
plot(anon)
nodelabels(round(branching.times(anon),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"Anonidium_mannii")

#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem

#attach patch clade
back<-bind.tree(back,anon,where=1)

#check attachment
plot(back)

###
# mona
###

fries<-read.nexus("starbeast/mona/combined_mona_star.tree")

#fries divergence time from outgroup
bt_back<-branching.times(back)
fries_out<-bt_back[names(bt_back) %in% "64"]
fries_out
plot(fries)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_fries<-branching.times(fries)
fries_root<-bt_fries[names(bt_fries) %in% "6"]
fries_root

fries_crown<-bt_fries[names(bt_fries) %in% "7"]
fries_crown

fries_rel<-(fries_crown/fries_root)
fries_rel

fries<-drop.tip(fries,c("o","ob"))

abs_scale_stem<-fries_out*fries_rel

#make root height 1
fries$edge.length<-fries$edge.length/fries_crown

#check node heights
plot(fries)
nodelabels(round(branching.times(fries),3))

#rescale patch clade node heights
fries$edge.length<-fries$edge.length*abs_scale_stem

#rename tips
fries$tip.label<-c("M. enghiana_a","M. enghiana_b","M. enghiana_c")

#check branching times
plot(fries)
nodelabels(round(branching.times(fries),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"Friesodielsia_enghiana")

#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem

#attach patch clade
back<-bind.tree(back,fries,where=6)

#check attachment
plot(back)


###
# Green
###

#fries divergence time from outgroup
bt_back<-branching.times(back)
plot(back)
nodelabels()
green_out<-bt_back[names(bt_back) %in% "48"]
green_out

green<-read.nexus("starbeast/green/combined_green_star.tree")
plot(green)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_green<-branching.times(green)
green_root<-bt_green[names(bt_green) %in% "8"]
green_root

green_crown<-bt_green[names(bt_green) %in% "9"]
green_crown

green_rel<-(green_crown/green_root)
green_rel

green<-drop.tip(green,c("o"))

abs_scale_stem<-green_out*green_rel

#make root height 1
green$edge.length<-green$edge.length/green_crown

#check node heights
plot(green)
nodelabels(round(branching.times(green),3))

#rescale patch clade node heights
green$edge.length<-green$edge.length*abs_scale_stem

#rename tips
green$tip.label<-c("G. suaveolens_a","G. suaveolens_b","G. suaveolens_c",
                   "G. suaveolens_d","G. suaveolens_e","G. suaveolens_f")

#check branching times
plot(green)
nodelabels(round(branching.times(green),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"Greenwayodendron_suaveolens")

#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem

#attach patch clade
back<-bind.tree(back,green,where=6)

#check attachment
plot(back)

#plot tree
par(mar=c(3,1,1,1))
back<-ladderize(back)
plot(back,cex=0.6,label.offset = 0.5,x.lim=c(0,130))
axisPhylo()


library(phyloch)
library(strap)
library(OutbreakTools)

back$root.time <- max(branching.times(back))

geoscalePhylo(tree=ladderize(back,right=FALSE), units=c("Epoch"), boxes="Epoch",cex.ts=0.75,cex.age=1,cex.tip=0.75,lwd=3, width=2,x.lim = c(-45,110))

plot(back)
back_anon<-back
write.tree(back_anon,"backbone_annonaceae/back_anon.tree")

####
####
# Palms
####
####

###
# P barteri
### 

back<-read.nexus("backbone_palms/combine_back_palms_mats_aug2020.tree")

library(phylotools)
asd<-read.csv("~/Dropbox/projects/AJH_AFRODYN/backbone_tree_palms/sub_labels.csv")
back<-sub.taxa.label(back,asd)

back<-drop.tip(back,c("SM094","S_mannii_mayoko_3","P_barteri_abeilles_4","Pa-DOUDOU6","PB-Campo 3"))

plot(back)
bt_back<-branching.times(back)
nodelabels()
nodelabels(round(bt_back,3),cex=0.7)
#anni divergence time from outgroup
anni_out<-bt_back[names(bt_back) %in% "25"]
anni_out


anni<-read.nexus("starbeast/podo/combined_podo_star.tree")
anni<-drop.tip(anni,"g")
plot(anni)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_anni<-branching.times(anni)
anni_root<-bt_anni[names(bt_anni) %in% "7"]
anni_root

anni_crown<-bt_anni[names(bt_anni) %in% "8"]
anni_crown

anni_rel<-(anni_crown/anni_root)
anni_rel

anni<-drop.tip(anni,"f")

abs_scale_stem<-anni_out*anni_rel

#make root height 1
anni$edge.length<-anni$edge.length/anni_crown

#check node heights
plot(anni)
nodelabels(branching.times(anni))

#rescale patch clade node heights
anni$edge.length<-anni$edge.length*abs_scale_stem

#rename tips
anni$tip.label<-c("P. barteri_a","P. barteri_b","P. barteri_c","P. barteri_d","P. barteri_e")

#check branching times
plot(anni)
nodelabels(round(branching.times(anni),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"P_barteri_Bamba_5")
back_edge
#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem
plot(back)
#attach patch clade
back<-bind.tree(back,anni,where=1)

#check attachment
plot(back)


###
# P acaulis
###

plot(back)
bt_back<-branching.times(back)
nodelabels()
nodelabels(round(bt_back,3),cex=0.7)
#anni divergence time from outgroup
anni_out<-bt_back[names(bt_back) %in% "29"]
anni_out


anni<-read.nexus("starbeast/podo/combined_podo_star.tree")
anni<-drop.tip(anni,c("b","c","d","e"))
plot(anni)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_anni<-branching.times(anni)
anni_root<-bt_anni[names(bt_anni) %in% "4"]
anni_root

anni_crown<-bt_anni[names(bt_anni) %in% "5"]
anni_crown

anni_rel<-(anni_crown/anni_root)
anni_rel

anni<-drop.tip(anni,"a")

abs_scale_stem<-anni_out*anni_rel

#make root height 1
anni$edge.length<-anni$edge.length/anni_crown

#check node heights
plot(anni)
nodelabels(round(branching.times(anni),4))

#rescale patch clade node heights
anni$edge.length<-anni$edge.length*abs_scale_stem

#rename tips
anni$tip.label<-c("P. acaulis_a","P. acaulis_b")

#check branching times
plot(anni)
nodelabels(round(branching.times(anni),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"Pa-MAYO1")
back_edge

#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem
plot(back)
#attach patch clade
back<-bind.tree(back,anni,where=1)

#check attachment
plot(back)

###
# Sclero
###

plot(back)
bt_back<-branching.times(back)
nodelabels()
nodelabels(round(bt_back,3),cex=0.7)
#anni divergence time from outgroup
anni_out<-bt_back[names(bt_back) %in% "39"]
anni_out


anni<-read.nexus("starbeast/sclero/combined_sclero_star.tree")
plot(anni)
nodelabels()

#get the relative branching times
#and divide the age of the ingroup crown by that of the root, obtaining the relative
#scale of the stem edge from outgroup to ingroup; prune outgroups

bt_anni<-branching.times(anni)
anni_root<-bt_anni[names(bt_anni) %in% "4"]
anni_root

anni_crown<-bt_anni[names(bt_anni) %in% "5"]
anni_crown

anni_rel<-(anni_crown/anni_root)
anni_rel

anni<-drop.tip(anni,"o")

abs_scale_stem<-anni_out*anni_rel

#make root height 1
anni$edge.length<-anni$edge.length/anni_crown

#check node heights
plot(anni)
nodelabels(round(branching.times(anni),4))

#rescale patch clade node heights
anni$edge.length<-anni$edge.length*abs_scale_stem

#rename tips
anni$tip.label<-c("S. mannii_a","S. mannii_b")

#check branching times
plot(anni)
nodelabels(round(branching.times(anni),3))

#get edge where patch will be attached
back_edge<-which.edge(back,"S_mannii_Chaillu_7")
back_edge

#reduce stem length to fit in patch clade
back$edge.length[back_edge]<-back$edge.length[back_edge]-abs_scale_stem
plot(back)
#attach patch clade
back<-bind.tree(back,anni,where=5)

#check attachment
plot(back)

back$root.time <- max(branching.times(back))

back_palms<-back
write.tree(back_palms,"~/Dropbox/afrodyn_data/backbone_palms/back_palms.tree")


library(strap)
geoscalePhylo(tree=ladderize(back,right = T), units=c("Epoch"), boxes="Epoch",   cex.tip=0.75, cex.age=0.75, cex.ts=0.75, label.offset=0.5,lwd=3)
geoscalePhylo(tree=ladderize(back,right=FALSE), units=c("Epoch"), boxes="Epoch",cex.ts=0.75,cex.age=1,cex.tip=0.75,lwd=3, width=2,x.lim = c(-45,110))

save.image("~/Dropbox/projects/AJH_AFRODYN/backbone_and_patch/backbone_and_patch.Rdata")

