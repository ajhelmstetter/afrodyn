rm(list = ls())
library(MASS)
library(pegas)
library(vcfR)
library(poppr)
library(scales)
library(RColorBrewer)
palette(brewer.pal(8, "Set2"))

setwd("~/Dropbox/afrodyn_data/amova/north_south/")

#load in VCFs
load("~/Dropbox/afrodyn_data/amova/vcfs.Rdata")

vcfs <-
  list(anni_vcf,
       anon_vcf,
       green_vcf,
       fries_vcf,
       podo_a_vcf,
       podo_b_vcf,
       sclero_vcf)

names(vcfs) <-
  c("anni", "anon", "green", "mona", "podo_a", "podo_b", "sclero")

#location data
anni_locs<-read.csv("../../stable_areas/locs/anni.csv")
anon_locs<-read.csv("../../stable_areas/locs/anon.csv")
green_locs<-read.csv("../../stable_areas/locs/green.csv")
mona_locs<-read.csv("../../stable_areas/locs/mona.csv")
podo_a_locs<-read.csv("../../stable_areas/locs/podo_a.csv")
podo_b_locs<-read.csv("../../stable_areas/locs/podo_b.csv")
sclero_locs<-read.csv("../../stable_areas/locs/sclero.csv")

###
# anni
###

anni_locs$subpop <- anni_locs$cluster

#plot pops and N/S inversion
plot(
  anni_locs$long,
  anni_locs$lat,
  pch = 15,
  cex = 1.5,
  col = alpha(anni_locs$subpop, 0.5)
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:4)),
       pch = 15,
       col = 1:length(unique(anni_locs$subpop)))


#assign pops to north or south of climatic inversion discontinuity
anni_locs$north_south <- anni_locs$cluster
anni_locs$north_south[anni_locs$north_south == 1] <- "south"
anni_locs$north_south[anni_locs$north_south == 2] <- "north"
anni_locs$north_south[anni_locs$north_south == 3] <- "south"
anni_locs$north_south[anni_locs$north_south == 4] <- "south"

anni_locs$pop_subpop <-
  paste(anni_locs$subpop, anni_locs$north_south, sep = "_")

###
# anon
###

anon_locs$subpop <- anon_locs$cluster
plot(x = anon_locs$long,
     y = anon_locs$lat,
     col = "white")
text(
  x = anon_locs$long,
  y = anon_locs$lat,
  labels = anon_locs$cluster
)

#plot pops and N/S inversion
plot(
  anon_locs$long,
  anon_locs$lat,
  pch = 2,
  cex = 1.5,
  col = alpha(anon_locs$subpop,0.5)
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:8)),
       pch = 2,
       col = 1:length(unique(anon_locs$subpop)))

#assign pops to north or south of climatic inversion discontinuity
anon_locs$north_south <- anon_locs$cluster
anon_locs$north_south[anon_locs$north_south == 1] <- "north"
anon_locs$north_south[anon_locs$north_south == 2] <- "south"
anon_locs$north_south[anon_locs$north_south == 3] <- "north"
anon_locs$north_south[anon_locs$north_south == 4] <- "north"
anon_locs$north_south[anon_locs$north_south == 5] <- "south"
anon_locs$north_south[anon_locs$north_south == 6] <- "south"
anon_locs$north_south[anon_locs$north_south == 7] <- "south"
anon_locs$north_south[anon_locs$north_south == 8] <- "north"

anon_locs$pop_subpop <-
  paste(anon_locs$subpop, anon_locs$north_south, sep = "_")

###
# green
###

##
# Sort IDs and names
##
setwd("~/Dropbox/projects/AJH_AFRODYN/clustering/green/")
id_ind <- read.csv("ID_index.csv")
id_vou <- read.csv("ID_voucher.csv")
only_seq <- id_vou[id_vou$ID %in% id_ind$ID,]
sort(id_ind$ID)
id_ind[order(id_ind$ID),]
foo <-
  cbind(only_seq[order(only_seq$ID),], id_ind[order(id_ind$ID),])
sb <- data.frame(as.character(foo$index), as.character(foo$voucher))
colnames(sb) <- c('index', 'voucher')
sb$index <- as.character(sb$index)
sb$voucher <- as.character(sb$voucher)
sb <- sb[match(sort(sb$index), sb$index),]
sb <- sb[order(sb$index),]
sb <- sb[c(1:34, 54:length(sb$index)),]
sb <- sb[match(sb$voucher, green_locs$x.keyName),]

green_locs$index <- sb$index
green_locs$subpop <- green_locs$cluster

#plot pops and N/S inversion
plot(
  green_locs$coords.ddlong,
  green_locs$coords.ddlat,
  pch = 15,
  cex = 1.5,
  col = green_locs$subpop
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:6)),
       pch = 15,
       col = 1:length(unique(green_locs$subpop)))


#assign poops to north or south of equator
green_locs$north_south <- green_locs$cluster

#not sure whether south (3) or north (4) (of 1.5N)
green_locs$north_south[green_locs$north_south == 1] <- "north"

green_locs$north_south[green_locs$north_south == 2] <- "south"
green_locs$north_south[green_locs$north_south == 3] <- "south"
green_locs$north_south[green_locs$north_south == 4] <- "north"
green_locs$north_south[green_locs$north_south == 5] <- "north"

#not sure whether south (25) or north (32) (of 1.5N)
green_locs$north_south[green_locs$north_south == 6] <- "north"

green_locs$pop_subpop <-
  paste(green_locs$subpop, green_locs$north_south, sep = "_")

###
# mona
###
setwd("~/Dropbox/afrodyn_data/amova/north_south/")

mona_locs$subpop <- mona_locs$cluster

#plot pops and N/S inversion
plot(
  mona_locs$long,
  mona_locs$lat,
  pch = 15,
  cex = 1.5,
  col = mona_locs$subpop
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:3)),
       pch = 15,
       col = 1:length(unique(mona_locs$subpop)))


#assign pops to north or south of climatic inversion discontinuity
mona_locs$north_south <- mona_locs$cluster
mona_locs$north_south[mona_locs$north_south == 1] <- "north"

#some north
mona_locs$north_south[mona_locs$north_south == 2] <- "south"

mona_locs$north_south[mona_locs$north_south == 3] <- "north"
mona_locs$pop_subpop <-
  paste(mona_locs$subpop, mona_locs$north_south, sep = "_")

###
# podo_a
###

podo_a_locs$subpop <- podo_a_locs$cluster

plot(
  podo_a_locs$long,
  podo_a_locs$lat,
  pch = 15,
  cex = 1.5,
  col = podo_a_locs$subpop
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:2)),
       pch = 15,
       col = 1:length(unique(podo_a_locs$subpop)))

#assign pops to north or south of climatic inversion discontinuity
#all south
podo_a_locs$north_south <- podo_a_locs$cluster
podo_a_locs$north_south[podo_a_locs$north_south == 1] <- "south"
podo_a_locs$north_south[podo_a_locs$north_south == 2] <- "south"

podo_a_locs$pop_subpop <-
  paste(podo_a_locs$subpop, podo_a_locs$north_south, sep = "_")

podo_a_locs <-
  podo_a_locs[match(indNames(vcfR2genind(podo_a_vcf)), podo_a_locs$foo.index),]

###
# podo_b
###

podo_b_locs$subpop <- podo_b_locs$cluster

plot(
  podo_b_locs$long,
  podo_b_locs$lat,
  pch = 15,
  cex = 1.5,
  col = podo_b_locs$subpop
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:5)),
       pch = 15,
       col = 1:length(unique(podo_b_locs$subpop)))

#assign pops to north or south of climatic inversion discontinuity
podo_b_locs$north_south <- podo_b_locs$cluster
podo_b_locs$north_south[podo_b_locs$north_south == 1] <- "south"
podo_b_locs$north_south[podo_b_locs$north_south == 2] <- "north"
podo_b_locs$north_south[podo_b_locs$north_south == 3] <- "north"
podo_b_locs$north_south[podo_b_locs$north_south == 4] <- "south"
podo_b_locs$north_south[podo_b_locs$north_south == 5] <- "south"

podo_b_locs$pop_subpop <-
  paste(podo_b_locs$subpop, podo_b_locs$north_south, sep = "_")

###
# sclero
###

sclero_locs$subpop <- sclero_locs$cluster

plot(
  sclero_locs$long,
  sclero_locs$lat,
  pch = 15,
  cex = 1.5,
  col = sclero_locs$subpop
)
abline(h = 3, lty = 2)
abline(h = 0, lty = 2)
legend("topright",
       as.character(seq(1:2)),
       pch = 15,
       col = 1:length(unique(sclero_locs$subpop)))

#assign pops to north or south of climatic inversion discontinuity
#as only two pops cannot be included in AMOVA
sclero_locs$north_south <- sclero_locs$cluster
sclero_locs$north_south[sclero_locs$north_south == 1] <- "south"
sclero_locs$north_south[sclero_locs$north_south == 2] <- "north"

sclero_locs$pop_subpop <-
  paste(sclero_locs$subpop, sclero_locs$north_south, sep = "_")

locs <-
  list(anni_locs,
       anon_locs,
       green_locs,
       mona_locs,
       podo_a_locs,
       podo_b_locs,
       sclero_locs)

#check labels
setdiff(locs[[1]]$grps.bar.index, indNames(vcfR2genind(vcfs[[1]])))
setdiff(indNames(vcfR2genind(vcfs[[1]])),locs[[1]]$grps.bar.index)
data.frame(indNames(vcfR2genind(vcfs[[4]])), locs[[4]]$grps.bar.index)

####
# AMOVA
####

setwd("~/Dropbox/afrodyn_data/amova/north_south/")

#empty lists
list_amova <- list()
list_amova_signif <- list()
list_rand_amova <- list()
list_rand_amova_signif <- list()

sp <-
  c("anni", "anon", "green", "mona", "podo_a", "podo_b", "sclero")

#get rid of species that have only two clusters as N/S comparison cannot be made
vcfs <- vcfs[c(1, 2, 3, 4, 6)]
sp <- sp[c(1, 2, 3, 4, 6)]
locs <- locs[c(1, 2, 3, 4, 6)]

par(mfrow = c(1, 3))

for (i in 1:length(vcfs)) {
  
  #open pdf
  pdf(paste(names(vcfs)[i], "_north_south.pdf", sep = ""))
  
  ### convert to genind
  aa.genind <- vcfR2genind(vcfs[[i]], return.alleles = TRUE)
  
  #calculate genetic distance
  distgenEUCL <-
    dist(
      aa.genind,
      method = "euclidean",
      diag = FALSE,
      upper = FALSE,
      p = 2
    )
  
  #plot hist
  hist(distgenEUCL)

  #make strata data frame
  df3 <-
    data.frame(locs[[i]]$pop_subpop, locs[[i]]$north_south, locs[[i]]$subpop)
  colnames(df3) <- c("pop_subpop", "north_south", "subpop")
  strata(aa.genind) <- df3
  
  #check strata numbers
  table(strata(aa.genind, ~ subpop))
  table(strata(aa.genind, ~ north_south / subpop, combine = FALSE))
  
  #run amova
  res_amova <-
    poppr.amova(
      aa.genind,
      ~ north_south / subpop,
      threads = 1,
      quiet = TRUE,
      within = F
    )
  
  list_amova[[i]] <- res_amova
  
  #perform rand test with 999 repetitions
  #swap individual labels
  res_amova_signif   <- randtest(res_amova, nrepet = 999)
  
  #plot significance results
  plot(res_amova_signif)
  
  #add to list
  list_amova_signif[[i]] <- res_amova_signif
  
  ###
  #repeat with randomized population structure
  ###
  
  #this swaps cluster / N/S assignments
  
  aa.genind.new <- aa.genind
  
  #randomise strata
  strata(aa.genind.new) <-
    strata(aa.genind)[sample(nInd(aa.genind)), -1]
  aa.genind.amova  <-
    poppr.amova(
      aa.genind.new,
      ~ north_south / subpop,
      threads = 1,
      quiet = TRUE,
      within = F
    )
  
  list_rand_amova[[i]] <- aa.genind.amova
  
  aa.genind.amova.test <- randtest(aa.genind.amova, nrepet = 999)
  list_rand_amova_signif[[i]] <- aa.genind.amova.test
  
  plot(aa.genind.amova.test, main = "list_rand_amova_signif")
  
  dev.off()
  
}


#make table with results 
for (p in 1:length(list_amova)) {
  foo <-
    cbind(list_amova[[p]]$componentsofcovariance[1:3, ],
          rep(sp[p], 3),
          list_amova_signif[[p]]$pvalue)
  if (p == 1) {
    res_cov <- foo
  } else {
    res_cov <- rbind(res_cov, foo)
  }
}

res_cov
colnames(res_cov) <- c("sigma", "percentage", "sp", "p")

res_cov <- res_cov[, c(3, 1, 2, 4)]

write.csv(res_cov, "res_cov_north_south.csv")

write.csv(
  res_cov,
  "/home/helmstet/Dropbox/projects/AJH_AFRODYN/AFRODYN_shared/afrodyn_paper/supp_tables/amova_res.csv"
)

save.image("north_south.Rdata")

###############################
# GROUPED BARPLOT
###############################


ref_type <-
  res_cov[grep("Between north_south", rownames(res_cov)),]
# library
library(ggplot2)

# Grouped
ggplot(ref_type, aes(fill = sp, y = percentage, x = sp)) +
  geom_bar(position = "dodge", stat = "identity") + xlab("Model type")  + ylab("% variance explained by refugia")

