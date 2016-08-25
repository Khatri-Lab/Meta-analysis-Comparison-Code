
require(methods)
source("meta_compare.functions.R")

## smaller object, ease of reading:
cardioResults <- readRDS("cardioResults.RDS")
cardio.results.gems <- cardioResults$gems
rm(cardioResults)

nstudies <- length(cardio.results.gems)
combos <- combn(nstudies, 5)
combos <- as.list(data.frame(combos)) ## convert to list for mclapply

##################  limit gene set  ########################
## if modify, change name of SAVERDS file, too
overlap <- Reduce(intersect, lapply(cardio.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
cardio.results.gems <- lapply(cardio.results.gems, function(gem) {
    gem$keys <- gem$keys[gem$keys %in% overlap]
    gem$expr <- gem$expr[names(gem$keys), ]
    gem})

library(parallel)
maxCores <- detectCores()
out.list <- mclapply(combos, mc.cores=maxCores, function(combo) {
    #Take each subset of 5 datasets, store class lengths
    cardio.subset <- cardio.results.gems[combo]
    classLengths <- unlist(lapply(cardio.subset , function(gem) length(gem$class)))
    
    #Get initial effect sizes, then combine probes -> genes
    all.ES <- lapply(cardio.subset, effect.sizes ) 
    all.effect.mini.cardio <- study.effects(all.ES)
    
    mini.cardio.DL <- combine.study.effects(all.effect.mini.cardio, between.method="DL", everything=F, verbose=T)
    mini.cardio.SJ <- combine.study.effects(all.effect.mini.cardio, between.method="SJ", everything=F, verbose=T)
    mini.cardio.HS <- combine.study.effects(all.effect.mini.cardio, between.method="HS", everything=F, verbose=T)
    
    list(DL=mini.cardio.DL, SJ=mini.cardio.SJ, 
         HS=mini.cardio.HS, sizes=classLengths)
})  
names(out.list) <- make.names(1:length(combos))

saveRDS(out.list, file=paste0("./data/mini.cardio.reps.out.5subsets.OVERLAPGENES.RDS"))
