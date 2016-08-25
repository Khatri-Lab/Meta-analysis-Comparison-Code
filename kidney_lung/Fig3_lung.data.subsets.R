
require(methods)
source("meta_compare.functions.R")


## smaller object, ease of reading:
lungResults <-readRDS("lungResults.RDS")
lung.results.gems <- lungResults$gems
rm(lungResults)

nstudies <- length(lung.results.gems)
combos <- combn(nstudies, 5)
combos <- as.list(data.frame(combos)) ## convert to list for mclapply

##################  limit gene set  ########################
## if modify, change name of SAVERDS file, too
overlap <- Reduce(intersect, lapply(lung.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
lung.results.gems <- lapply(lung.results.gems, function(gem) {
    gem$keys <- gem$keys[gem$keys %in% overlap]
    gem$expr <- gem$expr[names(gem$keys), ]
    gem})


out.list <- mclapply(combos, mc.cores=10, function(combo) {
    #Take each subset of 5 datasets, store class lengths
    lung.subset <- lung.results.gems[combo]
    classLengths <- unlist(lapply(lung.subset , function(gem) length(gem$class)))
    
    #Get initial effect sizes, then combine probes -> genes
    all.ES <- lapply(lung.subset, effect.sizes ) 
    all.effect.mini.lung <- study.effects(all.ES)
    
    mini.lung.DL <- combine.study.effects(all.effect.mini.lung, between.method="DL", everything=F, verbose=T)
    mini.lung.SJ <- combine.study.effects(all.effect.mini.lung, between.method="SJ", everything=F, verbose=T)
    mini.lung.HS <- combine.study.effects(all.effect.mini.lung, between.method="HS", everything=F, verbose=T)
    
    list(DL=mini.lung.DL, SJ=mini.lung.SJ, 
         HS=mini.lung.HS, sizes=classLengths)
})  
names(out.list) <- make.names(1:length(combos))

saveRDS(out.list, file=paste0("./data/mini.lung.reps.out.5subsets.OVERLAPGENES.RDS"))
