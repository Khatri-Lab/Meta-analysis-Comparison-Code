
require(methods)
source("meta_compare.functions.R")


## smaller object, ease of reading:
kidneyResults <-readRDS("kidneyResults.RDS")
kidney.results.gems <- kidneyResults$gems
rm(kidneyResults)

nstudies <- length(kidney.results.gems)
combos <- combn(nstudies, 5)
combos <- as.list(data.frame(combos)) ## convert to list for mclapply

##################  limit gene set  ########################
## if modify, change name of SAVERDS file, too
overlap <- Reduce(intersect, lapply(kidney.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
kidney.results.gems <- lapply(kidney.results.gems, function(gem) {
    gem$keys <- gem$keys[gem$keys %in% overlap]
    gem$expr <- gem$expr[names(gem$keys), ]
    gem})

library(parallel)
maxCores <- detectCores()
out.list <- mclapply(combos, mc.cores=maxCores, function(combo) {
    #Take each subset of 5 datasets, store class lengths
    kidney.subset <- kidney.results.gems[combo]
    classLengths <- unlist(lapply(kidney.subset , function(gem) length(gem$class)))
    
    #Get initial effect sizes, then combine probes -> genes
    all.ES <- lapply(kidney.subset, effect.sizes ) 
    all.effect.mini.kidney <- study.effects(all.ES)
    
    mini.kidney.DL <- combine.study.effects(all.effect.mini.kidney, between.method="DL", everything=F, verbose=T)
    mini.kidney.SJ <- combine.study.effects(all.effect.mini.kidney, between.method="SJ", everything=F, verbose=T)
    mini.kidney.HS <- combine.study.effects(all.effect.mini.kidney, between.method="HS", everything=F, verbose=T)
    
    list(DL=mini.kidney.DL, SJ=mini.kidney.SJ, 
         HS=mini.kidney.HS, sizes=classLengths)
})  
names(out.list) <- make.names(1:length(combos))

saveRDS(out.list, file=paste0("./data/mini.kidney.reps.out.5subsets.OVERLAPGENES.RDS"))
