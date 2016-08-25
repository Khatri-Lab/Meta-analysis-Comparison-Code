

require(methods)
source("./meta_compare.functions.R")

## smaller object, ease of reading:
lungResults <-readRDS("lungResults.RDS")
lung.results.gems <- lungResults$gems
rm(lungResults)
data.n <- unlist(lapply(lung.results.gems, function(gem) length(gem$class)))

## obtain all possible combinations of subsets
nstudies <- length(lung.results.gems)
combos <- lapply(1:nstudies, function(n) t(combn(nstudies, n, simplify=F)))
combos <- unlist(combos, recursive = F)

# examine sizes of all subsets
combo.n <- unlist(lapply(combos, function(comb) sum(data.n[comb])))
#hist(combo.n, plot=T)
cutoff <- sum(data.n)/2
cutoff <- ceiling(cutoff/50)*50
cuts.n <- cut(combo.n, breaks=seq(100, cutoff, 50)); rm(cutoff)
n.datasets <- unlist(lapply(combos, length))

##get median
combos.median <- lapply(levels(cuts.n), function(level){
    med.sets <- median(n.datasets[cuts.n==level], na.rm=T)
    combos[which(n.datasets==med.sets & cuts.n==level)]
})
names(combos.median) <- make.names(levels(cuts.n))

##################  limit gene set  ########################
## if modify, change name of SAVERDS file, too
overlap <- Reduce(intersect, lapply(lung.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
lung.results.gems <- lapply(lung.results.gems, function(gem) {
    gem$keys <- gem$keys[gem$keys %in% overlap]
    gem$expr <- gem$expr[names(gem$keys), ]
    gem})


#### 
out.list <- lapply(combos.median, function(combos){
    
    library(parallel)
    maxCores <- detectCores()  ## for SLURM, else use 10
    out <- mclapply(combos, mc.cores=12, function(combo){
        
        #Take each subset of 5 datasets, store class lengths
        lung.subset <- lung.results.gems[unlist(combo)]
        classLengths <- unlist(lapply(lung.subset , function(gem) length(gem$class)))
        
        #Get initial effect sizes, then combine probes -> genes
        #lapply(lung.subsets, function(gem) gem$class)
        all.ES <- lapply(lung.subset, effect.sizes ) 
        all.effect.mini.lung <- study.effects(all.ES)
        
        mini.lung.DL <- combine.study.effects(all.effect.mini.lung, between.method="DL", everything=F, verbose=T)
        mini.lung.SJ <- combine.study.effects(all.effect.mini.lung, between.method="SJ", everything=F, verbose=T)
        mini.lung.HS <- combine.study.effects(all.effect.mini.lung, between.method="HS", everything=F, verbose=T)
        
        list(DL=mini.lung.DL, 
             SJ=mini.lung.SJ, 
             HS=mini.lung.HS, 
             sizes=classLengths)
    })      
    
})

names(out.list) <- names(combos.median)

saveRDS(out.list, file=paste0("./data/mini.lung.reps.out.N.subsets_MEDIAN.OVERLAPGENES.RDS"))
