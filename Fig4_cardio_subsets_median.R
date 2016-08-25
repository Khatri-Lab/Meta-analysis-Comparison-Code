

require(methods)
source("meta_compare.functions.R")

## smaller object, ease of reading:
cardioResults <-readRDS("cardioResults.RDS")
cardio.results.gems <- cardioResults$gems
rm(cardioResults)
data.n <- unlist(lapply(cardio.results.gems, function(gem) length(gem$class)))

## obtain all possible combinations of subsets
nstudies <- length(cardio.results.gems)
combos <- lapply(1:nstudies, function(n) t(combn(nstudies, n, simplify=F)))
combos <- unlist(combos, recursive = F)

# examine sizes of all subsets
combo.n <- unlist(lapply(combos, function(comb) sum(data.n[comb])))
#hist(combo.n, plot=T)
cutoff <- sum(data.n)/2
cutoff <- ceiling(cutoff/50)*50
cuts.n <- cut(combo.n, breaks=seq(0, cutoff, 50)); rm(cutoff)
n.datasets <- unlist(lapply(combos, length))

##get median
combos.median <- lapply(levels(cuts.n), function(level){
    med.sets <- median(n.datasets[cuts.n==level], na.rm=T)
    combos[which(n.datasets==med.sets & cuts.n==level)]
})
names(combos.median) <- make.names(levels(cuts.n))

##################  limit gene set  ########################
## if modify, change name of SAVERDS file, too
overlap <- Reduce(intersect, lapply(cardio.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
cardio.results.gems <- lapply(cardio.results.gems, function(gem) {
    gem$keys <- gem$keys[gem$keys %in% overlap]
    gem$expr <- gem$expr[names(gem$keys), ]
    gem})


### Set number of cores. If don't want to use all resources, 
## need to enter something manually (ie maxCores <- 10)
library(parallel)
slurm <- Sys.getenv("SLURM_NPROCS")
maxCores <- ifelse(slurm=="", detectCores(), slurm)


#### 
out.list <- lapply(combos.median, function(combos){
    
    # parallelize inner
    out <- mclapply(combos, mc.cores=maxCores, function(combo){
        
        #Take each subset of 5 datasets, store class lengths
        cardio.subset <- cardio.results.gems[unlist(combo)]
        classLengths <- unlist(lapply(cardio.subset , function(gem) length(gem$class)))
        
        #Get initial effect sizes, then combine probes -> genes
        #lapply(cardio.subsets, function(gem) gem$class)
        all.ES <- lapply(cardio.subset, effect.sizes ) 
        all.effect.mini.cardio <- study.effects(all.ES)
        
        mini.cardio.DL <- combine.study.effects(all.effect.mini.cardio, between.method="DL", everything=F, verbose=T)
        mini.cardio.SJ <- combine.study.effects(all.effect.mini.cardio, between.method="SJ", everything=F, verbose=T)
        mini.cardio.HS <- combine.study.effects(all.effect.mini.cardio, between.method="HS", everything=F, verbose=T)
        
        list(DL=mini.cardio.DL, 
             SJ=mini.cardio.SJ, 
             HS=mini.cardio.HS, 
             sizes=classLengths)
    })      
    
})

names(out.list) <- names(combos.median)

saveRDS(out.list, file=paste0("./data/mini.cardio.reps.out.N.subsets_MEDIAN.OVERLAPGENES.RDS"))
