

require(methods)
source("./meta_compare.functions.R")

## smaller object, ease of reading:
kidneyResults <-readRDS("kidneyResults.RDS")
kidney.results.gems <- kidneyResults$gems
rm(kidneyResults)
data.n <- unlist(lapply(kidney.results.gems, function(gem) length(gem$class)))

## obtain all possible combinations of subsets
nstudies <- length(kidney.results.gems)
combos <- lapply(1:nstudies, function(n) t(combn(nstudies, n, simplify=F)))
combos <- unlist(combos, recursive = F)

# examine sizes of all subsets
combo.n <- unlist(lapply(combos, function(comb) sum(data.n[comb])))
#hist(combo.n, plot=T)
cutoff <- sum(data.n)/2
cutoff <- round(cutoff/50)*50
cuts.n <- cut(combo.n, breaks=seq(0, cutoff, 50)); rm(cutoff)
n.datasets <- unlist(lapply(combos, length))

table(cuts.n, n.datasets)
      
##get median
combos.median <- lapply(levels(cuts.n), function(level){
    med.sets <- median(n.datasets[cuts.n==level], na.rm=T)
    combos[which(n.datasets==med.sets & cuts.n==level)]
})
names(combos.median) <- make.names(levels(cuts.n))


##################  limit gene set  ########################
## if modify, change name of SAVERDS file, too
overlap <- Reduce(intersect, lapply(kidney.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
kidney.results.gems <- lapply(kidney.results.gems, function(gem) {
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
    
    out <- mclapply(combos, mc.cores=maxCores, function(combo){
        
        #Take each subset of 5 datasets, store class lengths
        kidney.subset <- kidney.results.gems[unlist(combo)]
        classLengths <- unlist(lapply(kidney.subset , function(gem) length(gem$class)))
        
        #Get initial effect sizes, then combine probes -> genes
        #lapply(kidney.subsets, function(gem) gem$class)
        all.ES <- lapply(kidney.subset, effect.sizes ) 
        all.effect.mini.kidney <- study.effects(all.ES)
        
        mini.kidney.DL <- combine.study.effects(all.effect.mini.kidney, between.method="DL", everything=F, verbose=T)
        mini.kidney.SJ <- combine.study.effects(all.effect.mini.kidney, between.method="SJ", everything=F, verbose=T)
        mini.kidney.HS <- combine.study.effects(all.effect.mini.kidney, between.method="HS", everything=F, verbose=T)
        
        list(DL=mini.kidney.DL, 
             SJ=mini.kidney.SJ, 
             HS=mini.kidney.HS, 
             sizes=classLengths)
    })      
    
})

names(out.list) <- names(combos.median)

saveRDS(out.list, file=paste0("./data/mini.kidney.reps.out.N.subsets_MEDIAN.OVERLAPGENES.RDS"))
