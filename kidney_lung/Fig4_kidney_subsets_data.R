

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
cuts.n <- cut(combo.n, breaks=seq(0, cutoff, 50))
n.datasets <- unlist(lapply(combos, length))

table(cuts.n, n.datasets)
#hist(combo.n[cuts.n=="(200,300]"])

## The counts of datasets above and below the median for a given threshold
above.below <- t(sapply(levels(cuts.n), function(level){
    med.sets <- median(n.datasets[cuts.n==level], na.rm=T)
    c(sum(n.datasets[cuts.n==level]<med.sets, na.rm=T),
      sum(n.datasets[cuts.n==level]>med.sets, na.rm=T))
}))

combos.to.test <- lapply(levels(cuts.n), function(level){
    med.sets <- median(n.datasets[cuts.n==level], na.rm=T)
    list(under.med = combos[which(n.datasets<med.sets & cuts.n==level)],
         over.med = combos[which(n.datasets>med.sets & cuts.n==level)])
})
names(combos.to.test) <- make.names(levels(cuts.n))


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


out.list <- list()
##For each dataset level
for(i in 1:length(combos.to.test)){
    
    ## for both above and below median
    for(j in 1:2){
        cat(i, j, "\t")
        combos <- combos.to.test[[i]][[j]]
        
        #parallelize actual analyses
        out <- mclapply(combos, mc.cores=maxCores, function(combo){
            #cat(". ")
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
            
            #clean memory between runs... who knows?
            rm(kidney.subset, all.ES, all.effect.mini.kidney)
            
            list(DL=mini.kidney.DL, SJ=mini.kidney.SJ, 
                 HS=mini.kidney.HS, sizes=classLengths)
        })
        names(out) <- make.names(0:length(combos))[-1]
        ab.bel <- c("below.med", "above.med")
        out.list[[ names(combos.to.test)[i] ]][[ ab.bel[j] ]] <- out
    }
}

saveRDS(out.list, file="./data/mini.kidney.reps.out.N.subsets.OVERLAPGENES.RDS")
