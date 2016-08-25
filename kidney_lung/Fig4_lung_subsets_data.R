

require(methods)
source("meta_compare.functions.R")

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
overlap <- Reduce(intersect, lapply(lung.results.gems, function(x) x$keys))
overlap <- overlap[!is.na(overlap)]
lung.results.gems <- lapply(lung.results.gems, function(gem) {
    gem$keys <- gem$keys[gem$keys %in% overlap]
    gem$expr <- gem$expr[names(gem$keys), ]
    gem})


out.list <- list()
##For each dataset level
for(i in 1:length(combos.to.test)){
    
    ## for both above and below median
    for(j in 1:2){
        cat(i, j, "\t")
        combos <- combos.to.test[[i]][[j]]
        
        #parallelize actual analyses
        library(parallel)
        maxCores <- detectCores()  ## for SLURM, else use 10
        out <- mclapply(combos, mc.cores=12, function(combo){
            #cat(". ")
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
            
            #clean memory between runs... who knows?
            rm(lung.subset, all.ES, all.effect.mini.lung)
            
            list(DL=mini.lung.DL, SJ=mini.lung.SJ, 
                 HS=mini.lung.HS, sizes=classLengths)
        })
        names(out) <- make.names(1:length(combos))
        ab.bel <- c("below.med", "above.med")
        out.list[[ names(combos.to.test)[i] ]][[ ab.bel[j] ]] <- out
    }
}

saveRDS(out.list, file="./data/mini.lung.reps.out.N.subsets.OVERLAPGENES.RDS")
