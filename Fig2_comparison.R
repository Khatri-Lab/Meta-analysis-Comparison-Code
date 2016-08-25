

source("meta_compare.functions.R")

overlapGenes <- function(gems){
    overlap <- Reduce(intersect, lapply(gems, function(x) x$keys))
    overlap <- overlap[!is.na(overlap)]
    overlap.gems <- lapply(gems, function(gem) {
        gem$keys <- gem$keys[gem$keys %in% overlap]
        gem$expr <- gem$expr[names(gem$keys), ]
        gem})
    return(overlap.gems)
}

checkNumList <- function(SJ.EB.HE.DL.REML.HS){
  genes <- sapply(SJ.EB.HE.DL.REML.HS, checkNums) 
  colnames(genes) <- c("Sidik-Jonkman", "Empiric Bayes", "Hedges-Olkin", "DerSimonian-Laird",  
                       "Restricted Maximum Likelihood", "Hunter-Schmidt")
  rbind(genes, paste("Intersect ", 6:1), filteredIntersect(SJ.EB.HE.DL.REML.HS))
}



####################   lung Ca   #############################
lungResults <- readRDS("lungResults.RDS")
lung.results.gems <- overlapGenes(lungResults$gems)

lungResults <- runMetaAnalysis(lung.results.gems)

all.effect.lung <- study.effects(lungResults$all.ES)
unlist(lapply(lungResults$gems, function(x) length(x$class)))

lung.DL <- combine.study.effects(all.effect.lung, between.method="DL", everything=F)
lung.HE <- combine.study.effects(all.effect.lung, between.method="HE", everything=F)
lung.SJ <- combine.study.effects(all.effect.lung, between.method="SJ", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
lung.REML <- combine.study.effects(all.effect.lung, between.method="REML", everything=F,
                                   verbose=F, control=list(stepadj=0.5, maxiter=1000))
lung.HS <- combine.study.effects(all.effect.lung, between.method="HS", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
lung.EB <- combine.study.effects(all.effect.lung, between.method="EB", everything=F,
                                 verbose=F, control=list(stepadj=0.25, maxiter=500))

save(lung.DL, lung.HE, lung.HS, lung.SJ,
     lung.REML, lung.EB, file="meta.compared.lung.OVERLAPGENES.RData")


lungList <- list(lung.SJ, lung.EB, lung.HE, lung.DL, lung.REML, lung.HS)
genes <- checkNumList(lungList)
write.csv(genes, file="./output/genes.compared.lung.csv", quote=F)


####################   kidney xplant   #############################
kidneyResults <- readRDS("kidneyResults.RDS")
kidney.results.gems <- overlapGenes(kidneyResults$gems)

kidneyResults <- runMetaAnalysis(kidney.results.gems)

all.effect.kidney <- study.effects(kidneyResults$all.ES)
unlist(lapply(kidneyResults$gems, function(x) length(x$class)))

kidney.DL <- combine.study.effects(all.effect.kidney, between.method="DL", everything=F)
kidney.HE <- combine.study.effects(all.effect.kidney, between.method="HE", everything=F)
kidney.SJ <- combine.study.effects(all.effect.kidney, between.method="SJ", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
kidney.REML <- combine.study.effects(all.effect.kidney, between.method="REML", everything=F,
                                   verbose=F, control=list(stepadj=0.5, maxiter=1000))
kidney.HS <- combine.study.effects(all.effect.kidney, between.method="HS", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
kidney.EB <- combine.study.effects(all.effect.kidney, between.method="EB", everything=F,
                                 verbose=F, control=list(stepadj=0.25, maxiter=500))

save(kidney.DL, kidney.HE, kidney.HS, kidney.SJ,
     kidney.REML, kidney.EB, file="meta.compared.kidney.OVERLAPGENES.RData")


kidneyList <- list(kidney.SJ, kidney.EB, kidney.HE, kidney.DL, kidney.REML, kidney.HS)
genes <- checkNumList(kidneyList)
write.csv(genes, file="./output/genes.compared.kidney.csv", quote=F)


####################   cardiomyopathy   #############################
cardioResults <- readRDS("cardioResults.RDS")
cardio.results.gems <- overlapGenes(cardioResults$gems)
cardioResults <- runMetaAnalysis(cardio.results.gems)

all.effect.cardio <- study.effects(cardioResults$all.ES)
unlist(lapply(cardioResults$gems, function(x) length(x$class)))

cardio.DL <- combine.study.effects(all.effect.cardio, between.method="DL", everything=F)
cardio.HE <- combine.study.effects(all.effect.cardio, between.method="HE", everything=F)
cardio.SJ <- combine.study.effects(all.effect.cardio, between.method="SJ", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
cardio.REML <- combine.study.effects(all.effect.cardio, between.method="REML", everything=F,
                                     verbose=F, control=list(stepadj=0.5, maxiter=1000))
cardio.HS <- combine.study.effects(all.effect.cardio, between.method="HS", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
cardio.EB <- combine.study.effects(all.effect.cardio, between.method="EB", everything=F,
                                   verbose=F, control=list(stepadj=0.25, maxiter=500))

save(cardio.DL, cardio.HE, cardio.HS, cardio.SJ,
     cardio.REML, cardio.EB, file="meta.compared.cardio.OVERLAPGENES.RData")

cardioList <- list(cardio.SJ, cardio.EB, cardio.HE, cardio.DL, cardio.REML, cardio.HS)
genes <- checkNumList(cardioList)
write.csv(genes, file="./output/genes.compared.cardio.csv", quote=F)



########################   Sepsis   ##################################
load("/projects/Sepsis_Tim/Data/output/JointDiscovery.earlylate.NoNeut_discovery.RData")
sepsis.results.gems <- overlapGenes(discoveryResults.joint$gems)
sepsisResults <- runMetaAnalysis(sepsis.results.gems)

all.effect.sepsis <- study.effects(sepsisResults$all.ES)
sum(unlist(lapply(sepsisResults$gems, function(x) length(x$class))))

sepsis.DL <- combine.study.effects(all.effect.sepsis, between.method="DL", everything=F)
sepsis.HE <- combine.study.effects(all.effect.sepsis, between.method="HE", everything=F)
sepsis.SJ <- combine.study.effects(all.effect.sepsis, between.method="SJ", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
sepsis.REML <- combine.study.effects(all.effect.sepsis, between.method="REML", everything=F,
                                     verbose=F, control=list(stepadj=0.5, maxiter=500))
sepsis.HS <- combine.study.effects(all.effect.sepsis, between.method="HS", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
sepsis.EB <- combine.study.effects(all.effect.sepsis, between.method="EB", everything=F,
                                   verbose=F, control=list(stepadj=0.5, maxiter=500))

save(sepsis.DL, sepsis.HE, sepsis.HS, sepsis.SJ,
     sepsis.REML, sepsis.EB, file="./meta.compared.sepsis.OVERLAPGENES.RData")

sepsisList <- list(sepsis.SJ, sepsis.EB, sepsis.HE, sepsis.DL, sepsis.REML, sepsis.HS)
genes <- checkNumList(sepsisList)
write.csv(genes, file="./output/genes.compared.sepsis.csv", quote=F)


########################   TB   ##################################
load("/projects/Sepsis_Tim/TB/data/output/TB_OD.LTB.vs.ATB.perSamp_discovery.RData") 
TB.results.gems <- overlapGenes(discoveryResults_OD.LTB.vs.ATB$gems)
TBResults <- runMetaAnalysis(TB.results.gems)

all.effect.TB <- study.effects(TBResults$all.ES)
sum(unlist(lapply(TBResults$gems, function(x) length(x$class))))

TB.DL <- combine.study.effects(all.effect.TB, between.method="DL", everything=F)
TB.HE <- combine.study.effects(all.effect.TB, between.method="HE", everything=F)
TB.SJ <- combine.study.effects(all.effect.TB, between.method="SJ", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
TB.REML <- combine.study.effects(all.effect.TB, between.method="REML", everything=F,
                                 verbose=F, control=list(stepadj=0.5, maxiter=500))
TB.HS <- combine.study.effects(all.effect.TB, between.method="HS", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
TB.EB <- combine.study.effects(all.effect.TB, between.method="EB", everything=F,
                               verbose=F, control=list(stepadj=0.5, maxiter=500))

save(TB.DL, TB.HE, TB.HS, TB.SJ,
     TB.REML, TB.EB, file="./meta.compared.TB.OVERLAPGENES.RData")

TBList <- list(TB.SJ, TB.EB, TB.HE, TB.DL, TB.REML, TB.HS)
genes <- checkNumList(TBList)
write.csv(genes, file="./output/genes.compared.TB.csv", quote=F)


####################   Influenza   #############################
influenzaResults <- readRDS("./input/influenzaResults.RDS")
influenza.results.gems <- overlapGenes(influenzaResults$gems)
influenzaResults <- runMetaAnalysis(influenza.results.gems)

all.effect.influenza <- study.effects(influenzaResults$all.ES)
sum(unlist(lapply(influenzaResults$gems, function(x) length(x$class))))

# max(influenzaResults$pooled.ES[rownames(influenza.DL), "pval.het"] - influenza.DL$pval.het)

influenza.DL <- combine.study.effects(all.effect.influenza, between.method="DL", everything=F)
influenza.HE <- combine.study.effects(all.effect.influenza, between.method="HE", everything=F)
influenza.SJ <- combine.study.effects(all.effect.influenza, between.method="SJ", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
influenza.REML <- combine.study.effects(all.effect.influenza, between.method="REML", everything=F,
                                        verbose=F, control=list(stepadj=0.5, maxiter=500))
influenza.HS <- combine.study.effects(all.effect.influenza, between.method="HS", everything=F)
## Fisher scoring algorithm did not converge without reducting step length + increasing iters
influenza.EB <- combine.study.effects(all.effect.influenza, between.method="EB", everything=F,
                                      verbose=F, control=list(stepadj=0.5, maxiter=500))

save(influenza.DL, influenza.HE, influenza.HS, influenza.SJ,
     influenza.REML, influenza.EB, file="./meta.compared.influenza.OVERLAPGENES.RData")

influenzaList <- list(influenza.SJ, influenza.EB, influenza.HE, influenza.DL, influenza.REML, influenza.HS)
genes <- checkNumList(influenzaList)
write.csv(genes, file="./output/genes.compared.influenza.csv", quote=F)

