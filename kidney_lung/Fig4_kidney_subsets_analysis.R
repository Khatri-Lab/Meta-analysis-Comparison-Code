

source("../meta_compare.functions.R")

load("meta.compared.kidney.OVERLAPGENES.RData")
kidney.true.q.01 <- Reduce(intersect, lapply(list(kidney.SJ, kidney.DL, kidney.HS), function(x){
    rownames(subset(x, ES.FDR<0.01 & n.studies==max(n.studies))) #n.studies==max(n.studies) redundant for OVERLAPGENES                     
}) )


### HUGE FILES  
### 30Gb uncompressed -- takes 5 min to load
kidney.N.subsets <- readRDS("./mini.kidney.reps.out.N.subsets.OVERLAPGENES.RDS")
kidney.N.median <- readRDS("./mini.kidney.reps.out.N.subsets_MEDIAN.OVERLAPGENES.RDS")  


##########  find sqrtprod sample size 
kidneyObject <- readRDS("./kidneyResults.RDS")

getSqrtProd <- function(X){
    GSEs <- names(X$sizes)
    sqrtprod <- sum(unlist(lapply(GSEs, function(GSEname) {
            counts <- table(kidneyObject$gems[[GSEname]]$class)
            sqrt(prod(counts))
    })))
    X$sqrtprod <- sqrtprod
    return(X)
}

kidney.N.subsets.all <- kidney.N.subsets
kidney.N.median.all <- kidney.N.median

kidney.N.subsets <- lapply(kidney.N.subsets.all, function(N){
    lapply(N, function(ab.bel){
        lapply(ab.bel, function(X) getSqrtProd(X))
    })
})

kidney.N.median <- lapply(kidney.N.median.all, function(N){
    lapply(N, function(X) getSqrtProd(X))
})


###############  get individual stats  #########################
getDLstats.indiv <- function(N.subsets){
    if(length(N.subsets)==0) return(NULL)
    stats <- applyTestChar(N.subsets, kidney.true.q.01, 
                           filt="ES", ES=1.3, fdr=0.01, all.studies=T)
    DL.stats <- stats[, c("DL.TP", "DL.FP")]
    FPR <- apply(DL.stats, 1, function(x) x[2]/sum(x))
    N.samples <- unlist(lapply(N.subsets, function(x) sum(x$sizes)))
    N.datasets <- unlist(lapply(N.subsets, function(x) length(x$sizes)))
        ## remove if not done above
    sqrtprod <- unlist(lapply(N.subsets, function(x) x$sqrtprod)) 
    cbind(N.samples=N.samples, N.datasets=N.datasets, DL.stats, FPR=FPR, sqrtprod=sqrtprod)
}

below.indiv <- lapply(kidney.N.subsets, function(Nrange){
    getDLstats.indiv(Nrange$below.med)
})

above.indiv <- lapply(kidney.N.subsets, function(Nrange){
    getDLstats.indiv(Nrange$above.med)
})
median.indiv <- lapply(kidney.N.median, getDLstats.indiv)

below.indiv <- do.call(rbind, below.indiv)
above.indiv <- do.call(rbind, above.indiv)
median.indiv <- do.call(rbind, median.indiv)

indiv <- rbind(below.indiv, median.indiv, above.indiv)
indiv <- subset(indiv, !is.na(DL.TP))
indiv <- indiv[order(indiv$N.datasets),]

indiv$DL.TN <- nrow(kidney.DL) - length(kidney.true.q.01) - indiv$DL.FP
indiv$DL.FN <- nrow(kidney.DL) - indiv$DL.TP - indiv$DL.TN - indiv$DL.FP
indiv$TPR <- indiv$DL.TP/(indiv$DL.TP + indiv$DL.FN)
indiv$accuracy <- (indiv$DL.TP + indiv$DL.TN)/nrow(kidney.DL)

# saveRDS(indiv, file="./output/kidney_indiv_OVERLAPGENES.RDS")

########## plot individual points by K datasets ##################

library(RColorBrewer)
indiv <- indiv[!indiv$N.datasets>6, ]

############  N vs accuracy  ##############################
# pdf("./output/kidney_K.vs.N_ACCURACY_lines.pdf", useDingbats=F, width=7, height=6)
indiv <- indiv[order(indiv$N.datasets),]
Kpal <- brewer.pal(7, "Spectral")[indiv$N.datasets]
with(indiv, plot(N.samples, accuracy, pch=20, col=Kpal,  
                 main="kidney dataset subsamples\nDersimonian-Laird q<0.01 & ES>1.3-fold"))
for(k in unique(indiv$N.datasets)){
  temp.K <- subset(indiv, N.datasets==k)
  model <- lm(accuracy~N.samples, temp.K)
  X <- c(min(temp.K$N.samples), max(temp.K$N.samples))
  Y <- predict(model, newdata=data.frame(N.samples=X))
  lines(x=X, y=Y, col="grey90", lwd=5)
  lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}
legend("topleft", legend=unique(indiv$N.datasets), fill=unique(Kpal), 
       title="K datasets", bg="white")

summary(lm(accuracy~N.datasets+N.samples, data=indiv))

dev.off()


######## ratio of N to geomMean - further right = more unbalanced
# pdf("./output/kidney_K.vs.NGMratio.pdf", useDingbats=F, width=7, height=6)
indiv$ratio <- with(indiv, N.samples/(sqrtprod*2))
with(indiv, plot(ratio, accuracy, pch=20, col=Kpal,  
                 main="kidney dataset subsamples\nDersimonian-Laird q<0.01 & ES>1.3-fold"))
for(k in unique(indiv$N.datasets)){
    temp.K <- subset(indiv, N.datasets==k)
    model <- lm(accuracy~ratio, temp.K)
    X <- with(temp.K, c(min(ratio), max(ratio)))
    Y <- predict(model, newdata=data.frame(ratio=X))
    lines(x=X, y=Y, col="grey90", lwd=5)
    lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}
legend("topright", legend=unique(indiv$N.datasets), fill=unique(Kpal), 
       title="K datasets", bg="white")

dev.off()

############  sqrtprod vs accuracy  ##############################
# pdf("./output/kidney_K.vs.sqrtprod_ACCURACY_lines_K2.pdf", useDingbats=F, width=7, height=6)
indiv <- indiv[order(indiv$N.datasets),]
Kpal <- brewer.pal(7, "Spectral")[indiv$N.datasets]
with(indiv, plot(sqrtprod, accuracy, pch=20, col=Kpal,  
                 main="kidney dataset subsamples\nDersimonian-Laird q<0.01 & ES>1.3-fold"))
for(k in unique(indiv$N.datasets)){
    temp.K <- subset(indiv, N.datasets==k)
    model <- lm(accuracy~sqrtprod, temp.K)
    X <- c(min(temp.K$sqrtprod), max(temp.K$sqrtprod))
    Y <- predict(model, newdata=data.frame(sqrtprod=X))
    lines(x=X, y=Y, col="grey90", lwd=5)
    lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}
legend("topleft", legend=unique(indiv$N.datasets), fill=unique(Kpal), 
       title="K datasets", bg="white")

summary(lm(accuracy~N.datasets+sqrtprod, data=indiv))

dev.off()


#############  FPR and TPR  ######################
# pdf("./output/kidney_K.vs.N_TP.FP_NOlines.pdf", useDingbats=F, width=7, height=9)
par(mfrow=c(2,1), mar=c(2,4,4,7), xpd=T)
indiv <- indiv[order(indiv$N.datasets),]
with(indiv, plot(N.samples, DL.FP, pch=20, col=Kpal, ylab="False Positives", 
                 main="kidney transplant dataset subsamples\nDersimonian-Laird q<0.01 & ES>1.3-fold"))
for(k in unique(indiv$N.datasets)){
    temp.K <- subset(indiv, N.datasets==k)
    model <- lm(DL.FP~N.samples, temp.K)
    X <- c(min(temp.K$N.samples), max(temp.K$N.samples))
    Y <- predict(model, newdata=data.frame(N.samples=X))
    #lines(x=X, y=Y, col="grey90", lwd=5)
    #lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}
legend("topright", inset=c(-0.2,0), legend=unique(indiv$N.datasets), fill=unique(Kpal), 
       title="K datasets", bg="white")

par(mar=c(4,4,2,7), xpd=F)
with(indiv, plot(N.samples, DL.TP, pch=20, col=Kpal, ylab="True Positives", xlab="N (samples)"))
for(k in unique(indiv$N.datasets)){
    temp.K <- subset(indiv, N.datasets==k)
    model <- lm(DL.TP~N.samples, temp.K)
    X <- c(min(temp.K$N.samples), max(temp.K$N.samples))
    Y <- predict(model, newdata=data.frame(N.samples=X))
    #lines(x=X, y=Y, col="grey90", lwd=5)
    #lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}

dev.off()


#############  FPR and TPR  ######################
# pdf("./output/kidney_K.vs.sqrtprod_TP.FP_NOlines.pdf", useDingbats=F, width=7, height=9)
par(mfrow=c(2,1), mar=c(2,4,4,7), xpd=T)
indiv <- indiv[order(indiv$N.datasets),]
with(indiv, plot(sqrtprod, DL.FP, pch=20, col=Kpal, ylab="False Positives", 
                 main="kidney transplant dataset subsamples\nDersimonian-Laird q<0.01 & ES>1.3-fold"))
for(k in unique(indiv$N.datasets)){
    temp.K <- subset(indiv, N.datasets==k)
    model <- lm(DL.FP~sqrtprod, temp.K)
    X <- c(min(temp.K$sqrtprod), max(temp.K$sqrtprod))
    Y <- predict(model, newdata=data.frame(sqrtprod=X))
    #lines(x=X, y=Y, col="grey90", lwd=5)
    #lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}
legend("topright", inset=c(-0.2,0), legend=unique(indiv$N.datasets), fill=unique(Kpal), 
       title="K datasets", bg="white")

par(mar=c(4,4,2,7), xpd=F)
with(indiv, plot(sqrtprod, DL.TP, pch=20, col=Kpal, ylab="True Positives", xlab="Sum of Geometric Mean of Class Size"))
for(k in unique(indiv$N.datasets)){
    temp.K <- subset(indiv, N.datasets==k)
    model <- lm(DL.TP~sqrtprod, temp.K)
    X <- c(min(temp.K$sqrtprod), max(temp.K$sqrtprod))
    Y <- predict(model, newdata=data.frame(sqrtprod=X))
    #lines(x=X, y=Y, col="grey90", lwd=5)
    #lines(x=X, y=Y, col=brewer.pal(7, "Spectral")[k], lwd=3)
}

dev.off()
