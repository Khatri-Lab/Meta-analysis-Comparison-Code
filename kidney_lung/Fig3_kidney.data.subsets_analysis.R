
setwd("/projects/Sepsis_Tim/metaComparison/data/")
source("../meta.compare.functions.R")

kidney.subsets <- readRDS("./mini.kidney.reps.out.5subsets.OVERLAPGENES.RDS")
# mean(unlist(lapply(kidney.subsets, function(x) sum(x$sizes))))

load("./meta.compared.kidney.OVERLAPGENES.RData")
kidney.true.q.01 <- Reduce(intersect, lapply(list(kidney.SJ, kidney.DL, kidney.HS), function(x){
  rownames(subset(x, ES.FDR<0.01))    #& n.studies==max(n.studies)
}))

 

#########  RANGEs ES, Het, Q  ##################
getTPFPpointsRange <- function(testRange){
    points <- sapply(testRange, function(test){
        with(test, c(SJ.FP=mean(SJ.FP), SJ.TP=mean(SJ.TP),
                     DL.FP=mean(DL.FP), DL.TP=mean(DL.TP),
                     HS.FP=mean(HS.FP), HS.TP=mean(HS.TP)) )  
    })  
    t(points)
}


kidney.hetrange <- mclapply(c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.5), mc.cores=4, function(het){
    applyTestChar(kidney.subsets, kidney.true.q.01, 
                  filt="het", het=het, all.studies=T)
})
kidney.het.points <- getTPFPpointsRange(kidney.hetrange)


kidney.ESrange <- mclapply(seq(1,2,0.1), mc.cores=5, function(ES){
    applyTestChar(kidney.subsets, kidney.true.q.01, 
                  filt="ES", ES=ES, all.studies=T)
})
kidney.ES.points <- getTPFPpointsRange(kidney.ESrange)    


kidney.Qrange <- mclapply(10^(-1*seq(2,10,1)), mc.cores=5, function(Q){
    applyTestCharSig(kidney.subsets, kidney.true.q.01, 
                     sig="fdr", level=Q, all.studies=T)
})
kidney.Q.points <- getTPFPpointsRange(kidney.Qrange)


kidney.ESrange.qrange <- mclapply(c(0.001, 0.01, 0.05, 0.1, 0.2), mc.cores=5, function(q){
  lapply(seq(1,1.5,0.1), function(ES){
    applyTestChar(kidney.subsets, kidney.true.q.01, 
                  filt="ES", ES=ES, fdr=q, all.studies=T)
  })
})
kidney.ES.Q.pointsList <- lapply(kidney.ESrange.qrange, getTPFPpointsRange)
 
#### WITH INDETERMINATE: 
# save(kidney.het.points, kidney.ES.points, kidney.Q.points, kidney.ES.Q.pointsList, 
#      file="kidney_5_samples_points.RData")
## skip calcs 
# load("kidney_5_samples_points.RData")


###########   plots   #################################
# pdf("./output/5-dataset_subsamples_kidney_noIndet.pdf", width=6, height=12, useDingbats=F)
par(mfrow=c(2,1))

ylim <- max(kidney.Q.points[,6])
xlim <- max(kidney.Q.points[,5])
plot(1, col="white", pch=19, #log="xy",
     ylim=c(1,ylim), xlim=c(1,xlim), 
     xlab="N, FP (log scale)", ylab="N, TP (log scale)", 
     main="kidney \n5-dataset sub-samples")
#lines(c(0.1,5000), c(0.1,5000), lty=4, col=1)
points(kidney.Q.points[1, c(1,3,5)], kidney.Q.points[1, c(2,4,6)], pch=19, col=c("gold", "green", "skyblue"))
lines(kidney.ES.points[,1:2], col="gold", lty=1, lwd=2)
lines(kidney.ES.points[,3:4], col="green", lty=1, lwd=2)
lines(kidney.ES.points[,5:6], col="skyblue", lty=1, lwd=2)
lines(kidney.het.points[,1:2], col="gold", lty=3, lwd=2)
lines(kidney.het.points[,3:4], col="green", lty=3, lwd=2)
lines(kidney.het.points[,5:6], col="skyblue", lty=3, lwd=2)
lines(kidney.Q.points[,1:2], col="gold", lty=2, lwd=2)
lines(kidney.Q.points[,3:4], col="green", lty=2, lwd=2)
lines(kidney.Q.points[,5:6], col="skyblue", lty=2, lwd=2)
leg <- legend("bottomright", legend=c("Effect Size (1-2 fold)","Significance (q<0.01 - q<1e-10)","Heterogeneity (p>0 - p>0.5)"),
              lty=1:3, lwd=2, col=1,
              bg="white", cex=1.1, bty="n", y.intersp=0.8) 
legend(x=leg$rect$left+leg$rect$w/3, y=leg$rect$top+leg$rect$h, legend=c("Sidik-Jonkman","DerSimonian-Laird","Hunter-Schmidt"),
       fill=c("gold", "green", "skyblue"), border=c("gold", "green", "skyblue"),
       bg="white", cex=1.1, bty="n", y.intersp=0.8) 
legend(x=leg$rect$left+leg$rect$w*.66, y=leg$rect$top+1.4*leg$rect$h, legend="q<0.01", pch=19,bg="white", cex=1.1, bty="n", y.intersp=0.8)

# dev.off()


#########  DL only
plot(1, col="white", pch=19, #log="xy",
     ylim=c(200,1600), xlim=c(0,800), 
     xlab="N, FP (log scale)", ylab="N, TP (log scale)", 
     main="kidney \n5-dataset sub-samples")
#lines(c(0.1,5000), c(0.1,5000), lty=5, col=1)
k<-1; lapply(kidney.ES.Q.pointsList, function(tmp){
  points(tmp[1, 3], tmp[1, 4], 
         pch=20, col=c("green"))
  text(tmp[1, 3], tmp[1, 4], 
       labels=paste0("q<", c(0.001, 0.01, 0.05, 0.1, 0.2))[k], 
       col=c("black"),
       cex=0.9, adj=c(-0.25, 0.7))
  lines(tmp[,3:4], col=rgb(0,1,0,0.2, alpha=0.3), lty=1, lwd=2.5)
  k <<- k+1
  for(i in 3:nrow(tmp)){
    row <- tmp[i,]
    ES <- paste(seq(1,1.5,0.1))[i]
    text(row[3], row[4], 
         labels=ES, col=c("black"),
         cex=0.8)        
  }
})
text(600, 900, "Numbers indicate\n effect size limits\n for given q value", cex=1.1)
legend("bottomright", legend="DerSimonian-Laird",
       fill="green",  cex=1.1, bty="n", y.intersp=0.7) 

dev.off()

