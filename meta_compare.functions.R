


source("meta_analysis_functions.R")
require(rmeta)
require(metafor)

study.effects <- function (list.of.effects, within.method="fixed.iv"){
    if( is.null(names(list.of.effects)) )
        names(list.of.effects) <- paste("data", 1:length(list.of.effects), sep="")
    study.effects <- lapply(list.of.effects, function(effects) {
        
        effects <- data.frame(effects)
        effects$keys <- as.character(effects$keys)
        
        ## remove probes that cannot be mapped or have insufficient observations to calculate effect size
        bad <- which( is.na(effects$g) | is.na(effects$keys) | effects$keys=="NA" )
        effects <- effects[ setdiff(1:nrow(effects), bad), ]
        
        ## expand probes that maps to multiple keys
        effects <- expand.df( effects )
        
        ## summarize multiple probes within a study
        effects <- summ.eff.within(effects, option = within.method)
    })
    
    return(study.effects)
}

combine.study.effects <- function(study.effects, between.method="DL", everything=F, ...){
    tmp <- multimerge(study.effects)
    g    <- tmp[, paste(names(study.effects), "_g", sep = ""), drop=FALSE]
    se.g <- tmp[, paste(names(study.effects), "_se.g", sep = ""), drop=FALSE]
    
    pooled.estimates <- data.frame( pool.inverseVarMod(g, se.g, method=between.method, ... ) )
    
    if (everything) {
        return(list(g=g, se.g=se.g, pooled.estimates=pooled.estimates))
    } else {
        return(pooled.estimates)
    }
}

pool.inverseVarMod <- function(g, se.g, method, ...){   
    cat("Using modified function\n")
    
    stopifnot( identical( rownames(g), rownames(se.g) ) )
    out <- matrix( nr=nrow(g), nc=8,
                   dimnames=list( rownames(g), c("n.studies", "summary", "se.summary", "tau2", "p.value", "Q", "df", "pval.het") ) )
    
    for(j in 1:nrow(g)){
        e  <- cleanNA(    g[j, ] )
        se <- cleanNA( se.g[j, ] )
        n  <- length(e)
        if(n==0){
            #for debug
            print(g[j, ])
            cat("\n",rownames(g)[j], "\nNo finite measurements... returning NA")
            summ <- se.summ <- pval.het <- df.het <- Q.het <- tau2 <- NA
        } else if(n==1){
            summ <- e   
            se.summ <- se
            pval.het <- df.het <- Q.het <- tau2 <- NA
        } else {
            fit <- rma(yi=e, sei=se, method = method, ...)
            summ <- fit$b
            se.summ <- fit$se
            tau2 <- fit$tau2
            Q.het = fit$QE
            df.het = n-1
            pval.het = fit$QEp
            rm(fit)
        }
        pval     <- 2*pnorm( abs(summ/se.summ), lower.tail=FALSE )
        out[j, ] <- c(n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
        rm(e, se, n, summ, se.summ, tau2, pval, Q.het, df.het, pval.het)
    }
    out <- data.frame(out)
    out$ES.FDR <- p.adjust(out$p.value, method="fdr")
    return(out)
}

checkNums <- function(pooled){
    n.p <- sum(pooled$p.value<0.01)
    pct.p <- signif(n.p/nrow(pooled), 4)*100
    n.q <- sum(pooled$ES.FDR<0.01)
    pct.q <- signif(n.q/nrow(pooled), 4)*100
    n.filt.ES <- nrow(subset(pooled, ES.FDR<0.01 & summary>log2(1.5)))
    n.filt.ES.het <- nrow(subset(pooled, ES.FDR<0.01 & summary>log2(1.5) & pval.het>0.01))
    n.filt.ES.het.0.1 <- nrow(subset(pooled, ES.FDR<0.01 & summary>log2(1.5) & pval.het>0.1))
    return(c(n.p=n.p, pct.p=pct.p, n.q=n.q, pct.q=pct.q, n.filt.ES=n.filt.ES, n.filt.ES.het=n.filt.ES.het, n.filt.ES.het.0.1=n.filt.ES.het.0.1))
}

ESplots <- function(pooled, name, ylim){
    tmp <- subset(pooled, n.studies==max(n.studies))
    colors <- colorRampPalette(c("navyblue", "purple4", "red", "darkorange2"))(nrow(tmp))
    ggplot(tmp[order(tmp$pval.het), ], aes(x=summary, y=-log10(ES.FDR))) +
        scale_x_continuous(limits=c(-1.1, 1.1)) + 
        scale_y_continuous(limits=c(-0.1, ylim)) + 
        ggtitle(name) + 
        xlab("Summary effect size") + ylab("-log10 FDR") +
        geom_point(colour=colors, size=1)
}


filteredIntersect <- function(pooledList){
    genes <- lapply(pooledList, function(pooled){
        rownames(subset(pooled, ES.FDR<0.01))
    })
    cat("Number of analyses in common for filtered genes:\n")
    rev(table(table(unlist(genes))))
    
}



#############  ANALYSIS  ##############################
.testChar <- function(mini.method, true, filt, all.studies, ES, het, fdr){
    
    ifelse(all.studies, studies <- max(mini.method$n.studies), studies <- 1)
    
    if(filt=="EShet") {
        pos <- rownames(subset(mini.method, n.studies>=studies & ES.FDR<fdr & abs(summary)>log2(ES) & pval.het>het ))
    } else if(filt=="ES") {
        pos <- rownames(subset(mini.method, n.studies>=studies & ES.FDR<fdr & abs(summary)>log2(ES)))
    } else if(filt=="het") {
        pos <- rownames(subset(mini.method, n.studies>=studies & ES.FDR<fdr & pval.het>het))
    } else {
        pos <- rownames(subset(mini.method, n.studies>=studies & ES.FDR<fdr))
    }
    
    TP <- sum(pos %in% true) 
    FP <- sum(!(pos %in% true))
    c(TP=TP, FP=FP) ###Note TP not TPR
}

applyTestChar <- function(mini.out, true, filt=F, all.studies, ES=1.5, het=0.1, fdr=0.01){
    outTestChar <- mclapply(mini.out, mc.cores=10,  function(mini){
        if(is.null(mini)) return(c(SJ.TP=NA, SJ.FP=NA, DL.TP=NA, DL.FP=NA, HS.TP=NA, HS.FP=NA))
        SJ <- .testChar(mini$SJ, true, filt=filt, all.studies, ES=ES, het=het, fdr=fdr)
        DL <- .testChar(mini$DL, true, filt=filt, all.studies, ES=ES, het=het, fdr=fdr)
        HS <- .testChar(mini$HS, true, filt=filt, all.studies, ES=ES, het=het, fdr=fdr)
        c(SJ=SJ, DL=DL, HS=HS)
    })
    outTestChar <- data.frame(t(data.frame(outTestChar)))
    return(outTestChar)
}

.testCharSig <- function(mini.method, kidney.true, sig, level, all.studies=F){
    ifelse(all.studies, studies <- max(mini.method$n.studies), studies <- 1)
    
    if(sig=="fdr") {
        pos <- rownames(subset(mini.method, n.studies>=studies & ES.FDR<level ))
    } else if(sig=="bon") {
        mini.method$bon <- p.adjust(mini.method$p.value, method = "bonferroni")
        pos <- rownames(subset(mini.method, n.studies>=studies & bon<level ))
    } 
    
    TP <- sum(pos %in% kidney.true) 
    FP <- sum(!(pos %in% kidney.true))
    c(TP=TP, FP=FP) ###Note TP not TPR
}

applyTestCharSig <- function(mini.out, true, sig, level, all.studies=F){
    outTestChar <- mclapply(mini.out, mc.cores=10, function(mini){
        SJ <- .testCharSig(mini$SJ, true, sig, level, all.studies)
        DL <- .testCharSig(mini$DL, true, sig, level, all.studies)
        HS <- .testCharSig(mini$HS, true, sig, level, all.studies)
        c(SJ=SJ, DL=DL, HS=HS)
    })
    outTestChar <- data.frame(t(data.frame(outTestChar)))
    return(outTestChar)
}

##automagically switch to log axes
plotTPFP <- function(outTestCharTPFP, pct=NULL, plotlim, legPos="bottomright", filt, labels=T, subsets=F){
    outTestCharTPFP <- outTestCharTPFP + 1
    ifelse(subsets,
           main <- sprintf("%s 5-datasets subsets\n %s", nrow(outTestCharTPFP), filt),
           main <- sprintf("100 reps of %s%% subsets\n %s", pct, filt) )
    with(outTestCharTPFP, 
         plot(SJ.FP, SJ.TP, col=rgb(1,0.8,0,0.4), pch=19, log="xy",
              ylim=c(1,plotlim), xlim=c(1,plotlim), 
              xlab="N, FP (log scale)", ylab="N, TP (log scale)", 
              main=main) +
            points(DL.FP, DL.TP, col=rgb(0,1,0,0.2), pch=19) +
            points(HS.FP, HS.TP, col=rgb(0.3,1,1,0.4), pch=19) +
            abline(0,1, lty=2) + grid(col="grey") +
            symbols(mean(SJ.FP), mean(SJ.TP), inches=0.08, 
                     stars=t(rep(0.0001,4)), add=T, bg="gold") +
            symbols(mean(DL.FP), mean(DL.TP), inches=0.08, 
                     stars=t(rep(0.0001,4)), add=T, bg="green") +
            symbols(mean(HS.FP), mean(HS.TP), inches=0.08, 
                     stars=t(rep(0.0001,4)), add=T, bg="skyblue") )
    with(outTestCharTPFP, 
         legend(legPos, legend=c(sprintf("Sidik-Jonkman\n  FDR=%s \n", .calcFDRLabel(SJ.TP, SJ.FP) ),
                            sprintf("DerSimonian-Laird\n  FDR=%s \n", .calcFDRLabel(DL.TP, DL.FP) ),
                            sprintf("Hunter-Schmidt\n  FDR=%s \n", .calcFDRLabel(HS.TP, HS.FP) ) ),
           fill=c("gold", "green", "skyblue"), border=c("gold", "green", "skyblue"),
           bg="white", cex=1.2, bty="n", y.intersp=0.7) )
}

.calcFDRLabel <- function(TP, FP){
    # regress FP on (TP+FP) - force int to 0 - slope == FDR, then get std err
    pos <- TP + FP 
    tmp <- summary(lm(FP ~ 0 + pos))
    out <- signif(tmp$coefficients[1:2], 2)
    return(paste(out[1], "+/-", out[2]))
}


