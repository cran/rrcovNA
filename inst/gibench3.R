gibench3 <- function(R=10, meth=c("mcd", "ogk", "S", "sde", "sign1",
                                  "EA", "TRC", "BEM",
                                  "rseq_ogk", "rseq_mcd", "rseq_sde", "rseq_S", "rseq_sign1", "rseq_sign2", "rseq_pcout"),
                                  miss=seq(0.0, 0.4, 0.1),
                                  seed=36)
{


    library(rrcovNA)

    set.seed(seed)
    data(bushfire)
    x <- as.matrix(bushfire)

    ## variiere missings:
    ptm <- proc.time()
    erg <- mySim(x=x, meth=meth, R=R, miss=miss)   ## increase R!
    FNrate <- (erg[,"TO"] - erg[,"TN"])/erg[,"TO"]
    FPrate <- erg[,"FP"]/(erg[,"N"] - erg[,"TO"])
    erg <- cbind(erg, FNrate, FPrate)

    erg1 <- aggregate(erg[,"FNrate"], list(factor(erg$method), factor(erg$miss)), mean)     # False negative
    erg1[,3] <- round(100*erg1[,3],2)
    names(erg1) <- c("Method", "Missing", "FNrate")
    erg2 <- aggregate(erg[,"FPrate"], list(factor(erg$method), factor(erg$miss)), mean)     # False positive
    erg2[,3] <- round(100*erg2[,3],2)
    names(erg2) <- c("Method", "Missing", "FPrate")

    ergm <- cbind(erg1, erg2[,3])
    names(ergm) <- c("Method", "Missing", "FNrate", "FPrate")
    save(erg, ergm, file="erg-bush.Rda")

    xtime1 <- (proc.time() - ptm)[1]
    cat("\nElapsed time (missings):", format(xtime1, nsmall=2),"\n\n")
    return(invisible(ergm))
}

## Uses:
##  - no, rbern() replaced by rbinom(): Rlab for  rbern() - Bernoulli random variables
##  - rrcovNA --> rrcov
##            --> norm
##  NOTE:
##  - EA, TRC and BEM are currently not implemented
##
##
##  Creates a data set with missing values (for testing purposes)
##  from a complete data set 'x'. The probability of
##  each item being missing is 'pr' (Bernoulli trials).
##
getmiss <- function(x, pr=0.1){
##    library(Rlab)
    n <- nrow(x)
    p <- ncol(x)
    done <- FALSE
    iter <- 0
    while(iter <= 50){
        bt <- rbinom(n*p, 1, pr)        # bt <- rbern(n*p, pr)
        btmat <- matrix(bt, nrow=n)
        btmiss <- ifelse(btmat==1, NA, 0)
        y <- x+btmiss
        if(length(which(rowSums(nanmap(y)) == p)) == 0)
            return (y)
        iter <- iter + 1
    }
    y
}

##  Returns a 1/0 map of the missing values in x
nanmap <- function(x)
    ifelse(is.na(x),1,0)

##
## Outliers bestimmen
##  ind - list of known outliers (usually returned by
##          genDat() as m$indexOutliers)
##  m1  - list of identified outliers (returned by applyMethod)
##
createMeasures <- function(ind, m1, xtime) {

    out <- length(ind)                       #  number of known outliers
    outIdent <- length(m1)                   #  number of identified outliers
    outTrue <- length(which(ind %in% m1))    # number of correctly identified outliers
    outFalse <- outIdent - outTrue           #
    Performance <- 9999                      # ev. noch machen, wenn Zeit uebrig

    return(matrix(c(out, outIdent, outTrue, outFalse, Performance, xtime), ncol=6))
}

### Matrix erstellen fuer Simulationsergebnisse
##  The first three columns are Method, Outlier Fraction, Missings percentage
##  The next 5 columns are for the corresponding results - will be filled in
##  by createMeasures()
##
createMatrix <- function(method, out, miss, n, R) {
    V1 <- rep(method, length=R)
    V2 <- rep(out, length=R)
    V3 <- rep(miss, length=R)
    V4 <- rep(n, length=R)
    V5 <- V6 <- V7 <- V8 <- V9 <- V10 <- rep(NA, length=R)
    d <- data.frame(V1,V2,V3,V4,V5,V6,V7,V8,V9, V10)

    ## N    - Number of observations
    ## TO   - number of true (known, generated) outlier
    ## PO   - number predicted outliers = TN + FN
    ## TN   - true negatives the number of outliers that were
    ##          correctly identified as outlier
    ## FN   - false negatives = TO-TN - the outliers that
    ##          were not identified, or masked outliers
    ## FP   - false positives - non-outliers that were classified as
    ##          outliers or swamped non-outliers
    ##
    ##  The following will be reported:
    ## ErrFN - outlier error rate = FN/TO
    ## ErrFP - inlier error rate  = FP/(n-TO)
    ##
    colnames(d) <- c("method", "out", "miss", "N", "TO",
                    "PO", "TN", "FP", "Perf", "Time")
    return(d)
}


## Simulieren
##  n - length of data
##  R - repeat R times (to average on)
##  d - input and output - the results matrix
##  meth, miss and out - arrays of methods, missing patterns
#3      and outlier fraction
doSim <- function(x, iout=c(7:11, 31, 32:38), d, R=1, meth, miss=0.1) {
    n <- nrow(x)
    p <- ncol(x)

    d2 <- createMatrix(meth=meth, out=0, miss=miss, n=n, R=R)
    for (i in 1:R) {
        m <- if(miss == 0) x else getmiss(x, pr=miss)                # add miss% MCAR missing values
        save(m, file="work-bush-data.rda")
##        cat("\niter=", i, "\n")
        xtime <- proc.time()
        m1 <- applyMethod(m, meth)

        if(!is.null(m1))
        {
            xtime <- (proc.time() - xtime)[1]
            d2[i, 5:ncol(d2)] <- createMeasures(iout, m1, xtime)
        }
    }
    d <- rbind(d, d2)
    return(d)
}

##  main loop: all methods, all missing patterns, all outlier patterns
##
mySim <- function(x, method, R=10, miss=0) {

    d <- NULL
    for (i in 1:length(method)) {

        cat("Wir beginnen mit Methode",method[i],"..... \n")
        flush.console()

        for (j in 1:length(miss)) {
            cat("\nSimulation fuer Methode", method[i], "miss=",miss[j],".....\n")
            d <- doSim(x, d=d, R=R, meth=method[i], miss=miss[j])
            cat("\nEnde!\n"); flush.console()
        }
        cat("\nFertig!\n"); flush.console()
    }
    cat("\nEnde der Simulation! \n")
    return(d)
}

## Apply the requested method on the data matrix x.
##  Returns an array of identified outliers
applyMethod <- function(x,
                    method=c("mcd", "sde", "ogk", "ogk-QC", "S",
                    "EA", "TRC", "BEM",
                    "sign1", "sign2", "pcout",
                    "rseq", "rseq_ogk", "rseq_mcd", "rseq_sde", "rseq_S", "rseq_sign1", "rseq_sign2", "rseq_pcout"))
{
    ## return the number of complete observations in x
    ncomplete <- function(x){
        na.x <- !is.finite(x %*% rep(1, ncol(x)))
        length(which(!na.x))
    }

    ## return a vector with the indices of the outliers respective to the Cov object 'x'
    getOutliers <- function(x)
        which(1.0 * getFlag(x) == 0)

    x <- as.matrix(x)
    method = match.arg(method)

    if(method=="EA")
    {
        warning(paste("Method currently not available: ", method))
        return(NULL)
## These are the default settin for EA which we obtained from the authors
##
##rrEA <- function(data,
##        weights = rep(1, nrow(data)),
##        reach="max",
##        transmission.function = "power",
##        power=ncol(data),
##        distance.type = "euclidean",
##        maxl = 5,
##        plotting = FALSE,
##        trace = FALSE,
##        prob.quantile = 0.9,
##        random.start = FALSE,
##       fix.start,
##        threshold=FALSE,
##        deterministic=FALSE)
##
        m1 <- rrEA(x,
            reach="max",                        # reach="max",
            transmission.function = "root",     # transmission.function = "power",
            power=ncol(x),                      # OK
            distance.type = "euclidean",        # OK
            maxl=5,                             # maxl = 5,
            prob.quantile = 0.9,
            random.start = FALSE,
            # fix.start,
            threshold=FALSE,
            deterministic=TRUE)
    }else if(method=="TRC")
    {
        warning(paste("Method currently not available: ", method))
        return(NULL)
        m1 <- try(rrTRC(x, prob.quantile=0.75, gamma=0.75, alpha=0.9, robust.regression="rank"))
        if(!is.numeric(m1))
        {
          cat("\nError in TRC - crash! Returning 0\n")
          save(x, file="work-err-TRC.rda")
          m1 <- c()
        }
    }else if(method == "BEM")
    {
        warning(paste("Method currently not available: ", method))
        return(NULL)
        m1 <- rrBEM(x)
    }else if(method == "mcd")
    {
        mcd <- CovNAMcd(x)
        m1 <- getOutliers(mcd)
    }else if(method == "ogk")
    {
        mcd <- CovNAOgk(x)
        m1 <- getOutliers(mcd)
    }else if(method == "ogk-QC")
    {
        ctrl <- CovControlOgk(smrob="s_mad", svrob="qc")
        mcd <- CovNAOgk(x, control=ctrl)
        m1 <- getOutliers(mcd)
    }else if(method == "sde")
    {
        mcd <- CovNASde(x)
        m1 <- getOutliers(mcd)
    }else if(method == "S")
    {
        mcd <- CovNASest(x, method="rocke")
        m1 <- getOutliers(mcd)
    }else if(method %in% c("sign1", "sign2", "pcout"))
    {
        ximp <- .imputation(x)        # impute missing data under the MLE
        out <- if(method == "sign1") sign1(ximp) else if(method == "sign2") sign2(ximp) else pcout(ximp)
        m1 <- which(out$wfinal01 == 0)
    }else if(method %in% c("rseq_sign1", "rseq_sign2", "rseq_pcout"))
    {
        ximp <- .imputation(x, impMeth="rseq", alpha=0.5)
        out <- if(method == "rseq_sign1") sign1(ximp) else if(method == "rseq_sign2") sign2(ximp) else pcout(ximp)
        m1 <- which(out$wfinal01 == 0)
    }else if(method == "rseq" | method == "rseq_mcd" )
    {
        mcd <- CovNAMcd(x, impMeth="rseq")
        m1 <- getOutliers(mcd)
    }else if(method == "rseq_ogk")
    {
        mcd <- CovNAOgk(x, impMeth="rseq")
        m1 <- getOutliers(mcd)
    }else if(method == "rseq_sde")
    {
        mcd <- CovNASde(x, impMeth="rseq")
        m1 <- getOutliers(mcd)
    }else if(method == "rseq_S")
    {
        mcd <- CovNASest(x, impMeth="rseq")
        m1 <- getOutliers(mcd)
    }else
        stop(paste("Method not defined: ", method))

    return(m1)
}
