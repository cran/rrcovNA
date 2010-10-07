setMethod("summary", "CovNARobust", function(object, ...){

    new("SummaryCovNARobust", covobj=object, evals=eigen(object@cov)$values)

})


setMethod("show", "SummaryCovNARobust", function(object){

    cat("\nCall:\n")
    print(object@covobj@call)

    digits = max(3, getOption("digits") - 3)
    cat("\nRobust Estimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nRobust Estimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(getEvals(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nRobust Distances: \n")
    dd <- getDistance(object)
    print.default(format(as.vector(dd$d), digits = digits), print.gap = 2, quote = FALSE)
})

setMethod("plot", signature(x="CovNARobust", y="missing"), function(x, y="missing",
                                which=c("all", "dd", "distance", "qqchi2", "tolEllipsePlot", "screeplot"),
                                classic= FALSE,
                                ask = (which=="all" && dev.interactive(TRUE)),
                                cutoff,
                                id.n,
                                tol = 1e-7, ...)
{
    data <- getData(x)
    ##  parameters and preconditions
    if(is.vector(data) || is.matrix(data)) {
        if(!is.numeric(data))
            stop( "x is not a numeric dataframe or matrix.")
    } else if(is.data.frame(data)) {
        if(!all(sapply(data,data.class) == "numeric"))
            stop( "x is not a numeric dataframe or matrix.")
    }

    n <- dim(data)[1]
    p <- dim(data)[2]

    if(length(getCenter(x))  == 0 ||  length(getCov(x)) == 0)
        stop( "Invalid object: attributes center and cov missing!")

    if(length(getCenter(x))  != p)
        stop( "Data set and provided center have different dimensions!")

    ## Check for singularity of the cov matrix
    if(rrcov:::isSingular(x))
        stop("The covariance matrix is singular!")

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(!missing(id.n) && !is.null(id.n)) {
        id.n <- as.integer(id.n)
        if(id.n < 0 || id.n > n)
            stop(sQuote("id.n")," must be in {1,..,",n,"}")
    }

    ccov <- CovNAClassic(data)
    print(getDistance(ccov))

    md <- rd <- NULL
    if(!rrcov:::isSingular(ccov))
        md <- sqrt(getDistance(ccov)$d)
    if(!rrcov:::isSingular(x))
        rd <- sqrt(getDistance(x)$d)

    which <- match.arg(which)
    op <- if (ask) par(ask = TRUE) else list()
    on.exit(par(op))

    mymiss <- .xmiss(x@X)

    ## distance-distance plot: here we need both robust and mahalanobis distances
    if((which == "all" || which == "dd") && !is.null(md) && !is.null(rd)) {
        .myddplotNA(md, rd, mymiss, cutoff=cutoff, id.n=id.n) # distance-distance plot
    }

    ## index plot of mahalanobis distances
    if((which == "all" || which == "distance")  && !is.null(rd)) {
        ylim <- NULL
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()

            ##VT::10.11.2007 - set same scale on both plots
            ylim <- c(min(rd,md), max(md,rd))
        }

        .mydistplotNA(rd, mymiss, cutoff, id.n=id.n, ylim=ylim)                   # index plot of robust distances
        if(classic && !is.null(md)) {
            .mydistplotNA(md, mymiss, cutoff, classic=TRUE, id.n=id.n, ylim=ylim) # index plot of mahalanobis distances
            par(opr)
        }
    }

    ## qq-plot of the mahalanobis distances versus the
    ## quantiles of the chi-squared distribution
    if((which == "all" || which == "qqchi2")  && !is.null(rd)) {
        if(classic && !is.null(md)) {
            opr <- if(prod(par("mfrow")) == 1) par(mfrow=c(1,2), pty="m") else list()
        }
        .qqplotNA(rd, p, mymiss, cutoff=cutoff, id.n=id.n) # qq-plot of the robust distances versus the
                                                 # quantiles of the chi-squared distribution
        if(classic && !is.null(md)) {
            .qqplotNA(md, p, mymiss, cutoff=cutoff, classic=TRUE, id.n=id.n)
                                                 # qq-plot of the mahalanobis distances
            par(opr)
        }
    }

##    if(which == "all" || which == "tolEllipsePlot") {
##        if(length(dim(data)) >= 2 && dim(data)[2] == 2){
##            if(!is.null(rd)){
##                if(classic &&  !is.null(md))
##                    .tolellipse(rcov=x, ccov = ccov, cutoff=cutoff, id.n=id.n, tol=tol, ...)
##                else
##                    .tolellipse(rcov=x, cutoff=cutoff, id.n=id.n, tol=tol, ...)
##            }
##        }
##        else if(which != "all")
##            warning("Warning: For tolerance ellipses the dimension must be 2!")
##    }

    if(which == "all" || which == "screeplot") {
        if(classic == TRUE)
            rrcov:::.myscreeplot(ccov=ccov, rcov=x)
        else
            rrcov:::.myscreeplot(rcov=x)
    }
}) ## end { plot("CovRobust") }

.labelNA <- function(x, y, xmiss, id.n=3) {
    if(id.n > 0) {
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        n <- length(y)
        ind <- sort(y, index.return=TRUE)$ix
        ind <- ind[(n-id.n+1):n]
        text(x[ind] + xrange/50, y[ind], ind, col=ifelse(xmiss[ind],"red","black"))
    }
}

.myddplotNA <- function(md, rd, xmiss, cutoff, id.n) {
    ##  Distance-Distance Plot:
    ##  Plot the vector y=rd (robust distances) against
    ##  x=md (mahalanobis distances). Identify by a label the id.n
    ##  observations with largest rd. If id.n is not supplied, calculate
    ##  it as the number of observations larger than cutoff. Use cutoff
    ##  to draw a horisontal and a vertical line. Draw also a dotted line
    ##  with a slope 1.
    n <- length(md)
    if(missing(id.n))
        id.n <- length(which(rd>cutoff))
    xlab <- "Mahalanobis distance"
    ylab <- "Robust distance"
    plot(md, rd, xlab=xlab, ylab=ylab, type="n")
    points(md, rd, type="p", col=ifelse(xmiss, "red", "black"))
   .labelNA(md, rd, xmiss, id.n)
    abline(0, 1, lty=2)
    abline(v=cutoff)
    abline(h=cutoff)

    title(main="Distance-Distance Plot")
}

.mydistplotNA <- function(x, xmiss, cutoff, classic = FALSE, id.n, ylim=NULL) {
    ## VT::10.11.2007 - change "Squared Robust distance" to "Robust distance"
    ## VT::10.11.2007 - Add parameter ylim to make possible robust and
    ##  classical plot to use the same y-scale

    ##  Index Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the observation indexes. Identify by a label the id.n
    ##  observations with largest value of x. If id.n is not supplied,
    ##  calculate it as the number of observations larger than cutoff.
    ##  Use cutoff to draw a horisontal line.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes

    n <- length(x)
    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    ylab <- paste(if(classic) "Mahalanobis" else "Robust", "distance")
    plot(x, ylab=ylab, xlab="Index", type="n", ylim=ylim)
    points(x, type="p", col=ifelse(xmiss,"red","black"))
    .labelNA(1:n, x, xmiss, id.n)
    abline(h=cutoff)

    title(main="Distance Plot")
}

.qqplotNA <- function(x, p, xmiss, cutoff, classic=FALSE, id.n) {
    ##  Chisquare QQ-Plot:
    ##  Plot the vector x (robust or mahalanobis distances) against
    ##  the square root of the quantiles of the chi-squared distribution
    ##  with p degrees of freedom.
    ##  Identify by a label the id.n observations with largest value of x.
    ##  If id.n is not supplied, calculate it as the number of observations
    ##  larger than cutoff.
    ##  Use classic=FALSE/TRUE to choose the label of the vertical axes


    ##  parameters and preconditions

    n <- length(x)

    if(missing(cutoff))
        cutoff <- sqrt(qchisq(0.975, p))

    if(missing(id.n))
        id.n <- length(which(x>cutoff))

    qq <- sqrt(qchisq(((1:n)-1/3)/(n+1/3), p))

    x <- sort(x, index.return=TRUE)
    ix <- x$ix
    x <- x$x

    if(classic)
        ylab="Mahalanobis distance"
    else
        ylab="Robust distance"

    ## xlab <- "Square root of the quantiles of the chi-squared distribution"
    xlab <- eval(substitute(expression(paste("Sqrt of the quantiles of the ", chi^2, " distribution"))))
    plot(qq, x, xlab=xlab, ylab=ylab, type="n")
    points(qq, x, type="p", col=ifelse(xmiss[ix], "red", "black"))

    if(id.n > 0) {
        ind <- (n-id.n+1):n
        xrange <- par("usr")
        xrange <- xrange[2] - xrange[1]
        text(qq[ind] + xrange/50, x[ind], ix[ind], col=ifelse(xmiss[ix], "red", "black"))
    }
    abline(0, 1, lty=2)
    title(main=eval(substitute(expression(paste(chi^2, " QQ-Plot")))))
}
