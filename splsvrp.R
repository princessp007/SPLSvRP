splsvrp <- function(x, y, K, select="simpls",scale.x=FALSE, scale.y=FALSE, trace=FALSE, fit="simpls", A1, B1, k1, d1)
{
    # always required to input: x, y, K, A1, B1, k1, d1
    # x: matrix
    # y: matrix or vector
    # K: number of latent components
    # select: "simpls" (update X) or "pls2" (update Y)
    # fit: "simpls"
    # scale.x: scale x
    # scale.y: scale y
    # A1: the number of groups of projections
    # B1: the number of projections in a group
    # k1: the sparsity level of the direction vector
    # d1: dimension of the random projections
   

    # initialization

    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    ip <- c(1:p)
    y <- as.matrix(y)
    q <- ncol(y)
    one <- matrix(1,1,n)

    # center & scale y & x

    mu <- one %*% y / n
    y <- scale( y, drop(mu), FALSE )
    meanx <- drop( one %*% x ) / n
    x <- scale( x, meanx, FALSE )

    if ( scale.x )
    {
        normx <- sqrt( drop( one %*% (x^2) ) / (n-1) )
        if ( any( normx < .Machine$double.eps ) )
        { stop("Some of the columns of the predictor matrix have zero variance.") }
        x <- scale( x, FALSE, normx )
    } else { normx <- rep( 1, p ) }

    if ( scale.y )
    {
        normy <- sqrt( drop( one %*% (y^2) ) / (n-1) )
        if ( any( normy < .Machine$double.eps ) )
        { stop("Some of the columns of the response matrix have zero variance.") }
        y <- scale( y, FALSE, normy )
    } else { normy <- rep( 1, q ) }

    # initilize objects

    betahat <- matrix( 0, p, q )
    betamat <- list()
    x1 <- x
    y1 <- y

    # main iteration

    if ( is.null(colnames(x)) )
    { xnames <- c(1:p) } else { xnames <- colnames(x) }

    new2As <- list()
    if ( trace )
    { cat("The variables that join the set of selected variables at each step:\n") }

    for (k in 1:K)
    {

        # fit direction vector
        set.seed(123)
        what <- SPCAvRP(data = x1, cov = FALSE, l = k1 , d = d1, A = A1, B = B1, center_data = FALSE)$vector
        
        # construct A

        A <- unique( ip[ what!=0 | betahat[,1]!=0 ] )
        new2A <- ip[ what!=0 & betahat[,1]==0 ]

        # fit pls with predictors in A

        xA <- x[,A,drop=FALSE]
        plsfit <- pls::plsr( y~xA, ncomp=min(k,length(A)),
                            method=fit, scale=FALSE )

        # update

        betahat <- matrix( 0, p, q )
        betahat[A,] <- matrix( coef(plsfit), length(A), q )
        betamat[[k]] <- betahat # for cv.spls
        pj <- plsfit$projection

        if ( select=="pls2" )
        {
            y1 <- y - x %*% betahat
        }
        if ( select=="simpls" )
        {
            pw <- pj %*% solve( t(pj) %*% pj ) %*% t(pj)
            x1 <- x
            x1[,A] <- x[,A,drop=FALSE] - x[,A,drop=FALSE] %*% pw
        }

        # print out variables that join the active set

        new2As[[k]] <- new2A

        if ( trace )
        {
            if ( length(new2A)<=10 )
            {
                cat( paste("- ",k,"th step (K=",k,"):\n",sep="") )
                cat( xnames[new2A] )
                cat( "\n" )
            } else
            {
                cat( paste("- ",k,"th step (K=",k,"):\n",sep="") )
                nlines <- ceiling(length(new2A)/10)
                for ( i in 0:(nlines-2) )
                {
                    cat( xnames[new2A[(10*i+1):(10*(i+1))]] )
                    cat( "\n" )
                }
                cat( xnames[new2A[(10*(nlines-1)+1):length(new2A)]] )
                cat( "\n" )
            }
        }
    }

    # return objects

    if(!is.null(colnames(x))) {
        rownames(betahat) <- colnames(x)
    }
    if(q > 1 & !is.null(colnames(y))) {
        colnames(betahat) <- colnames(y)
    }

    object <- list(x=x, y=y, betahat=betahat, A=A, betamat=betamat, new2As=new2As,
                    mu=mu, meanx=meanx, normx=normx, normy=normy, K=K, k1=k1, select=select, fit=fit,
                    projection=pj)
    class(object) <- "spls"
    object
}
