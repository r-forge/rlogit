rlogit <- function(formula, data, rpar = NULL, correlation = FALSE, weights = NULL, id = NULL,
                   start = NULL, R = 50, eta = NULL, halton = NULL, fixed = FALSE, norm = NULL,
                   ...){
  if (is.null(rpar)) stop("no random coefficients, use mlogit")
  class(formula) <- c("logitform","formula")
  start.time <- proc.time()
  cl <- match.call()

  type <- ifelse(correlation,"cor","diag")
  alt.name <- names(data)[2]
#  formula <- make.formula(formula,alt.name)
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  y <- mf[,3]
  choiceid <- data[[1]]
  if (!is.null(id)) id <- data[[id]]
  alt.name <- names(data)[[2]]
#  X <- make.X(formula,data)
  X <- model.matrix(formula,data)
  df.residual <- nrow(X)-ncol(X)
  n <- length(y)
  K <- ncol(X)
  rpar <- rpar.factor(rpar,data)                                          # change the names of rp in case of factor
  Vara <- sort(match(names(rpar),colnames(X)))
  Varc <- (1:K)[-Vara]
  nmean <- length(c(Varc,Vara))
  nvar <- length(Vara)
  K <- ifelse(correlation,nmean+.5*nvar*(nvar+1),nmean+nvar)
  Ka <- length(Vara)
  if (is.null(eta)) eta <- make.eta(R,Ka,halton)                          # create the random numbers matrix
  colnames(eta) <- colnames(X)[Vara]
  if (is.null(start)){
    start <- coef(mlogit(formula,data))
    cat("starting values computed\n")
    ln <- names(rpar[rpar=="ln"])
    start[ln] <- log(start[ln])
    start <- c(start,rep(.1,K-length(start)))
  }
  names(start) <- make.name(start,correlation,X,Vara,rpar)                # name the coefficients
  fixed <- parnamed(fixed,names(start),F,"fixed")
  start <- structure(start,fixed=fixed)
  if(!is.null(weights)){
    weights <- data[[weights]]/mean(data[[weights]])
#    print(head(choiceid));print(length(choiceid));stop()
    weights <- weights[!duplicated(choiceid)]
  }
  alt <- data[[alt.name]]
  Xl <- split(as.data.frame(X),alt)
  Xl <- lapply(Xl,as.matrix)
  Xal <- lapply(Xl, function(x) x[,Vara,drop=F])
  Xcl <- lapply(Xl, function(x) x[,Varc,drop=F])
  yl <- split(y,alt)
  if (!is.null(id)) idl <- split(id,alt)[[1]] else idl <- NULL

  use.maxlik <- (is.null(cl$method) || cl$method != 'mynlm')
  
  if (use.maxlik){
##     f <- function(param) compute.lnl(param, yl, Xal, Xcl, Varc, Vara, eta,
##                                      weights = weights, id = idl, rpar = rpar,
##                                      correlation = correlation)
    
##     g <- function(param) compute.g(param, yl, Xal, Xcl, Varc, Vara, eta,
##                                    weights = weights, id = idl, rpar = rpar,
##                                    correlation = correlation)

    f <- function(param)
      compute.tot(param, yl, Xal, Xcl, Varc, Vara, eta,
                  weights = weights, id = idl, rpar = rpar,
                  correlation = correlation,gradient = FALSE)
    g <- function(param)
      attr(compute.tot(param, yl, Xal, Xcl, Varc, Vara, eta,
                       weights = weights, id = idl, rpar = rpar,
                       correlation = correlation,gradient = TRUE),"gradi")

    result <- maxLik(f, g, start = start, ...)
    estimate <- result$estimate
    ll <- result$maximum
    gradient <- apply(g(estimate),2,sum)
    hessian <- result$hessian
    nb.iter <- result$iterations
    eps <- gradient%*%solve(-result$hessian)%*%gradient
    
  }
  else{
    za <- function(param, id,  weights, gradient)
      compute.tot(param, yl, Xal, Xcl, Varc, Vara, eta,
                  weights = weights, id = idl, rpar = rpar,
                  correlation = correlation,gradient)
    estimate <- mynlm(za, start, id = idl, weights = weights)
    nb.iter <- attr(estimate,"nb.iter")
    eps <- attr(estimate,"eps")
    attr(estimate, "nb.iter") <- attr(estimate,"eps") <- NULL
    ll <- compute.tot(estimate, yl, Xal, Xcl, Varc, Vara, eta,
                       weights = weights, id = idl, rpar = rpar,
                       correlation = correlation,gradient=TRUE)
    gradi <- attr(ll, 'gradi')
    gradient <- apply(gradi,2,sum)
    attributes(ll) <- NULL
    hessian <- -crossprod(gradi)
    rownames(hessian) <- colnames(hessian) <- names(estimate)

  }
  
  
  rpar <- make.rpar(rpar,correlation,estimate,norm)                      # create the rpar list
  elaps.time <- proc.time() - start.time
  est.stat <- list(elaps.time = elaps.time,                              # create the estimation summary
                   nb.iter = nb.iter,
                   nb.draws = R,
                   halton = halton,
                   eps = eps,
                   type = type
                   )
  class(est.stat) <- "est.stat"
  attr(estimate,"fixed") <- NULL
  rlogit <- list(coefficients = estimate,
              X = X, hessian = hessian, gradient = gradient,
              call = cl, model = mf, df.residual = df.residual,
              logLik = ll,type = type, fixed = fixed,
              rpar = rpar,est.stat = est.stat,eta = eta)
  class(rlogit) <- c("rlogit","mlogit")
  rlogit
}


  
