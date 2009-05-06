mynlm <- function(f, param,  ...){
  names.coef <- names(param)
                                        #  if is.list(maxit = 500, tol = 1E-05, print.level = 2, mult = FALSE, ...){
  mult.sup <- 4
  mult.inf <- 1/64
  maxit = 500
  tol = 3E-08
  print.level = 2
  mult = TRUE
#  mult = FALSE
  
  lim <- match.call(expand.dots=F)
  tot <- match.call()
  func <- tot[[2]]
  m <- match(c('param','id','weights'),names(tot), 0L)
  f <- tot[c(1L,m)]
  f[[1]] <- func
  f$gradient <- TRUE

  fixed <- attr(param,"fixed")
  if (is.null(fixed)) fixed <- rep(F,length(param))
  
  
  chi2 <- 1E+10
  i <- 0
  xini <- eval(f, parent.frame())
  
  while(chi2 > tol){
    f$gradient <- TRUE
    i <- i+1
    f$param <- param
    gradi <- attr(xini,"gradi")
    g <- if(is.matrix(gradi)) apply(gradi,2,sum) else sum(gradi)
    H <- crossprod(gradi)
    incr <- as.vector(crossprod(solve(H[!fixed,!fixed]),g[!fixed]))
    lambda <- 1
    param[!fixed] <- param[!fixed]+incr
    f$param <- param
    x <- eval(f, parent.frame())
    xx <- c(as.numeric(xini),as.numeric(x))
    if (!mult){
      param[!fixed] <- param[!fixed]+incr
    }
    else{
#      f$gradient=FALSE
      if (x < xini){
        xx <- c(as.numeric(xini),as.numeric(x))
        while(xx[length(xx)] < xx[1] && lambda >= mult.inf){
          lambda <- 1/2*lambda
          f$param <- param+lambda*incr
          x <- eval(f, parent.frame())
          xx <- c(xx,x)
          xi <- x[length(x)]
        }
        param <- param+lambda*incr
      }
      else{
        xi <- x
      }
    }
    chi2 <- crossprod(incr,g[!fixed])
    if (print.level>0){
      if (mult){
        chaine <- paste("iteration ",i,", lambda = ",lambda,", lnL = ",round(xi,8),", chi2 = ",round(chi2,8),"\n",sep="")
      }
      else{
        chaine <- paste("iteration ",i,", lnL = ",round(x,2),", chi2 = ",round(chi2,8),"\n",sep="")
      }        
      cat(chaine)
      cat("\n")
    }
    if (print.level>1){
      resdet <- rbind(param=param,gradient=g)
      print(round(resdet,3))
      cat("--------------------------------------------\n")
    }
    if (i>maxit){stop("maximum number of iterations reached\n")}
    xini <- x
  }
  attr(param,"nb.iter") <- i
  attr(param,"eps") <- chi2
  names(param) <- names.coef
  param
}
