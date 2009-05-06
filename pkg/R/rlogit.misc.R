makebeta <- function(mua, siga, rpar, eta, correlation){
  nr <- names(rpar)
  censored <- nr[rpar=="cn"]
  lognormal <- nr[rpar=="ln"]
  truncated <- nr[rpar=="tn"]
  normal  <- nr[rpar=="n"]
  uniform  <- nr[rpar=="u"]
  triangular  <- nr[rpar=="t"]

  Ka <- ncol(eta)
  R <- nrow(eta)
  
  betaa <- matrix(NA,R,Ka)
  betaa.mu <- betaa.sigma <- betaa

  colnames(betaa) <- colnames(betaa.mu) <- colnames(betaa.sigma) <- colnames(eta) <- names(mua)

  if (correlation){
    colnames(eta) <- NULL
    CC <- makeC(siga)
    sigeta <- tcrossprod(eta,CC)
    colnames(sigeta) <- nr
    mymua <- mua
    mysiga <- sigeta
    betaa <- t(mymua+t(mysiga))
    betaa.mu <- matrix(1,R,Ka)
    betaa.sigma <- eta[,rep(1:Ka,Ka:1)]
    for (i in 1:Ka){
      sigi <- i+cumsum(c(0,(Ka-1):1))[1:i]
      if (rpar[i]=="cn"){
        betaa[,i] <- pmax(betaa[,i],0)
        betaa.mu[,i] <- (betaa[,i] > 0)*1 +0
        betaa.sigma[,sigi] <- ( (betaa[,i] > 0)*1 + 0)*betaa.sigma[,sigi]
      }
      if (rpar[i]=="ln"){
        betaa[,i] <- exp(betaa[,i])
        betaa.mu[,i] <- betaa[,i]
        betaa.sigma[,sigi] <- betaa[,i]*betaa.sigma[,sigi]
      }
    }
  }
  else{
    if(length(censored)>0){
      sel <- censored
      betaa[,sel] <- pmax(t(mua[sel]+siga[sel]*t(eta[,sel,drop=F])),0)
      betaa.mu[,sel] <- as.numeric(betaa[,sel]>0)
      betaa.sigma[,sel] <- betaa.mu[,sel]*eta[,sel]
    }
    if(length(lognormal)>0){
      sel <- lognormal
      betaa[,sel] <- exp(t(mua[sel]+siga[sel]*t(eta[,sel,drop=F])))
      betaa.mu[,sel] <- betaa[,sel,drop=F]
      betaa.sigma[,sel] <- betaa.mu[,sel]*eta[,sel]
    }
    if(length(normal)>0){
      sel <- normal
      betaa[,sel] <- t(mua[sel]+siga[sel]*t(eta[,sel,drop=F]))
      betaa.mu[,sel] <- 1
      betaa.sigma[,sel] <- eta[,sel,drop=F]
    }

    if(length(uniform)>0){
      sel <- uniform
      etauni <- pnorm(eta[,sel,drop=F])
      betaa[,sel] <- t(mua[sel]-siga[sel]+2*t(etauni)*siga[sel])
      betaa.mu[,sel] <- 1
      betaa.sigma[,sel] <- 2*etauni-1
    }

    if(length(triangular)>0){
      sel <- triangular
      eta05 <- eta[,sel,drop=F]<0.5
      betaa.mu[,sel] <- 1
      betaa.sigma[,sel] <- eta05*(sqrt(2*pnorm(eta[,sel,drop=F]))-1)+
        !eta05*(1-sqrt(2*(1-pnorm(eta[,sel,drop=F]))))
      
      betaa[,sel] <- t(mua[sel]+siga[sel]*t(betaa.sigma[,sel]))
    }
  }
  list(betaa=betaa,betaa.mu=betaa.mu,betaa.sigma=betaa.sigma)
}


parnamed <- function(x, names, defdef, argum){
  if (is.null(x)) stop("ces nul\n")
  def <- FALSE
  K <- length(names)
  mod <- mode(defdef)
  if (is.null(names(x))){
    if ((length(x) != length(names)) & (length(x) > 1)){
      
      stop(paste("the length of the unnamed",
                 argum,
                 "argument should be",
                 length(names)))
    }
    if (length(x) == length(names)) names(x) <- names
    if (length(x)==1 || is.null(names(x))){
      x <- rep(x,K)
      names(x) <- names
    }
  }
  else{
    defpos <- names(x)==""
    namedcoef <- names(x)[!defpos]
    if (sum(defpos) > 1) stop("only one default value allowed")
    if (sum(defpos) == 0) def <- defdef
    else{
      unknown.names <- namedcoef[!(namedcoef %in% names)]
      if (length(unknown.names) > 0){
        unknown.names <- paste(unknown.names,collapse=", ")
        stop(paste("unknow name(s) in the",argum,":",unknown.names))
      }
      def <- x[defpos]
      if (mode(def)!=mod) stop(paste("the default value must be",mod))
    }
    ox <- x
    x <- rep(def,K)
    names(x) <- names
    x[namedcoef] <- ox[namedcoef]
  }
  x
}

gnrpoints <- function(low,up,n=100){
  low+(up-low)*(0:n)/n
}

halton <- function(prime=3,length=100,drop=10){
  halt <- 0
  t <- 0
  while(length(halt)<length+drop){
    t <- t+1
    halt <- c(halt,rep(halt,prime-1)+rep(seq(1,prime-1,1)/prime^t,each=length(halt)))
  }
  halt[(drop+1):(length+drop)]
}

makeC <- function(x){
  K <- (-1+sqrt(1+8*length(x)))/2
  mat <- matrix(0,K,K)
  mat[lower.tri(mat,diag=TRUE)] <- x
  mat
}

  

# Create the matrix of random numbers

make.eta <- function(R,Ka,halton){
  if (!is.null(halton)){
    length.halton <- rep(R,Ka)
    prime <- c(2,3,5,7,11,13,17,19,23)
    drop.halton <- rep(100,Ka)
    if (!is.na(halton) && !is.null(halton$prime)){
      if (length(halton$prime) != Ka){
        stop("wrong number of prime numbers indicated")
      }
      else{
        prime <- halton$prime
      }
      if (!is.na(halton) && !is.null(halton$drop)){
        if (!length(halton$drop) %in% c(1,Ka)) stop("wrong number of drop indicated")
        if (length(halton$drop) == 1){
          drop.halton <- rep(halton$drop,Ka)
        }
        else{
          drop.halton <- halton$drop
        }
      }
    }
    eta <- numeric(0)
    i <- 0
    for (i in 1:Ka){
      eta <- cbind(eta,qnorm(halton(prime[i],R,drop.halton[i])))
    }
  }
  else{
    eta <- matrix(rnorm(R*Ka),ncol=Ka,nrow=R)
  }
  eta
}


# Create the vector of names for the coefficients

make.name <- function(param,correlation,X,Vara,rpar){
  if (!correlation){
    names.param <- c(colnames(X),paste("sd.",colnames(X)[Vara],sep=""))
  }
  else{
    names.param <- c()
    Ka <- length(rpar)
    for (i in 1:Ka){
      names.param <- c(names.param,paste(names(rpar)[i],names(rpar)[i:Ka],sep="."))
    }
    names.param <- c(colnames(X),names.param)
  }
  if (length(param) != length(names.param)){
    stop(paste("the length of start is",
               length(param),
               "and the number of coefficients is",
               length(names.param)))
  }
  names.param
}

# Create the X matrix



make.rpar <- function(rpar,correlation,estimate,norm){
  K <- length(rpar)
  nr <- names(rpar)
  rpar <- lapply(rpar,function(x) list(dist=x))
  if (correlation){
    Ktot <- length(estimate)
    v <- estimate[(Ktot-0.5*K*(K+1)+1):Ktot]
    v <- tcrossprod(makeC(v))
    sv <- sqrt(diag(v))
    names(sv) <- nr
  }      
  for (i in (1:K)){
    m <- estimate[nr[i]]
    if (!correlation){
      s <- estimate[paste("sd.",nr[i],sep="")]
    }
    else{
      s <- sv[i]
    }
    names(m) <- names(s) <- NULL
    rpar[[i]]$mean <- m
    rpar[[i]]$sigma <- s
    rpar[[i]]$name <- nr[[i]]
    if (!is.null(norm)){
      vn <- estimate[norm]
      names(vn) <- NULL
      rpar[[i]]$norm <- vn
    }
  }
  lapply(rpar,function(x){attr(x,"class")="rpar";x})
}
  
rpar.factor <- function(rpar,data){
  for (i in 1:length(rpar)){
    n <- names(rpar)[i]
    clvar <- class(data[[n]])
    if (clvar=="factor"){
      lvar <- levels(data[[n]])[2]
      names(rpar)[i] <- paste(n,lvar,sep="")
    }
  }
  rpar
}

suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))){
    d <- dim(x[[1]])
    s <- matrix(0,d[1],d[2])
    for (i in 1:n){
      s <- s+x[[i]]
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n){
      s <- s+x[[i]]
    }
  }
  s
}
