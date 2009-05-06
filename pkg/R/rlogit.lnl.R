compute.P <- function(param, Xa, Xc, Varc, Vara, eta, rpar, correlation){
  betac <- param[Varc]
  mua <- param[Vara]
  K <- length(Varc)+length(Vara)
  siga <- param[-(1:K)]
  names(mua) <- names(siga) <- colnames(Xa[[1]])
  x <- makebeta(mua, siga, rpar, eta, correlation)
  betaa <- x$betaa
  A <- lapply(Xc, function(x) as.vector(crossprod(t(as.matrix(x)), betac)))
  B <- lapply(Xa, function(x) tcrossprod(as.matrix(x), betaa))
  AB <- mapply(function(x, y) exp(x+y), A, B, SIMPLIFY=FALSE)
  S <- suml(AB)
  P <- lapply(AB, function(x) x/S)
  P
}

compute.Pch <- function(P, y, id){
  Pch <- suml(mapply("*", P, y, SIMPLIFY=FALSE))
  if (!is.null(id)){
    Pch <- apply(Pch, 2, tapply, id, prod)
  }
  Pch
}

compute.lnl <- function(param, y, Xa, Xc, Varc, Vara, eta, weights, id, rpar, correlation){
  if (is.null(weights)) weights <- 1
  P <- compute.P(param, Xa, Xc, Varc, Vara, eta, rpar, correlation)
  Pch <- compute.Pch(P, y, id)
  pm <- apply(Pch, 1, mean)
  sum(weights*log(pm))
}

compute.g <- function(param, y, Xa, Xc, Varc, Vara, eta, weights, id, rpar, correlation){
  xac <- suml(mapply("*", Xa, y, SIMPLIFY=FALSE))
  xcc <- suml(mapply("*", Xc, y, SIMPLIFY=FALSE))
  Kc <- length(Varc)
  Ka <- length(Vara)
  K <- Kc+Ka
  mua <- param[Vara]
  siga <- param[-(1:K)]
  names(mua) <- names(siga) <- colnames(Xa[[1]])
  b <- makebeta(mua, siga, rpar, eta, correlation)
  P <- compute.P(param, Xa, Xc, Varc, Vara, eta, rpar, correlation)
  Pch <- compute.Pch(P, y, id)
  if (!is.null(id)) Pch <- Pch[as.character(id), ]
  pm <- apply(Pch, 1, mean)
  R <- nrow(eta)
  if (correlation){
    names.cor <- c()
    for (i in 1:Ka){
      names.cor <- c(names.cor, paste(names(rpar)[i], names(rpar)[i:Ka], sep=":"))
    }
    vecX <- c()
    for (i in 1:Ka){
      vecX <- c(vecX, i:Ka)
    }
    Xas <- lapply(Xa,  function(x) x[, vecX])
    xac <- suml(mapply("*", Xa, y, SIMPLIFY=FALSE))
    xacs <- suml(mapply("*", Xas, y, SIMPLIFY=FALSE))
    colnames(xacs) <- names(param)[-(1:(Ka+Kc))]
  }
  else{
    xacs <- xac
    Xas <- Xa
  }
  PCP <- lapply(P, function(x) Pch*x)
  PCPs <- lapply(PCP, function(x) apply(x, 1, sum))
  grad.cst <- xcc-suml(mapply("*", Xc, PCPs, SIMPLIFY=FALSE))/(R*pm)
  grad.mu <- (tcrossprod(Pch, t(b$betaa.mu))*xac -
    suml(mapply(function(x, y) x*tcrossprod(y, t(b$betaa.mu)), Xa, PCP, SIMPLIFY=FALSE)))/(R*pm)
  grad.sd <- (tcrossprod(Pch, t(b$betaa.sigma))*xacs -
    suml(mapply(function(x, y) x*tcrossprod(y, t(b$betaa.sigma)), Xas, PCP, SIMPLIFY=FALSE)))/(R*pm)
  if (is.null(weights)){
    gradi <- cbind(grad.cst, grad.mu, grad.sd)
  }
  else{
    gradi <- cbind(grad.cst, grad.mu, grad.sd)*weights
  }
  colnames(gradi) <- names(param)
  gradi
}

compute.tot <- function(param, y, Xa, Xc, Varc, Vara, eta, weights,
                        id, rpar, correlation, gradient=FALSE){

  if (is.null(weights)) weights <- rep(1, length(y[[1]]))
  P <- compute.P(param, Xa, Xc, Varc, Vara, eta, rpar, correlation)
  Pch <- compute.Pch(P, y, id)
  pm <- apply(Pch, 1, mean)
  lnL <- sum(weights[!duplicated(id)]*log(pm))
  if (gradient){
    xac <- suml(mapply("*", Xa, y, SIMPLIFY=FALSE))
    xcc <- suml(mapply("*", Xc, y, SIMPLIFY=FALSE))
    Kc <- length(Varc)
    Ka <- length(Vara)
    K <- Kc+Ka
    mua <- param[Vara]
    siga <- param[-(1:K)]
    names(mua) <- names(siga) <- colnames(Xa[[1]])
    b <- makebeta(mua, siga, rpar, eta, correlation)
#    P <- compute.P(param, Xa, Xc, Varc, Vara, eta, rpar, correlation)
#    Pch <- compute.Pch(P, y, id)
    if (!is.null(id)) Pch <- Pch[as.character(id), ]
    pm <- apply(Pch, 1, mean)
    R <- nrow(eta)
    if (correlation){
      names.cor <- c()
      for (i in 1:Ka){
        names.cor <- c(names.cor, paste(names(rpar)[i], names(rpar)[i:Ka], sep=":"))
      }
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX, i:Ka)
      }
      Xas <- lapply(Xa,  function(x) x[, vecX])
      xac <- suml(mapply("*", Xa, y, SIMPLIFY=FALSE))
      xacs <- suml(mapply("*", Xas, y, SIMPLIFY=FALSE))
      colnames(xacs) <- names(param)[-(1:(Ka+Kc))]
    }
    else{
      xacs <- xac
      Xas <- Xa
    }
    PCP <- lapply(P, function(x) Pch*x)
    PCPs <- lapply(PCP, function(x) apply(x, 1, sum))
    grad.cst <- xcc-suml(mapply("*", Xc, PCPs, SIMPLIFY=FALSE))/(R*pm)
    grad.mu <- (tcrossprod(Pch, t(b$betaa.mu))*xac -
                suml(mapply(function(x, y) x*tcrossprod(y, t(b$betaa.mu)),
                            Xa, PCP, SIMPLIFY=FALSE)))/(R*pm)
    grad.sd <- (tcrossprod(Pch, t(b$betaa.sigma))*xacs -
                suml(mapply(function(x, y) x*tcrossprod(y, t(b$betaa.sigma)),
                            Xas, PCP, SIMPLIFY=FALSE)))/(R*pm)
    if (is.null(weights)){
      gradi <- cbind(grad.cst, grad.mu, grad.sd)
    }
    else{
      gradi <- cbind(grad.cst, grad.mu, grad.sd)*weights
    }
    colnames(gradi) <- names(param)
    
    
    attr(lnL, "gradi") <- gradi
  }
  lnL
}
    
  
