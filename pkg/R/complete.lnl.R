
complete.lnl <- function(param,y,Xa,Xc,Varc,Vara,eta,weights=NULL,id,rpar,correlation,gradient=TRUE){
  if (is.null(weights)) weights <- 1
  betac <- param[Varc]
  mua <- param[Vara]
  Kc <- length(Varc)
  Ka <- length(Vara)
  K <- Kc+Ka
  siga <- param[-(1:K)]
  names(mua) <- names(siga) <- colnames(Xa[[1]])
  x <- makebeta(mua,siga,rpar,eta,correlation)
  betaa <- x$betaa
  A <- lapply(Xc,function(x) as.vector(crossprod(t(as.matrix(x)),betac)))
  B <- lapply(Xa,function(x) tcrossprod(as.matrix(x),betaa))
  AB <- mapply(function(x,y) exp(x+y),A,B,SIMPLIFY=FALSE)
  S <- suml(AB)
  P <- lapply(AB,function(x) x/S)
  Pch <- suml(mapply("*",P,y,SIMPLIFY=FALSE))
  if (!is.null(id)){
    Pch <- apply(Pch,2,tapply,id,prod)
  }
  pm <- apply(Pch,1,mean)
  lnl <- sum(weights*log(pm))

  if (gradient){
    xac <- suml(mapply("*",Xa,y,SIMPLIFY=FALSE))
    xcc <- suml(mapply("*",Xc,y,SIMPLIFY=FALSE))
    b <- makebeta(mua,siga,rpar,eta,correlation)
    if (!is.null(id)) Pch <- Pch[as.character(id),]
    pm <- apply(Pch,1,mean)
    # trois lignes précédentes commentées
    R <- nrow(eta)
    if (correlation){
      names.cor <- c()
      for (i in 1:Ka){
        names.cor <- c(names.cor,paste(names(rpar)[i],names(rpar)[i:Ka],sep=":"))
      }
      vecX <- c()
      for (i in 1:Ka){
        vecX <- c(vecX,i:Ka)
      }
      Xas <- lapply(Xa, function(x) x[,vecX])
      xac <- suml(mapply("*",Xa,y,SIMPLIFY=FALSE))
      xacs <- suml(mapply("*",Xas,y,SIMPLIFY=FALSE))
      colnames(xacs) <- names(param)[-(1:(Ka+Kc))]
    }
    else{
      xacs <- xac
      Xas <- Xa
    }
    PCP <- lapply(P,function(x) Pch*x)
    PCPs <- lapply(PCP,function(x) apply(x,1,sum))
    grad.cst <- xcc-suml(mapply("*",Xc,PCPs,SIMPLIFY=FALSE))/(R*pm)
    grad.mu <- (tcrossprod(Pch,t(b$betaa.mu))*xac -
                suml(mapply(function(x,y) x*tcrossprod(y,t(b$betaa.mu)),Xa,PCP,SIMPLIFY=FALSE)))/(R*pm)
    grad.sd <- (tcrossprod(Pch,t(b$betaa.sigma))*xacs -
                suml(mapply(function(x,y) x*tcrossprod(y,t(b$betaa.sigma)),Xas,PCP,SIMPLIFY=FALSE)))/(R*pm)
    if (is.null(weights)){
      gradi <- cbind(grad.cst,grad.mu,grad.sd)
    }
    else{
      gradi <- cbind(grad.cst,grad.mu,grad.sd)*weights
    }
    colnames(gradi) <- names(param)
    gradi
    attr(lnl,"gradient") <- gradi
  }
}
