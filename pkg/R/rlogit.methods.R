rpar <- function(x, par, ...){
  x$rpar[[par]]
}
  
print.est.stat <- function(x, ...){
  R <- x$nb.draws
  et <- x$elaps.time[3]
  i <- x$nb.iter
  halton <- x$halton
  eps <- as.numeric(x$eps)
  if (x$type != "simple") cat(paste("Simulated maximum likelihood with",R,"draws\n"))
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  tstr <- paste(h,"h:",m,"m:",s,"s",sep="")
#  cat(paste(i,"iterations,",round(et,0),"seconds\n"))
  cat(paste(i,"iterations,",tstr,"\n"))
  if (!is.null(halton)) cat("Halton's sequences used\n")
  cat(paste("g'(-H)^-1g =",sprintf("%5.3G",eps),"\n"))
}

print.rpar <- function(x, digits = max(3, getOption("digits") - 2), width = getOption("width"), ...){
  dist <- switch(x$dist,
                 "n"="normal",
                 "ln"="log-normal",
                 "cn"="censored normal",
                 "t"="triangular",
                 "u"="uniform"
                 )
  npar1 <- switch(x$dist,
                  "n"="mean",
                  "ln"="meanlog",
                  "cn"="mean",
                  "t"="center",
                  "u"="center"
                 )

  npar2 <- switch(x$dist,
                  "n"="sd",
                  "ln"="sdlog",
                  "cn"="sd",
                  "t"="span",
                  "u"="span"
                 )
  par1 <- x$mean
  par2 <- x$sigma
  cat(paste(dist," distribution with parameters ",round(par1,3)," (",npar1,")"," and ",round(par2,3)," (",npar2,")","\n",sep=""))
}

summary.rpar <- function(object, ...){
  rg <- rg.rpar(object)
  Q1 <- qrlogit(object)(0.25)
  M <-  qrlogit(object)(0.5)
  Q3 <- qrlogit(object)(0.75)
  m <- mean(object)
  r <- c('Min.'=rg[1],'1st Qu.'=Q1,'Median'=M,'Mean'=m,'3rd Qu.'=Q3,'Max.'=rg[2])
  r
}

summary.rlogit <- function(object, ...){
  std.err <- sqrt(diag(vcov(object)))
  b <- coefficients(object)
  z <- b/std.err
  p <- 2*(1-pnorm(abs(z)))
  CoefTable <- cbind(b,std.err,z,p)
  colnames(CoefTable) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")
  CoefTable[attr(object$coefficients,"fixed"),2:4] <- NA
  object$CoefTable <- CoefTable
  if (object$type != "simple"){
    rpar <- object$rpar
    object$summary.rpar <- t(sapply(rpar,summary))
  }
  class(object)="summary.rlogit"
  object
}  

print.summary.rlogit <- function(x, digits = max(3, getOption("digits") - 2),
                              width = getOption("width"), ...){
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  print(x$est.stat)
  cat("\nCoefficients :\n")
  print(x$CoefTable)
  cat(paste("\nlog Likelihood :",round(x$logLik,2),"\n"))
  if (x$type != "simple"){
    cat("\nrandom coefficients\n\n")
    print(x$summary.rpar)
  }
}


cor.rlogit <- function(x){
  cor.rlogit <- cov.rlogit(x)
  K <- nrow(cor.rlogit)
  sd.rlogit <- sqrt(diag(cor.rlogit))
  for (i in 1:K){
    for (j in 1:K){
      cor.rlogit[i,j] <- cor.rlogit[i,j]/sd.rlogit[i]/sd.rlogit[j]
    }
  }
  cor.rlogit
}
  
sd.rlogit <- function(x){
  sqrt(diag(cov.rlogit(x)))
}

cov.rlogit <- function(x){
  if (class(x)[1] != "rlogit") stop("cov.rlogit relevant only for rlogit objects")
  if (x$type != "cor") stop("cov.rlogit relevant only for random models with correlation")
  re <- names(x$call$rpar)[-1]  
  K <- length(re)
  Ktot <- length(coef(x))
  coef.cor <- coef(x)[(Ktot-.5*K*(K+1)+1):Ktot]
  x <- tcrossprod(makeC(coef.cor))
  dimnames(x) <- list(re,re)
  x
}

plot.rpar <- function(x, type = c("density","probability"), norm = TRUE, ...){
  type <- match.arg(type)
  if (type == "density") f <- drlogit
  if (type == "probability") f <- prlogit
  marg <- .05
  law <- x$dist
  rg <- rg.rpar(x)
  low <- rg[1]
  np <- x$name
  neg.values <- ifelse(low < 0,TRUE,FALSE)
  up <- rg[2]
  if (!is.finite(low)) low <- qrlogit(x,norm=norm)(0.005)
  if (!is.finite(up)) up <- qrlogit(x,norm=norm)(0.995)
  ptstot <- gnrpoints(low,up,1000)
  ytot <- do.call(f,list(x=x,norm=norm))(ptstot)
  ymax <- max(ytot)*(1+marg)
  plot(ptstot,ytot,type="n",ann=F,xaxs="i",yaxs="i",las=1,ylim=c(0,ymax),xlim=c(low-marg*(up-low),up+marg*(up-low)))
  ma <- paste("Distribution of",np)
  if (neg.values){
    pourc0 <- prlogit(x)(0)
    print(pourc0)
    ma <- paste(ma,":",round(pourc0*100,0),"% of 0")
    if (type == "density"){
      if (low<0){
        ptsneg <- gnrpoints(low,0,10)
        yneg <- do.call(f,list(x=x,norm=norm))(ptsneg)
        print(c(low,ptsneg))
        polygon(c(low,ptsneg,0),c(0,yneg,0),col="lightblue",border=NA)
      }
    }
    else{
      segments(low-marg*(up-low),pourc0,0,pourc0,lty="dotted")
      segments(0,0,0,pourc0,lty="dotted")
    }
  }
  lines(ptstot,ytot)
  if (law=="u" && type == "density"){
    segments(up,0,up,drlogit(x)(up))
    segments(low,0,low,drlogit(x)(low))
  }
  title(main=ma)
}

plot.rlogit <- function(x, ...){
  rpar <- x$rpar
  K <- length(rpar)
  nrow <- 1+(K>2)+(K>6)
  ncol <- 1+(K>1)+(K>4)
  opar <- par(mfrow=c(nrow,ncol))
  
  for (i in names(rpar)){
    plot(rpar(x,i), ...)
  }
  par(opar)
}
