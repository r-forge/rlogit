m2norm <- function(m,dist,norm){
  switch(dist,
         "n"=m/norm,
         "ln"=m-log(norm),
         "t"=m/norm,
         "cn"=m/norm,
         "u"=m/norm
         )
}

s2norm <- function(s,dist,norm){
  switch(dist,
         "n"=s/norm,
         "ln"=s,
         "t"=s/norm,
         "cn"=s/norm,
         "u"=s/norm
         )
}

# mean methods for rpar and rlogit objects
  
mean.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm
  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn)
  }
  switch(dist,
         "n"=m,
         "ln"=exp(m+0.5*s^2),
         "u"=m,
         "t"=m,
         "cn"=s*dnorm(-m/s)+m*(1-pnorm(-m/s))
         )
}

mean.rlogit <- function(x, par = 1, norm = TRUE, ...){
  x <- x$rpar[[par]]
  mean(x, norm = norm, ...)
}

# median methods for rpar and rlogit objects

med <- function(x, ...){
  UseMethod("med")
}

med.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm
  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn) 
  }
  switch(dist,
         "n"=m,
         "ln"=exp(m),
         "u"=m,
         "t"=m,
         "cn"=0
         )
}

med.rlogit <- function(x, par = 1, norm = TRUE, ...){
  x <- x$rpar[[par]]
  med(x, norm = norm, ...)
}

# stdev methods for rpar and rlogit objects

stdev <- function(x, ...){
  UseMethod("stdev")
}

stdev.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm
  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn)
  }
  switch(dist,
         "n"=s,
         "ln"=sqrt(exp(s^2)-1)*exp(m+0.5*s^2),
         "u"=s^2/3,
         "t"=s,
         "cn"=sqrt( s^2*(1-pnorm(-m/s))+m*(s*dnorm(-m/s)+m*(1-pnorm(-m/s)))-(s*dnorm(-m/s)+m*(1-pnorm(-m/s)))^2)
         )
}

stdev.rlogit <- function(x, par = 1, norm = TRUE, ...){
  x <- x$rpar[[par]]
  stdev(x, norm = norm, ...)
}

# qrlogit methods for rpar and rlogit objects

qrlogit <- function(x, ...){
  UseMethod("qrlogit")
}

qrlogit.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm
  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn)
  }
  switch(dist,
         "n"=function(x=(1:9)/10) qnorm(x,m,s),
         "ln"=function(x=(1:9)/10) qlnorm(x,m,s),
         "u"=function(x=(1:9)/10) qunif(x,m-s,m+s),#( m-s+2*s*x  )*(x > m-s & x < m+s) + 0,
         "t"=function(x=(1:9)/10) (m-s+sqrt(2*s^2*x))*(x<=0.5)+(m+s-sqrt(2*s^2*(1-x)))*(x>0.5),
         "cn"=function(x=(1:9)/10) qnorm(x,m,s)
         )
}

qrlogit.rlogit <- function(x, par = 1, y = NULL, norm = TRUE, ...){
  x <- x$rpar[[par]]
  if (is.null(y)){
    qrlogit(x, norm = norm, ...)
  }
  else{
    qrlogit(x, norm = norm, ...)(y)
  }
}

# prlogit methods for rpar and rlogit objects

prlogit <- function(x, ...){
  UseMethod("prlogit")
}

prlogit.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm
  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn)
  }
  switch(dist,
         "n"=function(x) pnorm(x,m,s),
         "ln"=function(x) plnorm(x,m,s),
         "u"=function(x) punif(x,m-s,m+s),#(x-m+s)/(2*s)*(x > m-s & x < m+s) + 0,
         "t"=function(x) (x >= (m-s) & x < m)*(x-m+s)^2/(2*s^2)+(x>=m & x <= (m+s))*(1-(m+s-x)^2/(2*s^2))+(x>(m+s))*1+0,
         "cn"=function(x) pnorm(x,m,s)
         )
}

prlogit.rlogit <- function(x, par = 1, y = NULL, norm = TRUE, ...){
  x <- x$rpar[[par]]
  if (is.null(y)){
    prlogit(x, norm = norm, ...)
  }
  else{
    prlogit(x, norm = norm, ...)(y)
  }
}

# drlogit methods for rlogit and rpar objects

drlogit <- function(x, ...){
  UseMethod("drlogit")
}

drlogit.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm
  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn)
  }
  switch(dist,
         "n"=function(x) dnorm(x,m,s),
         "ln"=function(x) dlnorm(x,m,s),
         "u"=function(x) (1/s+x*0)*(x >= m-s & x <= m+s) + 0,
         "t"=function(x) (x >= (m-s) & x < m)*(x-m+s)/s^2+(x>=m & x <= (m+s))*(s+m-x)/s^2+0,
         "cn"=function(x) dnorm(x,m,s),

         )
}

drlogit.rlogit <- function(x, par = 1, y = NULL, norm = TRUE, ...){
  x <- x$rpar[[par]]
  if (is.null(y)){
    drlogit(x, norm = norm, ...)
  }
  else{
    drlogit(x, norm = norm, ...)(y)
  }
}

rg.rpar <- function(x, norm = TRUE, ...){
  dist <- x$dist
  m <- x$mean
  s <- abs(x$sigma)
  vn <- x$norm

  if (norm && !is.null(vn)){
    s <- s2norm(s,dist,vn)
    m <- m2norm(m,dist,vn)
  }
  switch(dist,
         "n"=c(-Inf,+Inf),
         "ln"=c(0,+Inf),
         "u"=c(m-s,m+s),
         "t"=c(m-s,m+s),
         "cn"=c(0,+Inf)
         )
}
