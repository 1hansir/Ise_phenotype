sum.I <- function(yy,FUN,Yi,Vi=NULL)
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
  if (!is.null(Vi)) {
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

expit <- function(x){
  1/(1+exp(-x))
}

dexpit <- function(x){
  1/(1+exp(-x))/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}

loglik2 <- function(beta,h0C){
  sum(data$delta*log(expit(h0C+Z%*%beta)) + (1-data$delta)*log(1-expit(h0C+Z%*%beta)),na.rm=T)
}

ddnorm = function(x)
{
  -x*exp(-x^2/2)/sqrt(2*pi)
}

dddnorm = function(x)
{
  (x^2-1)*exp(-x^2/2)/sqrt(2*pi)
}

expexp = function(x)
{
  1-exp(-exp(x))
}


dexpexp = function(x)
{
  exp(x-exp(x))
}