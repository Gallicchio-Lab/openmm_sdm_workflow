# cat r*/*.out > data
# R CMD BATCH uwham_analysis.R

#.libPaths("/home/emilio/R/x86_64-pc-linux-gnu-library/3.0/")
library("UWHAM")

bias.fcn <- function(ebind, lam, lam1, lam2, alpha, u0, w0){
# This is for the bias ilogistic potential
# (lambda2-lambda1) ln[1+exp(-alpha (u-u0))]/alpha + lambda2 u + w0
    ebias1 <- 0*ebind
    if (alpha > 0) {
        ee <- 1 + exp(-alpha*(ebind-u0))
        ebias1 <- (lam2 - lam1)*log(ee)/alpha
    }
    ebias1 + lam2*ebind + w0
}

npot.fcn <- function(e0,ebind, bet, lam, lam1, lam2, alpha, u0, w0){ 
# This is the negative reduced energy 
# -beta*(U0+bias)
    -bet*(e0 + bias.fcn(ebind,lam, lam1, lam2, alpha, u0, w0))
}

uwham.r <- function(label,logQ,ufactormax,ufactormin=1){
  n <- dim(logQ)[1]
  m <- dim(logQ)[2]
  iniz <- array(0,dim=m) 
  uf <- ufactormax
  while(uf >= ufactormin & uf >= 1){
    mask <- seq(1,n,trunc(uf))
    out <- uwham(label=label.cross[mask], logQ=neg.pot[mask,],init=iniz)
    show(uf)
    iniz <- out$ze
    uf <- uf/2
  }
  out$mask <- mask
  out
}

histw <-
function (x, w, xaxis, xmin, xmax, ymax, bar = TRUE, add = FALSE, 
            col = "black", dens = TRUE) 
{
  nbin <- length(xaxis)
  xbin <- cut(x, breaks = xaxis, include.lowest = T, labels = 1:(nbin -  1))
  y <- tapply(w, xbin, sum)
  y[is.na(y)] <- 0
  y <- y/sum(w)
  if (dens) 
    y <- y/(xaxis[-1] - xaxis[-nbin])
  if (!add) {
    plot.new()
    plot.window(xlim = c(xmin, xmax), ylim = c(0, ymax))
    axis(1, pos = 0)
    axis(2, pos = xmin)
  }
  if (bar == 1) {
    rect(xaxis[-nbin], 0, xaxis[-1], y)
  }
  else {
    xval <- as.vector(rbind(xaxis[-nbin], xaxis[-1]))
    yval <- as.vector(rbind(y, y))
    lines(c(min(xmin, xaxis[1]), xval, max(xmax, xaxis[length(xaxis)])), 
          c(0, yval, 0), lty = "11", lwd = 2, col = col)
  }
  invisible()
  list(y = y, breaks = xaxis)
}



data.t <- read.table("repl.cycle.potE.temp.lambda.ebind.lambda1.lambda2.alpha.u0.w0.dat")

potEs = data.t$V3
temps = data.t$V4
lambdas = data.t$V5
ebinds = data.t$V6
lambda1s = data.t$V7
lambda2s = data.t$V8
alphas = data.t$V9
u0s = data.t$V10
w0s = data.t$V11

#parameters for ilogistic expected
lam     <-c( %s )
lambda1 <-c( %s )
lambda2 <-c( %s )
alpha   <-c( %s )
u0      <-c( %s )
w0coeff <-c( %s )


tempt <- c( %f )
bet <- 1.0/(0.001986209*tempt)
mtempt <- length(bet)
mlam <- length(lam)
m <- mlam*mtempt
N <- length(data.t$V1)

#extract U0 values as U-bias
#this is relevant only if the states are at different temperatures
e0 <- potEs
for (i in 1:N) {
    e0[i] <- e0[i] - bias.fcn(ebinds[i],lambdas[i],lambda1s[i],lambda2s[i],alphas[i],u0s[i],w0s[i])
}

neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in 1:mlam) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=e0,ebind=ebinds,bet[te],lam[be],lambda1[be],lambda2[be],alpha[be],u0[be],w0coeff[be])
             sid <- sid + 1
    }
}




# note levels
label.tempt <- factor(tempts, levels=tempt, labels=1:mtempt)
label.lam <- factor(lambdas, levels=lam, labels=1:mlam)
label.cross <- (as.numeric(label.lam)-1)*mtempt + as.numeric(label.tempt)
out <- uwham.r(label=label.cross, logQ=neg.pot,ufactormax=1,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
ve <- matrix(out$ve, nrow=mtempt, ncol=mlam)
sqrt(ve)/bet


Vsite <- 4.*pi*( %f )^3/3.0
dgsite <- -log(Vsite/1668.0)/bet[]
dgbind <- dgsite + (-ze[,mlam]/bet[]) - (-ze[,1]/bet[])
ddgbind <- sqrt(ve[,mlam]+ve[,1])/bet

#average energy
usl1 <- ebinds[label.lam == mlam]
de <- mean(usl1)
sde <- sd(usl1)
dde <- sde/sqrt(length(usl1))
ub <- de + bet*sde*sde

#DGbind as a function of temperature
dgbind
sink("result.log")
cat("DGb = ", dgbind[1]," +- ",ddgbind[1]," DE = ", de, " +- ",dde,"\n")
cat("ub=", ub, "sd=", sde, "\n")
cat("dgsite=", dgsite)
sink()

dglambda <- cbind(lam,-ze[1,]/bet)
plot(lam, -ze[1,]/bet, type="l")
write(t(dglambda),file="dglambda.dat",ncol=2)

#get plain be histograms
umin <- min(ebinds)
umax <- max(ebinds)
u <- ebinds[label.lam == mlam];
hs <- hist(u,plot=FALSE);
pmax = 1.2*max(hs$density)
for ( i in 1:mlam ){
    u <- ebinds[label.lam == i];
    hs <- hist(u,plot=FALSE);
    if ( i == 1 ) {
        plot(hs$mids,hs$density,type="l",xlim=c(umin,umax),ylim=c(0,pmax));
    }else{
        lines(hs$mids,hs$density);
    }
    outp <- cbind(hs$mids,hs$density);
    #the percent sign is escaped below
    write(t(outp),file=sprintf("p-l%%e.dat",lam[i]),ncol=2)
}
