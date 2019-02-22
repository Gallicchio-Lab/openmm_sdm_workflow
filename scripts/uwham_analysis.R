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



data.t <- read.table("repl.cycle.totE.potE.temp.lambda.ebind.lambda1.lambda2.alpha.u0.w0.dat")

#quadbias
#gamma <- c(0.000, 0.033, 0.067, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.089, 0.067, 0.044, 0.022, 0.000)
#bcoeff <- c(0.000,-0.997,-1.996,-3.000,-2.413,-1.827,-1.240,-0.653,-0.067, 0.520, 1.107, 1.355, 1.266, 1.177, 1.089, 1.000)
#w0coeff <- c(0.000,16.956,33.945,51.000,35.556,22.772,12.648, 5.183, 0.378,-1.768,-1.254, 1.084, 0.521, 0.067,-0.175, 0.000)

#linear
#gamma <- c(0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000)
#bcoeff <- c(0.000, 0.067, 0.133, 0.200, 0.267, 0.333, 0.400, 0.467, 0.533, 0.600, 0.667, 0.733, 0.800, 0.867, 0.933, 1.000)
#w0coeff <- c(0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000)

#ilogistic
lam     <-c( 0.000, 0.057, 0.114, 0.171, 0.229, 0.286, 0.343, 0.400, 0.457, 0.514, 0.571, 0.629, 0.686, 0.743, 0.800, 1.000 )
lambda1 <-c( 0.000,-0.080,-0.160,-0.240,-0.320,-0.350,-0.350,-0.350,-0.350,-0.295,-0.076, 0.143, 0.362, 0.581, 0.800, 1.000 )
lambda2 <-c( 0.000, 0.183, 0.366, 0.549, 0.731, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 0.800, 1.000 )
alpha   <-c( 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.200 )
u0      <-c(12.000,12.000,12.000,12.000,12.000,10.714, 8.657, 6.600, 4.543, 3.000, 3.000, 3.000, 3.000, 3.000, 3.000, 3.000 )
w0coeff <-c( 0.000,-2.413,-4.827,-7.240,-9.653,-9.530,-7.884,-6.238,-4.593,-3.313,-3.130,-2.948,-2.765,-2.583,-2.400, 0.000 )


tempt <- c(300)
bet <- 1.0/(0.001986209*tempt)
mtempt <- length(bet)
mlam <- length(lam)
m <- mlam*mtempt
N <- length(data.t$V1)

#extract U0 values as U-bias
#this is relevant only if the states are at different temperatures
data.t$e0 <- data.t$V3
for (i in 1:N) {
    data.t$e0[i] <- data.t$e0[i] - bias.fcn(data.t$V7[i],data.t$V6[i],data.t$V8[i],data.t$V9[i],data.t$V10[i],data.t$V11[i],data.t$V12[i])
}

neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in 1:mlam) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=data.t$e0,ebind=data.t$V7,bet[te],lam[be],lambda1[be],lambda2[be],alpha[be],u0[be],w0coeff[be])
             sid <- sid + 1
    }
}




# note levels
label.tempt <- factor(data.t$V5, levels=tempt, labels=1:mtempt)
label.lam <- factor(data.t$V6, levels=lam, labels=1:mlam)
label.cross <- (as.numeric(label.lam)-1)*mtempt + as.numeric(label.tempt)
out <- uwham.r(label=label.cross, logQ=neg.pot,ufactormax=1,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
ve <- matrix(out$ve, nrow=mtempt, ncol=mlam)
sqrt(ve)/bet


Vsite <- 4.*pi*(2.5)^3/3.0
dgsite <- -log(Vsite/1668.0)/bet[]
dgbind <- dgsite + (-ze[,mlam]/bet[]) - (-ze[,1]/bet[])
ddgbind <- sqrt(ve[,mlam]+ve[,1])/bet

#average energy
de = mean(data.t$V7[data.t$V6 > 0.99])
sde <- sd(data.t$V7[data.t$V6 > 0.99])
dde = sde/sqrt(length(data.t$V7[data.t$V6 > 0.99]))
ub = de + bet*sde*sde

#DGbind as a function of temperature
dgbind
sink("result.log")
cat("DGb = ", dgbind[1]," +- ",ddgbind[1]," DE = ", de, " +- ",dde,"\n")
cat("ub=", ub, "sd=", sde, "\n")
sink()

dglambda <- cbind(lam,-ze[1,]/bet)
plot(lam, -ze[1,]/bet, type="l")
write(t(dglambda),file="dglambda.dat",ncol=2)

#get plain be histograms
umin <- min(data.t$V7)
umax <- max(data.t$V7)
u <- data.t$V7[label.lam == mlam];
hs <- hist(u,plot=FALSE);
pmax = 1.2*max(hs$density)
for ( i in 1:mlam ){
    u <- data.t$V7[label.lam == i];
    hs <- hist(u,plot=FALSE);
    if ( i == 1 ) {
        plot(hs$mids,hs$density,type="l",xlim=c(umin,umax),ylim=c(0,pmax));
    }else{
        lines(hs$mids,hs$density);
    }
    outp <- cbind(hs$mids,hs$density);
    write(t(outp),file=sprintf("p-l%e.dat",lam[i]),ncol=2)
}
