### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
###############################################
### Script for generating Figure 3 - Supp 1 ###
###############################################
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")
par(mfrow=c(2,1))
par(mar = c(4.5,4.5,0.5,1))

## calculating optimal number of sites in homozygotes and hets, and the number as which they can be invaded
D = 300       # Number of DSBs per meiosis
Pt = 5000     # Number of expressed PRDM9 molecules per meiosis

### Weak binding sites
S2 = 200000   # Number of weaker binding sites
PbE = 0.99    # Proportion of PRDM9 bound to weaker sites when there are no stronger ones
H2 = PbE*Pt/(4*S2)
Pf = Pt - PbE*Pt
k2 = Pf*(1-H2)/H2

###############
### Panel A ###
###############
# when considering weak hotspots
k1 = 50

### homozygotes
pf = sapply(1:5000, FUN = function(s1) {
  calculate.pf(c(s1,S2),c(k1,k2),Pt=Pt)
})
H1 = pf/(pf+k1)
H2 = pf/(pf+k2)

Pb1 = H1*4*(1:5000)
Pb2 = H2*4*S2

p.DSB = Pb1/(Pb1+Pb2)
p.sym = (1-(1-H1)^2)^2

plot(p.DSB, t="l", col="blue",lwd=3, ylim=c(0,1),
     ylab = "",
     xlab = "")
title(ylab = "Probability", cex.lab = 1.2, line=3)
title(xlab = expression(paste("Number of hotspots (",italic('n' [1]),")",sep="")), cex.lab = 1.2, line=3)

points(p.sym, t="l", col="red", lwd=3)

### heterozygotes
pf = sapply(1:5000, FUN = function(s1) {
  calculate.pf(c(s1,S2),c(k1,k2),Pt=Pt/2)
})
H1 = pf/(pf+k1)
H2 = pf/(pf+k2)

Pb1 = H1*4*(1:5000)
Pb2 = H2*4*S2

p.DSB = Pb1/(Pb1+Pb2)
p.sym = (1-(1-H1)^2)^2

points(p.DSB, t="l", col="blue", lwd=3, lty=3)
points(p.sym, t="l", col="red", lwd=3, lty=3)

###############
### Panel B ###
###############
# when considering weak hotspots
k1 = 5

### homozygotes
pf = sapply(1:5000, FUN = function(s1) {
  calculate.pf(c(s1,S2),c(k1,k2),Pt=Pt)
})
H1 = pf/(pf+k1)
H2 = pf/(pf+k2)

Pb1 = H1*4*(1:5000)
Pb2 = H2*4*S2

p.DSB = Pb1/(Pb1+Pb2)
p.sym = (1-(1-H1)^2)^2

plot(p.DSB, t="l", col="blue",lwd=3, ylim=c(0,1),
     ylab = "",
     xlab = "")
title(ylab = "Probability", cex.lab = 1.2, line=3)
title(xlab = expression(paste("Number of hotspots (",italic('n' [1]),")",sep="")), cex.lab = 1.2, line=3)

points(p.sym, t="l", col="red", lwd=3)

### heterozygotes
pf = sapply(1:5000, FUN = function(s1) {
  calculate.pf(c(s1,S2),c(k1,k2),Pt=Pt/2)
})
H1 = pf/(pf+k1)
H2 = pf/(pf+k2)

Pb1 = H1*4*(1:5000)
Pb2 = H2*4*S2

p.DSB = Pb1/(Pb1+Pb2)
p.sym = (1-(1-H1)^2)^2

points(p.DSB, t="l", col="blue", lwd=3, lty=3)
points(p.sym, t="l", col="red", lwd=3, lty=3)




