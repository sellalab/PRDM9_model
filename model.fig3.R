### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 3 ###
######################################
# Fitness in an individual homozygous for a PRDM9 allele, or heterozygous for two PRDM9 alleles with the same binding distribution
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")
par(mfrow=c(2,1))
par(mar = c(4.5,4.5,0.5,1))

## calculating optimal number of sites in homozygotes and hets, and the number as which they can be invaded
D = 300       # Number of DSBs per meiosis
Pt = 5000     # Number of expressed PRDM9 molecules per meiosis
r = 1/40      # Proportion of genome on smallest chr

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

hom.W = sapply(1:5000, FUN = function(s1) {
  calculate.W.hom(c(s1,S2), c(k1,k2), Pt, D, r)
})

het.W = sapply(1:5000, FUN = function(s1) {
  calculate.W.het(c(s1,S2), c(s1,S2), c(k1,k2), Pt, D, r)
})

invade.W = sapply(1:5000, FUN = function(s1) {
  calculate.W.het(c(s1,S2), c(which(het.W == max(het.W)),S2), c(k1,k2), Pt, D, r)
})

plot(hom.W, t="l", col="blue",lwd=3, ylim=c(0,1),
     ylab = "",
     xlab = "")
title(ylab = substitute(paste("Fitness (",italic("W"),")")), cex.lab = 1.2, line=3)
title(xlab = expression(paste("Number of hotspots (",italic('n' [1]),")",sep="")), cex.lab = 1.2, line=3)

polygon(x = c(-1000,-1000,which(invade.W > hom.W)[length(which(invade.W > hom.W))],which(invade.W > hom.W)[length(which(invade.W > hom.W))]),
        y = c(2,-1,-1,2), col = "mistyrose", border = NA)
box()

points(hom.W, t="l", col="blue",lwd=3)
points(het.W, t="l", col="red",lwd=3)

abline(v = which(het.W == max(het.W)),
       col="red",lty =2,lwd=2)
abline(v = which(hom.W == max(hom.W)),
       col="blue",lty =2,lwd=2)
abline(v = which(invade.W > hom.W)[length(which(invade.W > hom.W))],
       lty = 2, col="mistyrose3", lwd = 2)

###############
### Panel B ###
###############
# when considering strong hotspots
k1 = 5

hom.W = sapply(1:5000, FUN = function(s1) {
  calculate.W.hom(c(s1,S2), c(k1,k2), Pt, D, r)
})

het.W = sapply(1:5000, FUN = function(s1) {
  calculate.W.het(c(s1,S2), c(s1,S2), c(k1,k2), Pt, D, r)
})

invade.W = sapply(1:5000, FUN = function(s1) {
  calculate.W.het(c(s1,S2), c(which(het.W == max(het.W)),S2), c(k1,k2), Pt, D, r)
})

plot(hom.W, t="l", col="blue",lwd=3, ylim=c(0,1),
     ylab = "",
     xlab = "")
title(ylab = substitute(paste("Fitness (",italic("W"),")")), cex.lab = 1.2, line=3)
title(xlab = expression(paste("Number of hotspots (",italic('n' [1]),")",sep="")), cex.lab = 1.2, line=3)

polygon(x = c(-1000,-1000,974,974),
        y = c(2,-1,-1,2), col = "mistyrose", border = NA)
polygon(x = c(4913,4913,6000,6000),
        y = c(2,-1,-1,2), col = "mistyrose", border = NA)
box()

points(hom.W, t="l", col="blue",lwd=3)
points(het.W, t="l", col="red",lwd=3)

abline(v = which(het.W == max(het.W)),
       col="red",lty =2,lwd=2)
abline(v = which(hom.W == max(hom.W)),
       col="blue",lty =2,lwd=2)
abline(v = 974,
       lty = 2, col="mistyrose3", lwd = 2)
abline(v = 4913,
       lty = 2, col="mistyrose3", lwd = 2)

legend(2500,.75,c("Homozygotes","Heterozygotes"),lt=c(1,1),lwd=c(3,3),col=c("blue","red"))



