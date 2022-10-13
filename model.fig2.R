### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 2 ###
######################################
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")
par(mfrow=c(2,1))
par(mar = c(4.5,4.5,0.5,1))

###############
### Panel A ###
###############
# Cartoon of fitness as a function of time under different assumptions.

plot(x = 1:100, y=rep(0,100),
     t= "l", col = "blue", lwd=3,
     xlim = c(1,100), ylim = c(-75,100),
     main = "", ylab = "", xlab = "",
     xaxt="n", yaxt="n")

points(x = c(1,100), y = c(0,40),
       t="l", col="red", lwd = 3)
points(x = c(1,100), y = c(0,-40),
       t="l", col="green", lwd = 3)

legend(1,100,c("non-competitive binding, variable DSBs",
               "non-competitive binding, constant DSBs",
               "competitive binding, constant DSBs"),
       lty = c(1,1,1), lwd = c(3,3,3), col=c("green","blue","red"))

title(ylab = substitute(paste("Fitness (",italic("W"),")",sep="")), cex.lab = 1.2, line=1)
title(xlab = "Time", cex.lab = 1.2, line=1)

###############
### Panel B ###
###############
# The probability of an individual binding site being bound (thicker lines) or
# of a given locus being bound symmetrically (thinner line) if all sites are
# very weak (blue) or very strong (red), i.e., with very large or small
# dissociation constants respectively. For sake of comparison, in both cases, we
# set the number of expressed PRDM9 molecules such that 5,000 sites would be
# bound in the presence of 20,000 binding sites, roughly consistent with
# observations in mice (Parvanov et al. 2017; Grey et al. 2017). When binding
# sites are very weak, most PRDM9 molecules are not bound and therefore our
# model behaves as if there were no competition between binding sites.

H = 5000/(20000*4)

### Modeling competitive binding

# To calculate the probability of being bound, or of being bound symmetrically,
# we first have to calculate what total concentration of PRDM9 is necessary to
# result in 5000 bound sites when there are 20000 possible binding sites.
k = 1
Pf = H*k/(1-H)
Pt = 5000+Pf

# Given the total concentration, we can calculate the free concentration of PRDM9
# as a function of the number of sites remaining.
pf.competitive = sapply(1:20000, FUN = function(x) {
  calculate.pf(x, k, Pt)
})

# From this, we can calculate the probability of a remaining site being bound
# as a function of the number of sites remaining
Hi.competitive = pf.competitive/(k+pf.competitive)

### Modeling non-competitive binding
# this follows the same steps taken about for k=1, except now k=1e6

k = 1e6
Pf = H*k/(1-H)
Pt = 5000+Pf
pf.noncompetitive = sapply(1:20000, FUN = function(x) {
  pf = calculate.pf(x, k, Pt)
})
Hi.noncompetitive = pf.noncompetitive/(k+pf.noncompetitive)


### Calculating probabilities of being symmetrically bound
sym.competitive = (1-(1-Hi.competitive)^2)^2
sym.noncompetitive = (1-(1-Hi.noncompetitive^2))^2

### Plotting the above
plot(Hi.competitive, t="l", col="red", lwd=3,
     xlab = "", ylab = "", ylim=c(0,1))
points(sym.competitive, t="l", col="red")
points(Hi.noncompetitive, t="l", col="blue", lwd=3)
points(sym.noncompetitive, t="l", col="blue")

legend(5000,1,c("of being bound","of being bound symmetrically",
                expression(paste("non-competitive (",italic("k"),"=10"^"6",")")),
                substitute(paste("competitive (",italic("k"),"=1)",sep=""))),
       lty = c(1,1,1,1), lwd = c(3,1,3,3), col=c("black","black","blue","red"))

title(ylab = "Probability", cex.lab = 1.2, line=3)
title(xlab = substitute(paste("Number of sites (",italic("n"),")",sep="")), cex.lab = 1.2, line=3)



