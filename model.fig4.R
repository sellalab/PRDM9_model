### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 4 ###
######################################
# Fitness in an individual homozygous for a PRDM9 allele, or heterozygous for two PRDM9 alleles with the same binding distribution
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")

par(mfrow=c(6,2))

D = 300       # Number of DSBs per meiosis
Pt = 5000     # Number of expressed PRDM9 molecules per meiosis
r = 1/40      # Proportion of genome on smallest chr
S2 = 200000   # Number of weaker binding sites
PbE = 0.99    # Proportion of PRDM9 bound to weaker sites when there are no stronger ones

### Calculating k2 from the above
H2 = PbE*Pt/(4*S2)
Pf = Pt - PbE*Pt
k2 = Pf*(1-H2)/H2

###############
### k1 = 50 ###
###############

### loading per generation results
per.gen.6.50 = read.table("source_data/per.gen.results.1e+06_50.txt")
colnames(per.gen.6.50) = c("S1","W","pi")

per.gen.3.50 = read.table("source_data/per.gen.results.1000_50.txt",
                          nrows = 2e6, skip = 2e6)
colnames(per.gen.3.50) = c("S1","W","pi")

temp.6.50 = (length(per.gen.6.50$S1)-1999):length(per.gen.6.50$S1)
temp.3.50 = 1:2e6

### calculating optimal numbers of hotspots...
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

opt.hom = which(hom.W == max(hom.W))
opt.het = which(het.W == max(het.W))
t = which(invade.W > hom.W)
invade.S = t[length(t)]

### mean number of sites for N = 10^6
par(mar = c(0,10,2,0))

plot(per.gen.6.50$S1[temp.6.50], ylim = c(0,3000), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,3000),las=1)
mtext(expression(paste("Hotspots (",italic("n" [1]),")",sep="")), side=2, line=1.5, cex=1.1, las=1)
mtext(expression(paste(italic("N")," = 10" ^6)), side=3, line=-0.1, cex=1.2)

polygon(y = c(opt.hom,opt.hom,opt.het,opt.het),
        x = c(-1000,4000,4000,-1000), col="azure",border = NA)
abline(h = opt.hom, lty = 2, col = "azure3")
abline(h = opt.het, lty = 2, col = "azure3")

polygon(y = c(invade.S,invade.S,0,0),
        x = c(-1000,4000,4000,-1000), col="mistyrose",border = NA)
abline(h = invade.S, lty = 2, col = "mistyrose3")

points(per.gen.6.50$S1[temp.6.50], t = "l")

### mean number of sites for N = 10^3
par(mar = c(0,0,2,10))

plot(per.gen.3.50$S1[temp.3.50], ylim = c(0,3000), t = "l",
     ylab = "", axes = "F", xlab = "")
mtext(expression(paste(italic("N")," = 10" ^3)), side=3, line=-0.1, cex=1.2)

polygon(y = c(opt.hom,opt.hom,opt.het,opt.het),
        x = c(-0,3e6,3e6,-0), col="azure",border = NA)
abline(h = opt.hom, lty = 2, col = "azure3")
abline(h = opt.het, lty = 2, col = "azure3")

polygon(y = c(invade.S,invade.S,0,0),
        x = c(-0,3e6,3e6,-0), col="mistyrose",border = NA)
abline(h = invade.S, lty = 2, col = "mistyrose3")

points(per.gen.3.50$S1[temp.3.50], t = "l")

## diversity for N = 10^6
par(mar = c(1,10,1,0))
plot(per.gen.6.50$pi[temp.6.50], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,1),las=1)
mtext(expression(paste("Diversity (",italic(pi),")",sep="")), side=2, line=1.5, cex=1.1, las=1)

abline(h = mean(per.gen.6.50$pi[temp.6.50]), col = "red", lty = 2)

## diversity for N = 10^3
par(mar = c(1,0,1,10))
plot(per.gen.3.50$pi[temp.3.50], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
abline(h = mean(per.gen.3.50$pi[temp.3.50]), col = "red", lty = 2)

## fitness for N = 10^6
par(mar = c(2,10,0,0))
plot(per.gen.6.50$W[temp.6.50], ylim = c(0.5,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0.5,1),las=1)
mtext(expression(paste("Fitness (",italic("W"),")",sep="")), side=2, line=1.5, cex=1.1, las=1)
abline(h = max(hom.W), col = "blue")
abline(h = max(het.W), col = "red")

axis(side=1, at=c(1,2000))
mtext("Time (generations)", side=1, line=0.8, cex=1.0)

## fitness for N = 10^3
par(mar = c(2,0,0,10))
plot(per.gen.3.50$W[temp.3.50], ylim = c(0.5,1), t = "l",
     ylab = "", axes = "F", xlab = "")
abline(h = max(hom.W), col = "blue")
abline(h = max(het.W), col = "red")

axis(side=1, at=c(1,2e6), labels = c("1","2,000,000"))
mtext("Time (generations)", side=1, line=0.8, cex=1.0)



###############
### k1 = 5 ###
###############

### loading per generation results
per.gen.6.5 = read.table("source_data/per.gen.results.1e+06_5.txt")
colnames(per.gen.6.5) = c("S1","W","pi")

per.gen.3.5 = read.table("source_data/per.gen.results.1000_5.txt",
                          nrows = 5e5, skip = 2e6)
colnames(per.gen.3.5) = c("S1","W","pi")

temp.6.5 = (length(per.gen.6.5$S1)-499):length(per.gen.6.5$S1)
temp.3.5 = 1:5e5

### calculating optimal numbers of hotspots...
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

opt.hom = which(hom.W == max(hom.W))
opt.het = which(het.W == max(het.W))
t = which(invade.W > hom.W)
invade.S = 974

### mean number of sites for N = 10^6
par(mar = c(0,10,2,0))

plot(per.gen.6.5$S1[temp.6.5], ylim = c(0,3000), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,3000),las=1)
mtext(expression(paste("Hotspots (",italic("n" [1]),")",sep="")), side=2, line=1.5, cex=1.1, las=1)

polygon(y = c(opt.hom,opt.hom,opt.het,opt.het),
        x = c(-1000,4000,4000,-1000), col="azure",border = NA)
polygon(y = c(invade.S,invade.S,0,0),
        x = c(-1000,4000,4000,-1000), col="mistyrose",border = NA)
polygon(y = c(invade.S,invade.S,opt.het,opt.het),
        x = c(-1000,4000,4000,-1000), col="plum2",border = NA)

abline(h = opt.hom, lty = 2, col = "azure3")
abline(h = opt.het, lty = 2, col = "azure3")
abline(h = invade.S, lty = 2, col = "mistyrose3")

points(per.gen.6.5$S1[temp.6.5], t = "l")

### mean number of sites for N = 10^3
par(mar = c(0,0,2,10))

plot(per.gen.3.5$S1[temp.3.5], ylim = c(0,3000), t = "l",
     ylab = "", axes = "F", xlab = "")

polygon(y = c(opt.hom,opt.hom,opt.het,opt.het),
        x = c(-1000,1e6,1e6,-1000), col="azure",border = NA)
polygon(y = c(invade.S,invade.S,0,0),
        x = c(-1000,1e6,1e6,-1000), col="mistyrose",border = NA)
polygon(y = c(invade.S,invade.S,opt.het,opt.het),
        x = c(-1000,1e6,1e6,-1000), col="plum2",border = NA)

abline(h = opt.hom, lty = 2, col = "azure3")
abline(h = opt.het, lty = 2, col = "azure3")
abline(h = invade.S, lty = 2, col = "mistyrose3")

points(per.gen.3.5$S1[temp.3.5], t = "l")

## diversity for N = 10^6
par(mar = c(1,10,1,0))
plot(per.gen.6.5$pi[temp.6.5], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,1),las=1)
mtext(expression(paste("Diversity (",italic(pi),")",sep="")), side=2, line=1.5, cex=1.1, las=1)

abline(h = mean(per.gen.6.5$pi[temp.6.5]), col = "red", lty = 2)

## diversity for N = 10^3
par(mar = c(1,0,1,10))
plot(per.gen.3.5$pi[temp.3.5], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
abline(h = mean(per.gen.3.5$pi[temp.3.5]), col = "red", lty = 2)

## fitness for N = 10^6
par(mar = c(2,10,0,0))
plot(per.gen.6.5$W[temp.6.5], ylim = c(0.95,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0.95,1),las=1)
mtext(expression(paste("Fitness (",italic("W"),")",sep="")), side=2, line=1.5, cex=1.1, las=1)
abline(h = max(hom.W), col = "blue")
abline(h = max(het.W), col = "red")

axis(side=1, at=c(1,500))
mtext("Time (generations)", side=1, line=0.8, cex=1.0)

## fitness for N = 10^3
par(mar = c(2,0,0,10))
plot(per.gen.3.5$W[temp.3.5], ylim = c(0.95,1), t = "l",
     ylab = "", axes = "F", xlab = "")
abline(h = max(hom.W), col = "blue")
abline(h = max(het.W), col = "red")

axis(side=1, at=c(1,5e5), labels = c("1","500,000"))
mtext("Time (generations)", side=1, line=0.8, cex=1.0)









