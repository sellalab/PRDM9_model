### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 6-2 ###
######################################
# Fitness in an individual homozygous for a PRDM9 allele, or heterozygous for two PRDM9 alleles with the same binding distribution
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")
par(mfrow=c(5,2))

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

per.gen.3.50 = read.table("source_data/per.gen.results.1000_50.txt",
                          nrows = 2e6, skip = 2e6)
colnames(per.gen.3.50) = c("S1","W","pi")

per.gen.3.5 = read.table("source_data/per.gen.results.1000_5.txt",
                         nrows = 5e5, skip = 2e6)
colnames(per.gen.3.5) = c("S1","W","pi")

par(mfrow=c(6,2))

par(mar = c(0,10,2,0))
plot(per.gen.3.50$pi, ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
mtext(expression(paste(italic("k" [1])," = 50")), side=3, line=-0.1, cex=1.2)
axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(1,2e6), labels = c("1","2000000"))

par(mar = c(0,0,2,10))
plot(per.gen.3.5$pi, ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
mtext(expression(paste(italic("k" [1])," = 5")), side=3, line=-0.1, cex=1.2)
axis(side=1, at=c(1,5e5),labels = c("1","500000"))


### mutation at prdm9 (v) only
if (TRUE) {
  per.gen.3.50 = read.table("source_data/per.gen.results.control_v_50.txt")
  colnames(per.gen.3.50) = c("S1","W","pi")

  per.gen.3.5 = read.table("source_data/per.gen.results.control_v_5.txt",
                           skip = 5e5)
  colnames(per.gen.3.5) = c("S1","W","pi")

  par(mar = c(0,10,2,0))
  plot(per.gen.3.50$pi[(1e5 + 1):(2.1e6)], ylim = c(0,1), t = "l",
       ylab = "", axes = "F", xlab = "")
  axis(side=2, at=c(0,1),las=1)
  axis(side=1, at=c(1,2e6), labels = c("1","2000000"))

  par(mar = c(0,0,2,10))
  plot(per.gen.3.5$pi[(1):(5e5)], ylim = c(0,1), t = "l",
       ylab = "", axes = "F", xlab = "")
  axis(side=1, at=c(1,5e5),labels = c("1","500000"))



}

### mutation at prdm9 binding sites (u) only

per.gen.3.50 = read.table("source_data/per.gen.results.control_u_50.txt",
                          nrows = 2e6, skip = 2e6)
colnames(per.gen.3.50) = c("S1","W","pi")

per.gen.3.5 = read.table("source_data/per.gen.results.control_u_5.txt",
                         nrows = 5e5, skip = 2e6)
colnames(per.gen.3.5) = c("S1","W","pi")

par(mar = c(0,10,2,0))
plot(per.gen.3.50$pi[1:4000], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(1,4e3))
mtext(expression(paste("Diversity (",italic(pi),")",sep="")), side=2, line=1.5, cex=1.1, las=1)

par(mar = c(0,0,2,10))
plot(per.gen.3.5$pi[1:2500], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=1, at=c(1,2500))


# mutation at both
per.gen.3.50 = read.table("source_data/per.gen.results.both_5.txt")
colnames(per.gen.3.50) = c("S1","W","pi")

per.gen.3.5 = read.table("source_data/per.gen.results.both_5.txt")
colnames(per.gen.3.5) = c("S1","W","pi")

par(mar = c(0,10,2,0))
plot(per.gen.3.50$pi[50001:54000], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(1,4000))


par(mar = c(0,0,2,10))
plot(per.gen.3.5$pi[50001:52500], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=1, at=c(1,2500))


# mutation at both
per.gen.6.50 = read.table("source_data/per.gen.results.1e+06_50.txt")
colnames(per.gen.6.50) = c("S1","W","pi")

per.gen.6.5 = read.table("source_data/per.gen.results.1e+06_5.txt")
colnames(per.gen.6.5) = c("S1","W","pi")

par(mar = c(0,10,2,0))
plot(per.gen.6.50$pi[50001:54000], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(1,4000))
mtext("Time (generations)", side=1, line=0.8, cex=1.0)

par(mar = c(0,0,2,10))
plot(per.gen.6.5$pi[50001:52500], ylim = c(0,1), t = "l",
     ylab = "", axes = "F", xlab = "")
axis(side=1, at=c(1,2500))
mtext("Time (generations)", side=1, line=0.8, cex=1.0)


