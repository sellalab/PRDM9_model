### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 6 ###
######################################
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")

sum.table = read.table("source_data/summary.table.txt",stringsAsFactors = FALSE)
pi.table = read.table("source_data/pi.table.txt", stringsAsFactors = FALSE)

#functions
Pf.finder = function(si, ki, pt, pf0 = 0, acc = -5:11) {
  pf = pf0
  for (i in 1:length(acc)) {
    accuracy = acc[i]
    while (pt - sum(si*4*pf/(pf+ki)) > pf) {
      pf = pf+1*10^(-accuracy)
    }
    pf = pf-1*10^(-accuracy)
  }
  dif = pf - (pt - sum(si*4*pf/(pf+ki)))
  pf = pf+1*10^(-acc[length(acc)])
  dif2 = pf - (pt - sum(si*4*pf/(pf+ki)))
  if (abs(dif) == min(abs(c(dif,dif2)))) {
    pf = pf-1*10^(-acc[length(acc)])
  }
  pf
}

W.hom = function(Si, Ki, Pt, D, r) {
  Pf = Pf.finder(Si,Ki,Pt)
  Hi = Pf/(Pf+Ki)
  c = D/(Pt-Pf)
  ai = 1-(1-c*(Hi^2)*(2-Hi)*(2-c*Hi))^2
  W = 1-prod((1-ai)^(Si*r))
  W
}
W.het.dup = function(Si, Ki, Pt, D, r) {
  Pf = Pf.finder(Si,Ki,Pt/2) # half PRDM9 expression
  Hi = Pf/(Pf+Ki)
  c = D/(Pt-(Pf*2)) # need to consider Pf from both alleles
  ai = 1-(1-c*(Hi^2)*(2-Hi)*(2-c*Hi))^2
  W = 1-prod((1-ai)^(Si*r))^2 # twice as many 'attempts' in terms of sites
  W
}
W.het.unique = function(SiA, SiB, Ki, Pt, D, r) {
  PfA = Pf.finder(SiA,Ki,Pt/2) # half PRDM9 expression
  PfB = Pf.finder(SiB,Ki,Pt/2) # half PRDM9 expression
  HiA = PfA/(PfA+Ki)
  HiB = PfB/(PfB+Ki)
  c = D/(Pt-PfA-PfB) # need to consider Pf from both alleles
  aiA = 1-(1-c*(HiA^2)*(2-HiA)*(2-c*HiA))^2
  aiB = 1-(1-c*(HiB^2)*(2-HiB)*(2-c*HiB))^2
  W = 1-prod((1-aiA)^(SiA*r))*prod((1-aiB)^(SiB*r)) # twice as many 'attempts' in terms of sites
  W
}

## calculating optimal number of sites in homozygotes and hets, and the number as which they can be invaded
Pt = 5000
D = 300       # Number of DSBs per meiosis
Pt = 5000     # Number of expressed PRDM9 molecules per meiosis
r = 1/40      # Proportion of genome on smallest chr

### Weak binding sites
S2 = 200000   # Number of weaker binding sites
PbE = 0.99    # Proportion of PRDM9 bound to weaker sites when there are no stronger ones

### Calculating k2 from the above
H2 = PbE*Pt/(4*S2)
Pf = Pt - PbE*Pt
k2 = Pf*(1-H2)/H2

opt.hom = NULL
opt.het = NULL
s2invasion = list()

for (i in 1:10) {
  k1 = i*5
  print(paste("Starting k1",k1))
  W.hom.temp = sapply(1:5000, FUN = function(S1) {
    W.hom(c(S1,S2),c(k1,k2),Pt,D,r)
  })
  opt.hom = c(opt.hom, which(W.hom.temp == max(W.hom.temp)))

  print(paste("..."))
  W.het.temp = sapply(1:5000, FUN = function(S1) {
    W.het.dup(c(S1,S2),c(k1,k2),Pt,D,r)
  })
  opt.het = c(opt.het, which(W.het.temp == max(W.het.temp)))

  print(paste("..."))
  W.het.invader = sapply(1:5000, FUN = function(S1) {
    W.het.unique(c(S1,S2),c(which(W.het.temp == max(W.het.temp)),S2),
                 c(k1,k2),Pt,D,r)
  })

  s2invasion[[i]] = which(W.het.invader > W.hom.temp)
}

par(mfrow=c(2,2))
par(mar = c(4,6,2,1))


#################################
### plot of S1i and S1e by k1 ###
#################################

temp.table = sum.table[which(sum.table$N == 1000),]

plot(x = temp.table$k1,
     y = temp.table$mean.S1i, t = "l", col="purple",
     xlab = expression(paste("Dissociation constant at hotspots (",italic("k" [1]),")")),
     ylab = expression(paste("Number of hotspots (",italic("n" [1]),")")),
     main = "",axes = "F",
     ylim = c(0,3000))

polygon(x = c((1:10)*5,(10:1)*5),
        y = c(opt.het,rev(opt.hom)), col = "azure", border = NA)

polygon(x = c((10:1)*5,(1:10)*5),
        y = c(rep(0,10),unlist(lapply(s2invasion, FUN = function(temp) {
          test = which(temp[2:length(temp)] - temp[1:(length(temp)-1)] != 1)
          if (length(test) > 0) {
            temp = test
          } else {
            temp = max(temp)
          }
          temp
        }))), col = "mistyrose", border = NA)

red.y.coords = unlist(lapply(s2invasion, FUN = function(temp) {
  test = which(temp[2:length(temp)] - temp[1:(length(temp)-1)] != 1)
  if (length(test) > 0) {
    temp = test
  } else {
    temp = max(temp)
  }
  temp
}))[1:2]
blue.y.coords = opt.het[1:2]
new.x = 318/46.4
new.y = -20.2*new.x + 1075

polygon(x = c(5, new.x, 5),
        y = c(opt.het[1], new.y, red.y.coords[1]), col="plum2", border=NA)

points(x = (1:10)*5,
       y = unlist(lapply(s2invasion, FUN = function(temp) {
         test = which(temp[2:length(temp)] - temp[1:(length(temp)-1)] != 1)
         if (length(test) > 0) {
           temp = test
         } else {
           temp = max(temp)
         }
         temp
       })), t = "l", lty=2, col="mistyrose3")
points(x = (1:10)*5, y = opt.het, col = "azure3",lty=2,t="l")
points(x = (1:10)*5, y = opt.hom, col = "azure3",lty=2,t="l")

points(x = temp.table$k1, y = temp.table$mean.S1i, t = "l", col = "blue")

temp.table = sum.table[which(sum.table$N == round(10^4)),]
points(x = temp.table$k1, y = temp.table$mean.S1i, t = "l", col = "green")

temp.table = sum.table[which(sum.table$N == round(10^5)),]
points(x = temp.table$k1, y = temp.table$mean.S1i, t = "l", col = "orange")

temp.table = sum.table[which(sum.table$N == round(10^6)),]
points(x = temp.table$k1, y = temp.table$mean.S1i, t = "l", col = "red")

temp.table = sum.table[which(sum.table$N == round(10^3)),]
points(x = temp.table$k1, y = temp.table$mean.S1e, t = "l", lty = 2, col = "blue")

temp.table = sum.table[which(sum.table$N == round(10^4)),]
points(x = temp.table$k1, y = temp.table$mean.S1e, t = "l", lty = 2, col = "green")

temp.table = sum.table[which(sum.table$N == round(10^5)),]
points(x = temp.table$k1, y = temp.table$mean.S1e, t = "l", lty = 2, col = "orange")

temp.table = sum.table[which(sum.table$N == round(10^6)),]
points(x = temp.table$k1, y = temp.table$mean.S1e, t = "l", lty = 2, col = "red")

legend(5,3000,c(expression(10 ^6), expression(10 ^5), expression(10 ^4), expression(10 ^3)),
       title = expression(paste("Pop size (",italic("N"),")")),cex = 1,bty="n",
       lty = rep(1,7),col=c("red","orange","green","blue"))
legend(20,3000,c("invading alleles","exiting alleles"),lty=c(1,2),title=expression(paste("Mean ",italic("n" [1])," among...")),bty="n",cex=1)

legend(30,1500,c("Susceptible \n to invasion","Optimal hotspots","Both"),bty="n",fill=c("mistyrose","azure","plum2"),border=NA, cex = 1)

box()
axis(side=2, at=c(0,3000),las=1)
axis(side=1, at=c(5,20,35,50),las=1)

#################################
### plot turnover by log10(N) ###
#################################

k1 = 5

temp.table = sum.table[which(sum.table$k1 == k1),]
N.options = c(3,3.5,4,4.5,5,5.5,6)

plot(x = N.options,
     y = log10(temp.table$mean.tau),
     t = "l",col="red",ylim=c(1.5,5.3),
     main = "",axes = "F",
     xlab = expression(paste("Population Size (",italic("N"),")")),
     ylab = expression(paste("Turnover time (log" [10]," scale)",sep="")))

k1 = 20
temp.table = sum.table[which(sum.table$k1 == k1),]
points(x = N.options, y = log10(temp.table$mean.tau),
       t = "l", col="orange")


k1 = 35
temp.table = sum.table[which(sum.table$k1 == k1),]
points(x = N.options, y = log10(temp.table$mean.tau),
       t = "l", col="green")

k1 = 50
temp.table = sum.table[which(sum.table$k1 == k1),]
points(x = N.options, y = log10(temp.table$mean.tau),
       t = "l", col="blue")

abline(a = 8, b=-1, lty=2)

legend(5.3,5.3,c(5,20,35,50), title = expression(italic("k" [1])),
       lty = rep(1,10), col=c("red","orange",
                              "green","blue"), bty = "n")

box()
axis(side=2, at=c(2,3,4,5),las=1, labels = c(expression(10 ^2),
                                             expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5)))
axis(side=1, at=c(3,4,5,6),las=1, labels = c(expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5),
                                             expression(10 ^6)))

##############################
### plot of pi by log10(N) ###
##############################

k1 = 5
temp = pi.table[which(pi.table$k1 == k1),]
plot(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$pi, t= "l", col = "red", ylim = c(0,1),lty=2,
     ylab = expression(paste("PRDM9 Diversity (",pi,")")), xlab = expression(paste("Population size (",italic("N"),")",sep="")),
     main = "", axes = "F")
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$pi.25,rev(temp$pi.75)), col = rgb(1,0,0,0.1), border = NA)

k1 =20
temp = pi.table[which(pi.table$k1 == k1),]
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$pi.25,rev(temp$pi.75)), col = rgb(1,0.5,0,0.1), border = NA)

k1 =35
temp = pi.table[which(pi.table$k1 == k1),]
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$pi.25,rev(temp$pi.75)), col = rgb(0,1,0,0.1), border = NA)

k1 =50
temp = pi.table[which(pi.table$k1 == k1),]
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$pi.25,rev(temp$pi.75)), col = rgb(0,0,1,0.1), border = NA)

k1 = 5
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$pi, t= "l", col = "red",lwd=1,lty=1)
k1 = 20
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$pi, t= "l", col = "orange",lwd=1,lty=1)
k1 = 35
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$pi, t= "l", col = "green",lwd=1,lty=1)
k1 = 50
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$pi, t= "l", col = "blue",lwd=1,lty=1)



legend(3,0.8,c(5,20,35,50),lty=rep(1,4),lwd=rep(1,4),
       col=c("red","orange","green","blue"),
       title=expression(italic("k" [1])),bty="n")

legend(3,1, c("mean diversity"),
       lty=1,lwd=1,bty="n")
legend(3,0.93, c("Interquantile range"),
       fill = rgb(1,0,0,0.1), bty="n", border=NA)

box()
axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(3,4,5,6),las=1, labels = c(expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5),
                                             expression(10 ^6)))

##############################
### plot of W by log10(N) ###
##############################

k1 = 5
temp = pi.table[which(pi.table$k1 == k1),]

plot(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "red", ylim = c(.65,1),lty=2,
     ylab = "Fitness (W)", xlab = expression(paste("Population size (",italic("N"),")",sep="")),
     main = "", axes = "F")
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$W.25,rev(temp$W.75)), col = rgb(1,0,0,0.1), border = NA)

k1 =20
temp = pi.table[which(pi.table$k1 == k1),]
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$W.25,rev(temp$W.75)), col = rgb(1,0.5,0,0.1), border = NA)
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "orange",lwd=1,lty=1)

k1 =35
temp = pi.table[which(pi.table$k1 == k1),]
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$W.25,rev(temp$W.75)), col = rgb(0,1,0,0.1), border = NA)
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "green",lwd=1,lty=1)

k1 =50
temp = pi.table[which(pi.table$k1 == k1),]
polygon(x = c(c(3,3.5,4,4.5,5,5.5,6), rev(c(3,3.5,4,4.5,5,5.5,6))),
        y = c(temp$W.25,rev(temp$W.75)), col = rgb(0,0,1,0.1), border = NA)
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "blue",lwd=1,lty=1)

k1 = 5
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "red",lwd=1,lty=1)
k1 = 20
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "orange",lwd=1,lty=1)
k1 = 35
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "green",lwd=1,lty=1)
k1 = 50
temp = pi.table[which(pi.table$k1 == k1),]
points(x = c(3,3.5,4,4.5,5,5.5,6), y = temp$W, t= "l", col = "blue",lwd=1,lty=1)


legend(3,0.75,c(5,20,35,50),lty=rep(1,4),lwd=rep(1,4),
       col=c("red","orange","green","blue"),
       title=expression(italic("k" [1])),bty="n")

legend(3.75,.7, c("mean fitness"),
       lty=1,lwd=1,bty="n")
legend(3.75,0.68, c("Interquantile range"),
       fill = rgb(1,0,0,0.1), bty="n", border=NA)

box()
axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(3,4,5,6),las=1, labels = c(expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5),
                                             expression(10 ^6)))



