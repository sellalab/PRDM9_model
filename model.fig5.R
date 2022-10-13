### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 5 ###
######################################
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")

## calculating optimal number of sites in homozygotes and hets, and the number as which they can be invaded
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

for (k1 in c(5,20,35,50)) {
  #print(paste("Starting k1",k1))
  W.hom.temp = sapply(1:5000, FUN = function(S1) {
    calculate.W.hom(c(S1,S2),c(k1,k2),Pt,D,r)
  })
  opt.hom = c(opt.hom, which(W.hom.temp == max(W.hom.temp)))

  W.het.temp = sapply(1:5000, FUN = function(S1) {
    calculate.W.het(c(S1,S2),c(S1,S2),c(k1,k2),Pt,D,r)
  })
  opt.het = c(opt.het, which(W.het.temp == max(W.het.temp)))

  W.het.invader = sapply(1:5000, FUN = function(S1) {
    calculate.W.het(c(S1,S2),c(which(W.het.temp == max(W.het.temp)),S2),
                 c(k1,k2),Pt,D,r)
  })
  s2invasion[[length(s2invasion)+1]] = which(W.het.invader > W.hom.temp)
}

par(mfrow=c(4,1))

S1i.table = read.table("source_data/S1i.table.txt")
S1e.table = read.table("source_data/S1e.table.txt")

#######################################################################
### plot ###
#######################################################################
par(mfrow=c(4,1))

##### K1 = 5
par(mar = c(0,5,0,1))

for (i in 1:4) {
  k1 = c(5,20,35,50)[i]
  t.invade = S1i.table[which(S1i.table[,2] == k1),]
  t.exit = S1e.table[which(S1e.table[,2] == k1),]

  plot(y = t.invade[1,3:(dim(t.invade)[2])],
       x = (1:50)*100 -50,
       xlim=c(1,4000),
       ylim=c(0,0.3),t="l",xaxt="n",cex.lab=1.3,bty="n",yaxt="n",
       ylab = "",xlab = "",col="blue",srt=90)

  if (i == 1) {
    text(x = -250,y = 0.15,
         labels = expression(paste(italic('k' [1]),"=",5,sep="")),
         xpd = NA,srt = 0,cex = 1.6)
  }
  if (i == 2) {
    text(x = -250,y = 0.15,
         labels = expression(paste(italic('k' [1]),"=",20,sep="")),
         xpd = NA,srt = 0,cex = 1.6)
  }
  if (i == 3) {
    text(x = -250,y = 0.15,
         labels = expression(paste(italic('k' [1]),"=",35,sep="")),
         xpd = NA,srt = 0,cex = 1.6)
  }
  if (i == 4) {
    text(x = -250,y = 0.15,
         labels = expression(paste(italic('k' [1]),"=",50,sep="")),
         xpd = NA,srt = 0,cex = 1.6)
  }

  polygon(x = c(opt.hom[i],opt.hom[i],opt.het[i],opt.het[i]),
          y = c(-1,1,1,-1), col="azure",border = NA)
  abline(v = opt.hom[i], lty = 2, col = "azure3")
  abline(v = opt.het[i], lty = 2, col = "azure3")

  temp = s2invasion[[i]]
  test = which(temp[2:length(temp)] - temp[1:(length(temp)-1)] != 1)
  if (length(test) > 0) {
    temp = test
  } else {
    temp = max(temp)
  }
  polygon(x = c(0,0,temp,temp),
          y = c(-1,1,1,-1), col="mistyrose",border = NA)
  abline(v = temp, lty = 2, col = "mistyrose3")

  if (temp > opt.het[i]) {
    polygon(x = c(opt.het[i],opt.het[i],temp,temp),
            y = c(-1,1,1,-1), col="plum2",border = NA)
    abline(v = temp, lty = 2, col = "mistyrose3")
    abline(v = opt.het[k1/5], lty = 2, col = "azure3")
  }

  points(y = t.invade[1,3:(dim(t.invade)[2])],x = (1:50)*100-50,t="l",col="blue")
  points(y = t.invade[2,3:(dim(t.invade)[2])],x = (1:50)*100-50,t="l",col="green")
  points(y = t.invade[3,3:(dim(t.invade)[2])],x = (1:50)*100-50,t="l",col="orange")
  points(y = t.invade[4,3:(dim(t.invade)[2])],x = (1:50)*100-50,t="l",col="red")

  points(y = t.exit[1,3:(dim(t.exit)[2])]/3,x = (1:100)*50-25,t="l",col="blue", lty=6)
  points(y = t.exit[2,3:(dim(t.exit)[2])]/3,x = (1:100)*50-25,t="l",col="green", lty=6)
  points(y = t.exit[3,3:(dim(t.exit)[2])]/3,x = (1:100)*50-25,t="l",col="orange", lty=6)
  points(y = t.exit[4,3:(dim(t.exit)[2])]/3,x = (1:100)*50-25,t="l",col="red", lty=6)

  if (i == 1) {
    legend(1600,.23,c("Invading alleles","Exiting alleles"),lty = c(1,6), bty="n",cex = 1.6)
    legend(2600,.29,c(expression(10 ^6),expression(10 ^5),expression(10 ^4),expression(10 ^3)),
           lty = c(1,1,1,1), bty = "n",col = c("red","orange","green","blue"),title =expression(paste("Pop size (",italic("N"),")",sep="")), cex = 1.6)
    legend(3200,.25,c("Susceptible \n to invasion","Optimal hotspots","Both"),bty="n",fill=c("mistyrose","azure","plum2"),border=NA, cex = 1.6)
  }


}








