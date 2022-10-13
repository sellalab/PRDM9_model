### Baker et al. 2022 - Down the Penrose stairs: How selection for fewer recombination hotspots maintains their existence
######################################
### Script for generating Figure 6-1 ###
######################################
setwd("/Users/zbaker/Documents/Columbia Work/Modeling_final/github")
source("simPRDM9_functions.R")

k1 = 20

temp = read.table("source_data/pi.table2.txt")
normal.pi = as.numeric(temp[1,])
Nu.pi = as.numeric(temp[2,])
Nv.pi = as.numeric(temp[3,])
Nuv.pi = as.numeric(temp[4,])

temp = read.table("source_data/t.table.txt")
normal.t = as.numeric(temp[1,])
Nu.t = as.numeric(temp[2,])
Nv.t = as.numeric(temp[3,])
Nuv.t = as.numeric(temp[4,])

par(mfrow=c(1,2))
plot(y = log10(normal.t), x = c(3,3.5,4,4.5,5,5.5,6),
     t = "l", lwd = 3, xlab = expression(paste("Population size (",italic("N"),")",sep="")),
     ylab = expression(paste("Turnover time (log scale)",sep="")),
     axes = "F")


points(y = log10(Nu.t), x = c(3,3.5,4,4.5,5,5.5,6),
       t = "l", lwd = 3, col = "blue")
points(y= log10(Nv.t), x = c(3,3.5,4,4.5,5,5.5,6),
       t = "l", lwd = 3, col = "red")
points(y= log10(Nuv.t), x = c(3,3.5,4,4.5,5,5.5,6),
       t = "l", lwd = 3, col = "purple")
legend(3,3,c(expression(paste("PRDM9 (",nu,")",sep="")),
               expression(paste("binding sites (",mu,")",sep="")),
               expression(paste("both (",mu," and ",nu,")",sep="")),
               "neither"),
       title = "Scaled mutation rate at",
       lty=rep(1,4), lwd = rep(3,4),
       col = c("blue","red","purple","black"))

axis(side=2, at=c(2,3,4,5),las=1, labels = c(expression(10 ^2),
                                             expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5)))
axis(side=1, at=c(3,4,5,6),las=1, labels = c(expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5),
                                             expression(10 ^6)))
box()

plot(y = normal.pi, x = c(3,3.5,4,4.5,5,5.5,6), ylim = c(0,1),
     t="l", lwd=3, xlab = expression(paste("Population size (",italic("N"),")",sep="")),
     ylab = expression(paste("PRDM9 Diversity (",pi,")",sep="")), axes = "F")

points(y = Nu.pi, x = c(3,3.5,4,4.5,5,5.5,6),
       t = "l", lwd = 3, col = "blue")
points(y= Nv.pi, x = c(3,3.5,4,4.5,5,5.5,6),
       t = "l", lwd = 3, col = "red")
points(y= Nuv.pi, x = c(3,3.5,4,4.5,5,5.5,6),
       t = "l", lwd = 3, col = "purple")

axis(side=2, at=c(0,1),las=1)
axis(side=1, at=c(3,4,5,6),las=1, labels = c(expression(10 ^3),
                                             expression(10 ^4),
                                             expression(10 ^5),
                                             expression(10 ^6)))
box()




