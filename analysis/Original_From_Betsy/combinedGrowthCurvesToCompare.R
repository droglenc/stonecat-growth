#combining two methods to estimate growth curve and visually comparing them 
# First I will add the Fabens method using recap data then Von Beralanffy
#then i will graph all of the estimates to compare them. 


#Here I use the Fabens method to estimate a growth curve using the mark recapture data.
#I create 90% confidence intervals on the growth curve using a bootstrapping method (1000 iterations)
#I then overlay the estimated ages obtained from reading spines, to see how they compare
#The two different methods line up rather well even though there is a lot of length area covered
# Please make sure that the location of the needed csv file matches where you saved it. 
library(FSA)
library(nlstools)
xlbl <- "Age (yrs)"
ylbl <- "Total Length (mm)"
clr <- rgb(0,0,0,1/5)
setwd("C:/Users/betsy/Desktop/Betsy/Stonecats_Final/R_Stonecats")
d <- read.csv("GrowthCurve_Fabens.csv", header=TRUE)
#added a column for change in length in yrs
d$changeYrs <- d$ChangeTime/365.25
#creating a starting list for Least Sums Squares
svb1 <- list('Linf'=150, 'K'=0.6)

#fitting the model and getting the results
fit1 <- nls(ChangeLength~(Linf-StartingLength)*(1-exp(-K*changeYrs)),data=d,start=svb1)
summary(fit1)
confint(fit1)
(cf <- coef(fit1))

#now its time to bootstrap the data
boot1 <- nlsBoot(fit1,niter=1000)
est1 <- boot1$coefboot
est1[1:5,]

# Confidence level is currently set at .9
cfboot <- confint(boot1,level=0.95)


#now for the VonBert.
#Adding the spine age data to have the recap data to graph over it
setwd("C:/Users/betsy/Desktop/Betsy/Stonecats_Final/R_Stonecats")
age <- read.csv("stonecatages.csv", header=TRUE)
age$MocFinalAge <- age$FinalAge+.5
library(nlstools)
#getting the starting points for von Bertalanffy
(svb2 <- vbStarts(length~FinalAge,data=age,type="typical"))
fit2 <- nls(length~Linf*(1-exp(-K*(FinalAge-t0))),data=age,start=svb2)
summary(fit2)
confint(fit2)
(cf2 <- coef(fit2))

#now its time to bootstrap the data
boot2 <- nlsBoot(fit2,niter=1000)
est2 <- boot2$coefboot
est2[1:5,]
#just looking at what the bootstapped values look like
cfboot2 <- confint(boot2,level=0.95)

#graphing everything at the same time.
#making the plot given the bootstrapping results. 
plot(length~MocFinalAge,data=age,xlab=xlbl,ylab=ylbl,font.lab=10,pch=16,col=clr,xlim=c(0,5.5),ylim=c(0,225))

#Black represent the recap data infomraiton 
curve(cf["Linf"]*(1-exp(-cf["K"]*(x-0))),from=0,to=5,n=500,lwd=4,col="black",lty=1,add=TRUE)
#adding lower CI
curve(cfboot[1,1]*(1-exp(-cfboot[2,1]*(x-0))),from=0,to=5,n=500,lwd=3,col="black",lty=2,add=TRUE)
#adding upper CI
curve(cfboot[1,2]*(1-exp(-cfboot[2,2]*(x-0))),from=0,to=5,n=500,lwd=3,col="black",lty=2,add=TRUE)

#Grey represent the outputs from the spien information 
curve(cf2["Linf"]*(1-exp(-cf2["K"]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=4,col="grey",lty=1,add=TRUE)
#adding lower CI
curve(cfboot2[1,1]*(1-exp(-cfboot2[2,1]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=3,col="grey",lty=2,add=TRUE)
#adding upper CI
curve(cfboot2[1,2]*(1-exp(-cfboot2[2,2]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=3,col="grey",lty=2,add=TRUE)


#Graphing the two plots but using T0=0 for the VonBertalanffy estimation to compare graphs 
plot(length~MocFinalAge,data=age,xlab=xlbl,ylab=ylbl,pch=16,col=clr,xlim=c(0,5.5),ylim=c(0,225))

#black represent the recap data infomraiton 
curve(cf["Linf"]*(1-exp(-cf["K"]*(x-0))),from=0,to=5,n=500,lwd=4,col="black",lty=1,add=TRUE)
#adding lower CI
curve(cfboot[1,1]*(1-exp(-cfboot[2,1]*(x-0))),from=0,to=5,n=500,lwd=3,col="black",lty=2,add=TRUE)
#adding upper CI
curve(cfboot[1,2]*(1-exp(-cfboot[2,2]*(x-0))),from=0,to=5,n=500,lwd=3,col="black",lty=2,add=TRUE)

#Grey represent the outputs from the spien information 
curve(cf2["Linf"]*(1-exp(-cf2["K"]*(x-0))),from=0,to=5,n=500,lwd=4,col="grey",lty=1,add=TRUE)
#adding lower 90
curve(cfboot2[1,1]*(1-exp(-cfboot2[2,1]*(x-0))),from=0,to=5,n=500,lwd=3,col="grey",lty=2,add=TRUE)
#adding upper 90
curve(cfboot2[1,2]*(1-exp(-cfboot2[2,2]*(x-0))),from=0,to=5,n=500,lwd=3,col="grey",lty=2,add=TRUE)

#using the t0 estimated from vonBert spines data imported into the estimates for recap data using Fabens
plot(length~MocFinalAge,data=age,xlab=xlbl,ylab=ylbl,pch=16,col=clr,xlim=c(0,5.5),ylim=c(0,225))

#black represent the recap data infomraiton 
curve(cf["Linf"]*(1-exp(-cf["K"]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=4,col="Black",lty=1,add=TRUE)
#adding lower CI
curve(cfboot[1,1]*(1-exp(-cfboot[2,1]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=3,col="black",lty=2,add=TRUE)
#adding upper CI
curve(cfboot[1,2]*(1-exp(-cfboot[2,2]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=3,col="black",lty=2,add=TRUE)

#Grey represent the outputs from the spien information 
curve(cf2["Linf"]*(1-exp(-cf2["K"]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=4,col="grey",lty=1,add=TRUE)
#adding lower CI
curve(cfboot2[1,1]*(1-exp(-cfboot2[2,1]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=3,col="grey",lty=2,add=TRUE)
#adding upper CI
curve(cfboot2[1,2]*(1-exp(-cfboot2[2,2]*(x-cf2["t0"]))),from=0,to=5,n=500,lwd=3,col="grey",lty=2,add=TRUE)
