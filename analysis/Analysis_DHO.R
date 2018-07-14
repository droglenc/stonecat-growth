################################################################################
## Load libraries
library(FSA)
library(readxl)
library(dplyr)
library(magrittr)
library(lubridate)
library(minpack.lm)
library(nlstools)
library(ggplot2)
library(plotrix)


################################################################################
## =============================================================================
##
## New York
##
## =============================================================================
################################################################################

## Load NY data and reduce to just those variables of interest
## Note that the observed annuli were modified to create fractional ages based
##   on a presumed growing season from June 1st to November 1st.

nyAges <- read_xlsx("data/NY.xlsx",sheet="Ages")

nyFish <- read_xlsx("data/NY.xlsx",sheet="Fish") %>%
  rename(date=`Date collected`,tl=`Total Length`,
         sl=`Standard Length`,wt=`Weight`) %>%
  mutate(mon=month(date,label=TRUE),yr=year(date),day=yday(date)) %>%
  select(-Location,-Sex,-`Excised Date`,-date) %>%
  left_join(nyAges[,c("ID","FinalAge")],by="ID") %>%
  rename(annuli=FinalAge) %>%
  mutate(age= case_when(
    day<152 ~ annuli,                        # Before June 1st
    day<305 ~ annuli + (day-152)/(305-152),  # Btwn June 1st and Nov 1st
    TRUE ~ annuli+1))                        # After Nov 1st
headtail(nyFish)

## Examine sample sizes by month of capture and annuli
addmargins(xtabs(~mon+annuli,data=nyFish))
## Summarize lengths and weights
Summarize(~tl,data=nyFish)
Summarize(~wt,data=nyFish)
## Examine length frequency histogram ... WEIRD ... non-random sampling??
hist(~tl,data=nyFish,w=10)
## Examine "growth trajectory" plot
plot(tl~age,data=nyFish,pch=19,col=col2rgbt("black",1/10))

## Fit traditional VBGF ... used several starting values, two algorithms to
##   make sure that the results were robust.
vbT <- vbFuns()
( vbTs1 <- vbStarts(tl~age,data=nyFish) )
vbTf1 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=vbTs1,algorithm="port")
### try different starting values, same algorithm
vbTs2 <- list(Linf=200,K=0.1,t0=0)
vbTf2 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=vbTs2,algorithm="port")
vbTs3 <- list(Linf=100,K=0.1,t0=-0.5)
vbTf3 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=vbTs3,algorithm="port")
### try different algorithms, same starting values
vbTf4 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=vbTs1)
vbTf5 <- nlsLM(tl~vbT(age,Linf,K,t0),data=nyFish,start=vbTs1)
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(vbTf1),coef(vbTf2),coef(vbTf3),coef(vbTf4),coef(vbTf5)),5)
### Quick-and-dirty visual
fitPlot(vbTf1,col.pt=col2rgbt("black",1/10),col.mdl="black")
residPlot(vbTf1)

## Bootstrap first fit
bootT1 <- nlsBoot(vbTf1)
cbind(Est=coef(vbTf1),confint(bootT1))

## Bootstrap CIs for increments
pLenAtAge <- apply(bootT1$coefboot,MARGIN=1,FUN=vbT,t=1:5)
pIncAtAge <- apply(pLenAtAge,MARGIN=2,FUN=diff)
mnLenAtAge <- apply(pLenAtAge,MARGIN=1,FUN=mean)
ciLenAtAge <- apply(pLenAtAge,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
mnIncAtAge <- apply(pIncAtAge,MARGIN=1,FUN=mean)
ciIncAtAge <- apply(pIncAtAge,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
nyPredictions <- data.frame(tl=mnLenAtAge[-5],tlCI=t(ciLenAtAge[,-5]),
                            inc=mnIncAtAge,incCI=t(ciIncAtAge))


windows(5,5); par(mar=c(3,3,0.7,0.7),mgp=c(1.9,0.5,0),tcl=-0.2)
plot(tl~age,data=nyFish,pch=19,col=col2rgbt("black",1/10),
     xlim=c(0,6),ylim=c(0,200),
     xlab="Consensus Spine Age",ylab="Total Length (mm)")
curve(vbT(x,coef(vbTf1)),from=0,to=6,add=TRUE)
x <- seq(0,6,length.out=99)
pLenAtAge <- apply(bootT1$coefboot,MARGIN=1,FUN=vbT,t=x)
ciLenAtAge <- apply(pLenAtAge,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
lines(x,ciLenAtAge["2.5%",],lty=2)
lines(x,ciLenAtAge["97.5%",],lty=2)

################################################################################
## =============================================================================
##
## Vermont
##
## =============================================================================
################################################################################

## Load VT data and reduce to just those variables of interest
## Removed one observation where the length was an obvious outlier
## Removed some bad tag names
## Removed three fish that were recaptured on same day (??)
vtraw <- read_xlsx("data/VT.xlsx") %>%
  rename(tag=`PIT #/Unique ID`,tl=`Length (mm)`) %>%
  select(-Site,-`Fish #`,-`Sampling Type`,-Temp,-`Weight (g)`,
         -Conductvity,-`VIE Pattern`,-`Mark/Recap`,-NOTES) %>%
  filterD(!(tag=="599C1E0" & tl==174)) %>%
  filterD(!tag %in% c("NA","No Tag","NT","RELS001")) %>%
  filterD(!tag %in% c("595F1(?)","595F1AE","65E7B32"))

## Total tagged by river
length(unique(vtraw$tag[vtraw$River=="LaPlatte"]))
length(unique(vtraw$tag[vtraw$River=="Missisquoi"]))

## Get a data.frame of just recaptured fish
### Find tag numbers for recaptured fish
recap.tags <- xtabs(~tag,data=vtraw)
recap.tags <- names(recap.tags)[recap.tags>1]  # tags that were recaptured
### Reduce data.frame to just fish with recaptured tag numbers
vtrecaps <- filterD(vtraw,tag %in% recap.tags) %>%
  arrange(tag,Date)
### Add an "event" variable (1=first seen, 2=next seen, etc.)
numrecaps <- xtabs(~tag,data=vtrecaps)
names(numrecaps) <- NULL
vtrecaps %<>% mutate(event=unlist(apply(numrecaps,1,seq_len))) %>%
  select(tag,River,event,Date,tl)

## Separate into recapture times ... Note that 1 is always marking, >1 is always
##   a recapture. If a fish is captured more than once, then >1 may be a marking
##   for a later recapture (i.e., event=2 can be a marking for event=3)
xtabs(~River+event,data=vtrecaps)
tag1 <- filterD(vtrecaps,event==1)
tag2 <- filterD(vtrecaps,event==2)
tag3 <- filterD(vtrecaps,event==3)
tag4 <- filterD(vtrecaps,event==4)

### Joins in order to compute change in langth and time-at-liberty
tmp4 <- right_join(tag3,tag4,by="tag")
tmp3 <- right_join(tag2,tag3,by="tag")
tmp2 <- right_join(tag1,tag2,by="tag")
### Put them all together, give some better variable names, change the dates
###   to fractions of a year since Jan 1, 2012, compute the time (in yrs) at
###   large, compute the change in length, and order the data.
vtrecaps2 <- rbind(tmp2,tmp3,tmp4) %>%
  rename(River=River.x,mDate=Date.x,rDate=Date.y,mtl=tl.x,rtl=tl.y) %>%
  mutate(mDate12=as.double(difftime(mDate,as.POSIXct("2012-01-01")))/365,
         rDate12=as.double(difftime(rDate,as.POSIXct("2012-01-01")))/365,
         dt=rDate12-mDate12,dtl=rtl-mtl) %>%
  select(-event.x,-event.y,-River.y) %>%
  arrange(River,tag,mtl) %>%
  as.data.frame()
### Delete all observations that were recaptured within 7 days (change in length
###   is likely just measuring error)
### Also computed number of years between marking and recapture
vtrecaps3 <- filterD(vtrecaps2,dt*365>7) %>%
  mutate(YAL=year(rDate)-year(mDate))
### Number of fish deleted
nrow(vtrecaps2)-nrow(vtrecaps3)

### remove unneeded objects
rm(tag1,tag2,tag3,tag4,tmp2,tmp3,tmp4,numrecaps)

## Quick explorations
### Few fish from Missisquoi ... thus don't separate
xtabs(~River,data=vtrecaps3)
### Number of years between captures
prop.table(xtabs(~YAL,data=vtrecaps3))*100
### Summarize lengths
Summarize(~mtl,data=vtrecaps3)
Summarize(~rtl,data=vtrecaps3)

### Look at captures over time by tag
###   ordered by initial capture data
###   51 fish at a time
ggplot(data=vtrecaps[vtrecaps$tag %in% recap.tags[1:51],],aes(x=Date,y=tag)) +
  geom_line(size=0.5) +
  geom_point(size=0.75,color="red")
ggplot(data=vtrecaps[vtrecaps$tag %in% recap.tags[52:103],],aes(x=Date,y=tag)) +
  geom_line(size=0.5) +
  geom_point(size=0.75,color="red")
ggplot(data=vtrecaps[vtrecaps$tag %in% recap.tags[104:153],],aes(x=Date,y=tag)) +
  geom_line(size=0.5) +
  geom_point(size=0.75,color="red")

### Histogram of times-at-large
windows(5,5); par(mar=c(3,3,0.7,0.7),mgp=c(1.9,0.5,0),tcl=-0.2)
hist(~dt,data=vtrecaps3,w=7/365,xlim=c(0,2),
     xlab="Time-at-Large (years)",ylab="Frequency of Recapture Events")

## Fit the Francis version of the tag-recapture VB model with a seasonal component
##  This will predict the mean annual growth rate at two lengths ... L1 and L2
##  Along with the time of peak growth (w) and a level of seasonality (u)

### Choose two values for L1 and L2
Ls1 <- c(100,150)
### Create a function with the model in it
vbFT <- function(Lm,T1,T2,L1,L2,g1,g2,w,u) {
  phi <- function(Ti,u,w) u*sin(2*pi*(Ti-w))/(2*pi) 
  ((L2*g1-L1*g2)/(g1-g2)-Lm)*(1-(1+(g1-g2)/(L1-L2))^((T2-T1)+(phi(T2,u,w)-phi(T1,u,w))))
}
### Fit the model with several starting values and different algs to test robustness
vbFTs1 <- list(g1=40,g2=10,w=0.5,u=0.5)
vbFTf1 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=vtrecaps3,
              start=vbFTs1,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
residPlot(vbFTf1)  # a little heteroscedastic

vbFTs2 <- list(g1=25,g2=20,w=0.25,u=0.5)
vbFTf2 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=vtrecaps3,
              start=vbFTs2,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
vbFTs3 <- list(g1=25,g2=20,w=0.75,u=0.1)
vbFTf3 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=vtrecaps3,
              start=vbFTs3,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
vbFTf4 <- nlsLM(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=vtrecaps3,
                start=vbFTs1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(vbFTf1),coef(vbFTf2),coef(vbFTf3),coef(vbFTf4)),5)

## Bootstrap the first fit
bootFT1 <- nlsBoot(vbFTf1)
cf1 <- coef(vbFTf1)
ci1 <- confint(bootFT1)
cbind(Est=cf1,ci1)

## Choose different L1s and L2s so as to get intervals for growth increments
Ls2 <- c(75,175)
vbFTf1a <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls2[1],Ls2[2],g1,g2,w,u),data=vtrecaps3,
               start=vbFTs1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
bootFT1a <- nlsBoot(vbFTf1a)
cf2 <- coef(vbFTf1a)
ci2 <- confint(bootFT1a)

Ls3 <- c(125,200)
vbFTf1b <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls3[1],Ls3[2],g1,g2,w,u),data=vtrecaps3,
               start=vbFTs1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
bootFT1b <- nlsBoot(vbFTf1b)
cf3 <- coef(vbFTf1b)
ci3 <- confint(bootFT1b)

vtPredictions <- data.frame(tl=c(Ls1,Ls2,Ls3),
                            inc=c(cf1[c("g1","g2")],cf2[c("g1","g2")],
                                      cf3[c("g1","g2")]),
                            ci=rbind(ci1[c("g1","g2"),],ci2[c("g1","g2"),],
                                     ci3[c("g1","g2"),])) %>%
  arrange(tl)



################################################################################
## =============================================================================
##
## Plot growth increment data across all populations
##
## =============================================================================
################################################################################

## Compute a SL to TL conversion from NY data
lmSL2TL <- lm(tl~sl,data=nyFish)
residPlot(lmSL2TL)  ## 3 possible outliers
fitPlot(lmSL2TL)
coef(lmSL2TL)
rSquared(lmSL2TL)

## Compute a TL to W conversion from NY data
nyFish1 <- mutate(nyFish,logw=log(wt),logtl=log(tl))
nyFish1 <- nyFish1[-c(66,172),] # 2 big outliers
lmTL2W <- lm(logw~logtl,data=nyFish1)
residPlot(lmTL2W)
fitPlot(lmTL2W)
### Curved, likely due to different years and times of year
### DON'T USE

## Results from historical studies
## Carlson (1966) Lake Vermillion, SD ... Total Length at formation, vertebrae
Carlson66 <- data.frame(age=1:7,
                        n=c(4,6,22,6,2,6,1),
                        mntl=c(78.5,96.8,113.8,137.6,155.0,175.6,193.0),
                        bctl=c(68.8,99.7,117.0,136.9,157.8,171.5,179.8)) %>%
  mutate(inc=c(diff(mntl),NA),incbc=c(diff(bctl),NA))
## Paruch (1979) Wisconsin ... Total Length at capture, pectoral spines
Paruch79 <- data.frame(age=0:5,
                        n=c(27,20,19,5,1,2),
                        mntl=c(46,102,148,159,199,187),
                        bctl=c(NA,51,95,124,152,162)) %>%
  mutate(inc=c(diff(mntl),NA),incbc=c(diff(bctl),NA))
## Gilbert (1953) Streams ... Standard Length, vertebrae
Gilbert53s <- data.frame(age=1:6,
                         n=NA,
                         mnsl=c(54,73,89,104,116,129)) %>%
  mutate(mntl=predict(lmSL2TL,data.frame(sl=mnsl)),
         inc=c(diff(mntl),NA))
## Gilbert (1953) Lake Erie ... Standard Length, vertebrae
Gilbert53le <- data.frame(age=1:9,
                          n=NA,
                          mnsl=c(68,121,162,181,195,203,208,224,237)) %>%
  mutate(mntl=predict(lmSL2TL,data.frame(sl=mnsl)),
         inc=c(diff(mntl),NA))



windows(5,5); par(mar=c(3,3,0.7,0.7),mgp=c(1.9,0.5,0),tcl=-0.2)
plot(NULL,xlab="Initial TL (mm)",ylab="Annual TL Increment (mm)",
     xlim=c(0,260),ylim=c(0,65))
with(vtPredictions,plotCI(tl,inc,li=ci.95..LCI,ui=ci.95..UCI,add=TRUE,pch=19))
lines(inc~tl,data=vtPredictions,col="black",lwd=2)
text(0,50,"Vermont",pos=4)
arrows(51,50,vtPredictions$tl[1]-4,vtPredictions$inc[1]+1,length=0.1)
with(nyPredictions,plotCI(tl,inc,li=tlCI.2.5.,ui=tlCI.97.5.,err="x",add=TRUE,col="gray60"))
with(nyPredictions,plotCI(tl,inc,li=incCI.2.5.,ui=incCI.97.5.,add=TRUE,col="gray60",pch=19))
lines(inc~tl,data=nyPredictions,col="gray60",lwd=2)
text(0,55,"New York",pos=4)
arrows(53,55,nyPredictions$tl[1]-4,nyPredictions$inc[1]+1,length=0.1)
lines(incbc~bctl,data=Carlson66,lty=2)
with(Carlson66,text(bctl[1],incbc[1],"Vermillion (SD)",pos=2))
lines(incbc~bctl,data=Paruch79,lty=2)
with(Paruch79,text(bctl[2],incbc[2],"Wisconsin",pos=2))
lines(inc~mntl,data=Gilbert53s,lty=2)
with(Gilbert53s,text(mntl[1],inc[1],"Ohio",pos=2))
lines(inc~mntl,data=Gilbert53le,lty=2)
with(Gilbert53le,text(mntl[1],inc[1],"Lake Erie",pos=2))

