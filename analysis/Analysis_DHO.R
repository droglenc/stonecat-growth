################################################################################
## Load libraries
library(FSA)
library(readxl)
library(dplyr)
library(magrittr)
library(lubridate)
library(minpack.lm)
library(nlstools)
library(plotrix)

## Set the random number seed so that the bootstrap results remain contants
set.seed(93834)

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

## Fit traditional VBGF ... used several starting values, two algorithms to
##   make sure that the results were robust.
vbT <- vbFuns()
( nystarts1 <- vbStarts(tl~age,data=nyFish) )
nyfit1 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=nystarts1,algorithm="port")
fitPlot(nyfit1,col.pt=col2rgbt("black",1/10),col.mdl="black")
residPlot(nyfit1)

### try different starting values, same algorithm
nystarts2 <- list(Linf=200,K=0.1,t0=0)
nyfit2 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=nystarts2,algorithm="port")
nystarts3 <- list(Linf=100,K=0.1,t0=-0.5)
nyfit3 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=nystarts3,algorithm="port")
### try different algorithms, same starting values
nyfit4 <- nls(tl~vbT(age,Linf,K,t0),data=nyFish,start=nystarts1)
nyfit5 <- nlsLM(tl~vbT(age,Linf,K,t0),data=nyFish,start=nystarts1)
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(nyfit1),coef(nyfit2),coef(nyfit3),coef(nyfit4),coef(nyfit5)),5)

## Bootstrap first fit
nyboot1 <- nlsBoot(nyfit1)
cbind(Est=coef(nyfit1),confint(nyboot1))

## Bootstrap CIs for increments (used in Figure 3 plot below)
nyLenPred <- apply(nyboot1$coefboot,MARGIN=1,FUN=vbT,t=1:5)
nyLenMean <- apply(nyLenPred,MARGIN=1,FUN=mean)
nyLenCI <- apply(nyLenPred,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
nyIncPred <- apply(nyLenPred,MARGIN=2,FUN=diff)
nyIncMean <- apply(nyIncPred,MARGIN=1,FUN=mean)
nyIncCI <- apply(nyIncPred,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
( nyPredictions <- data.frame(tl=nyLenMean[-5],tlCI=t(nyLenCI[,-5]),
                              inc=nyIncMean,incCI=t(nyIncCI)) )





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
length(unique(vtraw$tag[vtraw$River=="LaPlatte"]))   ### IN MANUSCRIPT
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
xtabs(~River+event,data=vtrecaps)   ### In manuscrpt
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

### remove unneeded objects
rm(tag1,tag2,tag3,tag4,tmp2,tmp3,tmp4,numrecaps)

## Francis version of the tag-recapture VB model with a seasonal component
##  This will predict the mean annual growth rate at two lengths ... L1 and L2
##  Along with the time of peak growth (w) and a level of seasonality (u)
### Choose two values for L1 and L2
Ls1 <- c(100,150,75,175,125,200)
### Create a function with the model in it
vbFT <- function(Lm,T1,T2,L1,L2,g1,g2,w,u) {
  phi <- function(Ti,u,w) u*sin(2*pi*(Ti-w))/(2*pi) 
  ((L2*g1-L1*g2)/(g1-g2)-Lm)*(1-(1+(g1-g2)/(L1-L2))^((T2-T1)+(phi(T2,u,w)-phi(T1,u,w))))
}



## =============================================================================
## -----------------------------------------------------------------------------
##
## Missisquoi River (alone)
##
## -----------------------------------------------------------------------------
## =============================================================================

## Isolate LaPlatte data
MISrecaps <- filterD(vtrecaps2,River=="Missisquoi")
nrow(MISrecaps)

### Delete all observations that were recaptured within 7 days (change in length
###   is likely just measuring error)
### Also computed number of years between marking and recapture
MISrecaps2 <- filterD(MISrecaps,dt*365>7) %>%
  mutate(YAL=year(rDate)-year(mDate))
nrow(MISrecaps2)
### Number of fish deleted
nrow(MISrecaps)-nrow(MISrecaps2)
### Number of years between captures
prop.table(xtabs(~YAL,data=MISrecaps2))*100
### Summarize lengths
Summarize(~mtl,data=MISrecaps2)
Summarize(~rtl,data=MISrecaps2)

### Fit the model with several starting values and different algs to test robustness
MISstarts1 <- list(g1=35,g2=15,w=0.5,u=2)
MISfit1 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=MISrecaps2,
               start=MISstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
residPlot(MISfit1)

MISstarts2 <- list(g1=25,g2=20,w=0.25,u=1)
MISfit2 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=MISrecaps2,
               start=MISstarts2,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
MISstarts3 <- list(g1=45,g2=10,w=0.75,u=1)
MISfit3 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=MISrecaps2,
               start=MISstarts3,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
MISfit4 <- nlsLM(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=MISrecaps2,
                 start=MISstarts1,algorithm="port",
                 lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(MISfit1),coef(MISfit2),coef(MISfit3),coef(MISfit4)),5)

## Bootstrap the first fit
MISboot1 <- nlsBoot(MISfit1)
MIScf1 <- coef(MISfit1)
MISci1 <- confint(MISboot1)
cbind(Est=MIScf1,MISci1)

## Tet intervals for growth increments by using different Ls
MISfit1a <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[3],Ls1[4],g1,g2,w,u),data=MISrecaps2,
                start=MISstarts1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
MISboot1a <- nlsBoot(MISfit1a)
MIScf2 <- coef(MISfit1a)
MISci2 <- confint(MISboot1a)

MISfit1b <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[5],Ls1[6],g1,g2,w,u),data=MISrecaps2,
                start=MISstarts1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
MISboot1b <- nlsBoot(MISfit1b)
MIScf3 <- coef(MISfit1b)
MISci3 <- confint(MISboot1b)

misPredictions <- data.frame(tl=Ls1,
                             inc=c(MIScf1[c("g1","g2")],MIScf2[c("g1","g2")],
                                   MIScf3[c("g1","g2")]),
                             ci=rbind(MISci1[c("g1","g2"),],MISci2[c("g1","g2"),],
                                      MISci3[c("g1","g2"),])) %>%
  arrange(tl)

#### These results appear highly variable, likely due to a small n. For example,
#### the predicted growth rate for a 100 mm fish is from 24 to 40 mm, whereas
#### the same for the LaPlatte R. is from 40 to 45 mm. Additionally, the 
#### decline in predicted growth increaments suggests a Linf in the neighborhood
#### of 300 mm, which seems too large. Finally, the growth incrments CIs overlap
#### for, generally, three 25-mm length classes (where they don't at all for the
#### LaPlatte R.). All-in-all, these data seem too weak to continue with. Also,
#### I would not lump them with the LaPlatte R fish ... why weaken that good data
#### with poorer data or fish that may (or may not) have a different growth pattern.




## =============================================================================
## -----------------------------------------------------------------------------
##
## LaPlatte River (alone)
##
## -----------------------------------------------------------------------------
## =============================================================================

## Isolate LaPlatte data
LPrecaps <- filterD(vtrecaps2,River=="LaPlatte")
nrow(LPrecaps)

### Delete all observations that were recaptured within 7 days (change in length
###   is likely just measuring error)
### Also computed number of years between marking and recapture
LPrecaps2 <- filterD(LPrecaps,dt*365>7) %>%
  mutate(YAL=year(rDate)-year(mDate))
nrow(LPrecaps2)
### Number of fish deleted
nrow(LPrecaps)-nrow(LPrecaps2)
### Number of years between captures
prop.table(xtabs(~YAL,data=LPrecaps2))*100
### Summarize lengths
Summarize(~mtl,data=LPrecaps2)
Summarize(~rtl,data=LPrecaps2)

### Fit the model with several starting values and different algs to test robustness
LPstarts1 <- list(g1=35,g2=15,w=0.5,u=2)
LPfit1 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=LPrecaps2,
                start=LPstarts1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
residPlot(LPfit1)  # a little heteroscedastic

LPstarts2 <- list(g1=25,g2=20,w=0.25,u=1)
LPfit2 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=LPrecaps2,
              start=LPstarts2,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPstarts3 <- list(g1=45,g2=10,w=0.75,u=1)
LPfit3 <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=LPrecaps2,
              start=LPstarts3,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPfit4 <- nlsLM(dtl~vbFT(mtl,mDate12,rDate12,Ls1[1],Ls1[2],g1,g2,w,u),data=LPrecaps2,
                start=LPstarts1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(LPfit1),coef(LPfit2),coef(LPfit3),coef(LPfit4)),5)

## Bootstrap the first fit
LPboot1 <- nlsBoot(LPfit1)
LPcf1 <- coef(LPfit1)
LPci1 <- confint(LPboot1)
cbind(Est=LPcf1,LPci1)

## Tet intervals for growth increments by using different Ls
LPfit1a <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[3],Ls1[4],g1,g2,w,u),data=LPrecaps2,
               start=LPstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPboot1a <- nlsBoot(LPfit1a)
LPcf2 <- coef(LPfit1a)
LPci2 <- confint(LPboot1a)

LPfit1b <- nls(dtl~vbFT(mtl,mDate12,rDate12,Ls1[5],Ls1[6],g1,g2,w,u),data=LPrecaps2,
               start=LPstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPboot1b <- nlsBoot(LPfit1b)
LPcf3 <- coef(LPfit1b)
LPci3 <- confint(LPboot1b)

LPPredictions <- data.frame(tl=Ls1,
                            inc=c(LPcf1[c("g1","g2")],LPcf2[c("g1","g2")],
                                  LPcf3[c("g1","g2")]),
                            ci=rbind(LPci1[c("g1","g2"),],LPci2[c("g1","g2"),],
                                     LPci3[c("g1","g2"),])) %>%
  arrange(tl)



################################################################################
## =============================================================================
##
## Manuscript Figures
##
## =============================================================================
################################################################################

## Fitted Line Plot ... IN THE MANUSCRIPT (Figure 1)
windows(5,5); par(mar=c(3,3,0.7,0.7),mgp=c(1.9,0.5,0),tcl=-0.2)
plot(tl~age,data=nyFish,pch=19,col=col2rgbt("black",1/5),
     xlim=c(0,6),ylim=c(0,200),
     xlab="Adjusted Consensus Spine Age",ylab="Total Length (mm)")
curve(vbT(x,coef(nyfit1)),from=0,to=6,add=TRUE)
x <- seq(0,6,length.out=99)
nyLenPred <- apply(nyboot1$coefboot,MARGIN=1,FUN=vbT,t=x)
nyLenCI <- apply(nyLenPred,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
lines(x,nyLenCI["2.5%",],lty=2)
lines(x,nyLenCI["97.5%",],lty=2)

### Histogram of times-at-large ... IN THE MANUSCRIPT (Figure 2)
hist(~dt,data=LPrecaps2,w=14/365,xlim=c(0,2),ylim=c(0,30),
     xlab="Time-at-Large (years)",ylab="Frequency of Capture-Recapture Events")


## Plot growth increment data across all populations ... IN THE MANUSCRIPT (Figure 3)
## Compute a SL to TL conversion from NY data
lmSL2TL <- lm(tl~sl,data=nyFish)
residPlot(lmSL2TL)  ## 3 possible outliers
fitPlot(lmSL2TL)
coef(lmSL2TL)
rSquared(lmSL2TL)

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


## Make the plot
windows(6,5); par(mar=c(3,3,0.7,0.7),mgp=c(1.9,0.5,0),tcl=-0.2)
plot(NULL,xlab="Initial Total Length (mm)",ylab="Annual Total Length Increment (mm)",
     xlim=c(-20,260),ylim=c(0,65))
with(LPPredictions,plotCI(tl,inc,li=ci.95..LCI,ui=ci.95..UCI,add=TRUE,pch=19))
lines(inc~tl,data=LPPredictions,col="black",lwd=2)
arrows(53,48,LPPredictions$tl[1]-4,LPPredictions$inc[1]+1,length=0.1)
text(57,48,"LaPalette R. (VT)",pos=2)
with(nyPredictions,plotCI(tl,inc,li=tlCI.2.5.,ui=tlCI.97.5.,err="x",add=TRUE,col="gray40"))
with(nyPredictions,plotCI(tl,inc,li=incCI.2.5.,ui=incCI.97.5.,add=TRUE,col="gray40",pch=19))
lines(inc~tl,data=nyPredictions,col="gray40",lwd=2)
arrows(66,53,nyPredictions$tl[1]-4,nyPredictions$inc[1]+1,length=0.1)
text(70,53,"Great Chazy R. (NY)",pos=2)
lines(incbc~bctl,data=Carlson66,lty=2)
with(Carlson66,text(bctl[1],incbc[1],"Vermillion R. (SD)",pos=2))
lines(incbc~bctl,data=Paruch79,lty=2)
with(Paruch79,text(bctl[2],incbc[2],"Wisconsin",pos=2))
lines(inc~mntl,data=Gilbert53s,lty=2)
with(Gilbert53s,text(mntl[1],inc[1],"OH Streams",pos=2))
lines(inc~mntl,data=Gilbert53le,lty=2)
with(Gilbert53le,text(mntl[1],inc[1],"Lake Erie (OH)",pos=2))






## Compute a TL to W conversion from NY data
nyFish1 <- mutate(nyFish,logw=log(wt),logtl=log(tl))
nyFish1 <- nyFish1[-c(66,172),] # 2 big outliers
lmTL2W <- lm(logw~logtl,data=nyFish1)
residPlot(lmTL2W)
fitPlot(lmTL2W)
### Curved, likely due to different years and times of year
### DON'T USE
