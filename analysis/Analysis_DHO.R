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
library(ggplot2)

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
##   on a presumed growing season from June 1st (152) to November 1st (305).
gstart <- 152 ## June 1st
gend <- 305   ## Nov 1st

nyAges <- read_xlsx("data/NY.xlsx",sheet="Ages")
nyFish <- read_xlsx("data/NY.xlsx",sheet="Fish") %>%
  rename(date=`Date collected`,tl=`Total Length`,
         sl=`Standard Length`,wt=`Weight`) %>%
  mutate(mon=month(date,label=TRUE),yr=year(date),day=yday(date)) %>%
  mutate(tl2=tl/(1-0.024)) %>%  # 2.4% freezing effect on TL (avg from in Ogle)
  select(-Location,-Sex,-`Excised Date`,-date) %>%
  left_join(nyAges[,c("ID","FinalAge")],by="ID") %>%
  rename(annuli=FinalAge) %>%
  mutate(age= case_when(
    day<gstart ~ annuli,                             # Before gstart
    day<gend ~ annuli + (day-gstart)/(gend-gstart),  # Btwn gstart and gend
    TRUE ~ annuli+1))                                # After gend
headtail(nyFish)

## Examine sample sizes by month of capture and annuli
addmargins(xtabs(~mon+annuli,data=nyFish))                     ### IN MANUSCRIPT
## Summarize lengths and weights
Summarize(~tl,data=nyFish,digits=1)                            ### IN MANUSCRIPT
## Examine length frequency histogram ... WEIRD ... non-random sampling??
hist(~tl,data=nyFish,w=10)


## Compute a SL to TL conversion from NY data
lmSL2TL <- lm(tl~sl,data=nyFish)
residPlot(lmSL2TL)  ## 3 possible outliers
fitPlot(lmSL2TL)
round(coef(lmSL2TL),3)                                         ### IN MANUSCRIPT
rSquared(lmSL2TL,digits=3)                                     ### IN MANUSCRIPT


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
### Simple prediction ... used when testing robustness of seasonal dates
predict(nyfit1,data.frame(age=2:4))

## Bootstrap first fit
nyboot1 <- nlsBoot(nyfit1)
round(cbind(Est=coef(nyfit1),confint(nyboot1)),2)              ### IN MANUSCRIPT




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
length(unique(vtraw$tag[vtraw$River=="LaPlatte"]))             ### IN MANUSCRIPT
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
xtabs(~River+event,data=vtrecaps)                              ### IN MANUSCRIPT
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
Ls1 <- c(100,150,75)
### Create a function with the model in it
vbFT <- vbFuns("Francis3")



## =============================================================================
## -----------------------------------------------------------------------------
##
## Missisquoi River (alone) ... NOT USED IN THE MANUSCRIPT SEE COMMENT AT END
##
## -----------------------------------------------------------------------------
## =============================================================================

## Isolate Missisquoi data
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
MISfit1 <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=MISrecaps2,
               start=MISstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
residPlot(MISfit1)

MISstarts2 <- list(g1=25,g2=20,w=0.25,u=1)
MISfit2 <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=MISrecaps2,
               start=MISstarts2,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
MISstarts3 <- list(g1=45,g2=10,w=0.75,u=1)
MISfit3 <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=MISrecaps2,
               start=MISstarts3,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
MISfit4 <- nlsLM(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=MISrecaps2,
                 start=MISstarts1,algorithm="port",
                 lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(MISfit1),coef(MISfit2),coef(MISfit3),coef(MISfit4)),5)

## Bootstrap the first fit
MISboot1 <- nlsBoot(MISfit1)
MIScf1 <- coef(MISfit1)
MISci1 <- confint(MISboot1)
cbind(Est=MIScf1,MISci1)

#### These results appear highly variable, likely due to a small n. For example,
#### the predicted growth rate for a 100 mm fish is from 24 to 40 mm, whereas
#### the same for the LaPlatte R. is from 40 to 45 mm. Results not shown here,
#### but previously calcuated ... the decline in predicted growth increaments
#### suggests a Linf in the neighborhood of 300 mm, which seems too large.
#### Finally, the growth incrments CIs overlap for, generally, three 25-mm
#### length classes (where they don't at all for the LaPlatte R.). All-in-all,
#### these data seem too weak to continue with. Also, I would not lump them with
#### the LaPlatte R fish ... why weaken that good data with poorer data or fish
#### that may (or may not) have a different growth pattern.




## =============================================================================
## -----------------------------------------------------------------------------
##
## LaPlatte River (alone)
##
## -----------------------------------------------------------------------------
## =============================================================================

## Isolate raw LaPlatte data
LPraw <- filterD(vtraw,River=="LaPlatte") %>%
  mutate(yr=factor(year(Date)),mon=month(Date,label=TRUE,abbr=FALSE))

Summarize(~tl,data=LPraw,digits=1)                             ### IN MANUSCRIPT

## Isolate recapture LaPlatte data
LPrecaps <- filterD(vtrecaps2,River=="LaPlatte")
nrow(LPrecaps)

### Delete all observations that were recaptured within 7 days (change in length
###   is likely just measuring error)
### Also computed number of years between marking and recapture
LPrecaps2 <- filterD(LPrecaps,dt*365>7) %>%
  mutate(YAL=year(rDate)-year(mDate))
nrow(LPrecaps2)
### Number of fish deleted
nrow(LPrecaps)-nrow(LPrecaps2)                                 ### IN MANUSCRIPT
### Number of years between captures
round(prop.table(xtabs(~YAL,data=LPrecaps2))*100,0)            ### IN MANUSCRIPT
### Summarize lengths
Summarize(~mtl,data=LPrecaps2,digits=1)                        ### IN MANUSCRIPT
Summarize(~rtl,data=LPrecaps2,digits=1)                        ### IN MANUSCRIPT

### Fit the model with several starting values and different algs to test robustness
LPstarts1 <- list(g1=35,g2=15,w=0.5,u=2)
LPfit1 <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=LPrecaps2,
                start=LPstarts1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
residPlot(LPfit1)  # a little heteroscedastic

LPstarts2 <- list(g1=25,g2=20,w=0.25,u=1)
LPfit2 <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=LPrecaps2,
              start=LPstarts2,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPstarts3 <- list(g1=45,g2=10,w=0.75,u=1)
LPfit3 <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=LPrecaps2,
              start=LPstarts3,algorithm="port",
              lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPfit4 <- nlsLM(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,Ls1[1],Ls1[2]),data=LPrecaps2,
                start=LPstarts1,algorithm="port",
                lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
### See if the coefficients are similar ... THEY ARE!!!
round(cbind(coef(LPfit1),coef(LPfit2),coef(LPfit3),coef(LPfit4)),5)

## Bootstrap the first fit
LPboot1 <- nlsBoot(LPfit1)
LPcf1 <- coef(LPfit1)
LPci1 <- confint(LPboot1)
round(cbind(Est=LPcf1,LPci1),2)                                ### IN MANUSCRIPT



################################################################################
## =============================================================================
##
## Manuscript Figures
##
## =============================================================================
################################################################################

## Setup a theme
theme_stonecat <- function (base_size=12,base_family="") {
  theme_bw(base_size=base_size,base_family=base_family) +
    theme(panel.grid.major=element_line(colour="gray90",linetype="dashed",size=0.25),
          panel.grid.minor=element_line(colour="gray90",linetype="dashed",size=0.20),
          panel.border = element_rect(color="black",size=0.3),
          strip.background = element_rect(color="black",size=0.3),
          strip.text=element_text(size=8,color="black"),
          axis.text=element_text(size=8,color="black"),
          axis.title=element_text(size=10,color="black"),
          axis.title.x=element_text(margin=margin(t=3,r=0,b=0,l=0),color="black"),
          axis.title.y=element_text(margin=margin(t=0,r=6,b=0,l=0),color="black"),
          axis.line=element_line(size=0.3),
          axis.ticks=element_line(size=0.3))
}

## Setup sizes
singlewide <- 3.50
doublewide <- 7.25
maxheight <- 7.5


### Length Frequency of all fish from LaPlatte R. (Figure 1)
LPraw <- filterD(vtraw,River=="LaPlatte") %>%
  mutate(yr=factor(year(Date)),mon=month(Date,label=TRUE,abbr=FALSE))
LPraw56789 <- filterD(LPraw,mon %in% c("May","June","July","August","September"))

lf <- ggplot(LPraw56789,aes(x=tl)) +
  geom_histogram(position="identity",binwidth=5,
                 fill="gray70",color="black",size=0.2) +
  scale_x_continuous(name="Total length (mm)",
                     limits=c(50,210)) +
  scale_y_continuous(name="Number of stonecats",
                     limits=c(0,40),expand=c(0,0),
                     breaks=seq(0,37,10)) +
  facet_grid(mon~yr)
lf + theme_stonecat()
ggsave("doc/figures/Figure1.TIFF",width=doublewide,height=maxheight)



### Histogram of times-at-large ... IN THE MANUSCRIPT (Figure 2)
tal <- ggplot(LPrecaps2,aes(x=dt)) +
  geom_histogram(position="identity",binwidth=14/365,
                 fill="gray70",color="black",size=0.2) +
  scale_x_continuous(name="Time-at-large (years)",
                     limits=c(0,2.05),expand=c(0,0)) +
  scale_y_continuous(name="Capture-recapture event frequency",
                     limits=c(0,25),expand=c(0,0))
tal + theme_stonecat()
ggsave("doc/figures/Figure2.TIFF",width=singlewide,height=singlewide)


## Fitted Line Plot ... IN THE MANUSCRIPT (Figure 3)
x <- seq(0,6,length.out=999)
nyLenPred <- predict(nyfit1,data.frame(age=x))
nyLenCI <- apply(nyboot1$coefboot,MARGIN=1,FUN=vbT,t=x)
nyLenCI <- apply(nyLenCI,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
nyLenCI <- data.frame(x,nyLenPred,t(nyLenCI))
names(nyLenCI) <- c("age","Est","LCI","UCI")

flp <- ggplot(nyFish,aes(x=age,y=tl)) +
  geom_point(alpha=1/5) +
  scale_x_continuous(name="Age (years)",limits=c(0,6),expand=c(0,0)) +
  scale_y_continuous(name="Total length (mm)",limits=c(0,200),expand=c(0,0)) +
  geom_line(data=nyLenCI,aes(x=age,y=Est)) +
  geom_line(data=nyLenCI,aes(x=age,y=LCI),linetype="dashed") +
  geom_line(data=nyLenCI,aes(x=age,y=UCI),linetype="dashed")
flp + theme_stonecat()
ggsave("doc/figures/Figure3.TIFF",width=singlewide,height=singlewide)



## Plot growth increment data across all populations ... IN THE MANUSCRIPT (Figure 4)
## Get predictions for LaPlatte R. Start with a length of 75 and 75 plus its
##   growth increment. Continue like that for several ages
LPfit1a <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,75,146),data=LPrecaps2,
               start=LPstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPboot1a <- nlsBoot(LPfit1a)
LPcf1 <- cbind(TL=c(75,146),Est=coef(LPfit1a)[c("g1","g2")],
               confint(LPboot1a)[c("g1","g2"),])

LPfit1b <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,117,165),data=LPrecaps2,
               start=LPstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPboot1b <- nlsBoot(LPfit1b)
LPcf2 <- cbind(TL=c(117,165),Est=coef(LPfit1b)[c("g1","g2")],
               confint(LPboot1b)[c("g1","g2"),])

LPfit1c <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,100,178),data=LPrecaps2,
               start=LPstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPboot1c <- nlsBoot(LPfit1c)
LPcf3 <- cbind(TL=c(100,178),Est=coef(LPfit1c)[c("g1","g2")],
               confint(LPboot1c)[c("g1","g2"),])["g2",,drop=FALSE]

LPfit1d <- nls(dtl~vbFT(mtl,mDate12,rDate12,g1,g2,w,u,100,187),data=LPrecaps2,
               start=LPstarts1,algorithm="port",
               lower=c(g1=0,g2=0,w=0,u=0),upper=c(g1=100,g2=100,w=1,u=10))
LPboot1d <- nlsBoot(LPfit1d)
LPcf4 <- cbind(TL=c(100,187),Est=coef(LPfit1d)[c("g1","g2")],
               confint(LPboot1d)[c("g1","g2"),])["g2",,drop=FALSE]

( LPPredictions <- data.frame(rbind(LPcf1,LPcf2,LPcf3,LPcf4)) %>%
    arrange(TL) %>%
    mutate(age=1:6,loc="LaPlatte R. (VT)") )

## Bootstrap CIs for increments
nyLenPred <- apply(nyboot1$coefboot,MARGIN=1,FUN=vbT,t=1:5)
nyLenMean <- apply(nyLenPred,MARGIN=1,FUN=mean)
nyLenCI <- apply(nyLenPred,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
nyIncPred <- apply(nyLenPred,MARGIN=2,FUN=diff)
nyIncMean <- apply(nyIncPred,MARGIN=1,FUN=mean)
nyIncCI <- apply(nyIncPred,MARGIN=1,FUN=quantile,probs=c(0.025,0.975))
( nyPredictions <- data.frame(age=1:5,TL=nyLenMean,tlCI=t(nyLenCI),
                              inc=c(nyIncMean,NA),
                              incCI=rbind(t(nyIncCI),c(NA,NA)),
                              loc="Great Chazy R. (NY)") )

## Results from historical studies
## Carlson (1966) Lake Vermillion, SD ... Total Length at formation, vertebrae
Carlson66 <- data.frame(age=1:7,
                        n=c(4,6,22,6,2,6,1),
                        mntl=c(78.5,96.8,113.8,137.6,155.0,175.6,193.0),
                        bctl=c(68.8,99.7,117.0,136.9,157.8,171.5,179.8),
                        loc="Vermillion R. (SD)") %>%
  mutate(inc=c(diff(mntl),NA),incbc=c(diff(bctl),NA),TL=bctl)
## Paruch (1979) Wisconsin ... Total Length at capture, pectoral spines
Paruch79 <- data.frame(age=0:5,
                        n=c(27,20,19,5,1,2),
                        mntl=c(46,102,148,159,199,187),
                        bctl=c(NA,51,95,124,152,162),
                       loc="Wisconsin") %>%
  mutate(inc=c(diff(mntl),NA),incbc=c(diff(bctl),NA),TL=bctl)
## Gilbert (1953) Streams ... Standard Length, vertebrae
Gilbert53s <- data.frame(age=1:6,
                         n=NA,
                         mnsl=c(54,73,89,104,116,129),
                         loc="Ohio streams") %>%
  mutate(TL=predict(lmSL2TL,data.frame(sl=mnsl)),
         inc=c(diff(TL),NA))
## Gilbert (1953) Lake Erie ... Standard Length, vertebrae
Gilbert53le <- data.frame(age=1:9,
                          n=NA,
                          mnsl=c(68,121,162,181,195,203,208,224,237),
                          loc="Lake Erie (OH)") %>%
  mutate(TL=predict(lmSL2TL,data.frame(sl=mnsl)),
         inc=c(diff(TL),NA))


vars <- c("age","TL","loc")
predGrowth <- rbind(LPPredictions[,vars],nyPredictions[,vars],Carlson66[,vars],
                      Paruch79[,vars],Gilbert53s[,vars],Gilbert53le[,vars])
predGrowth$loc <- ordered(predGrowth$loc,levels=c("LaPlatte R. (VT)","Great Chazy R. (NY)","Lake Erie (OH)","Wisconsin","Vermillion R. (SD)","Ohio streams"))

PG <- ggplot(data=predGrowth,aes(x=age,y=TL,group=loc)) +
  geom_line(aes(linetype=loc),size=0.5) +
  scale_linetype_manual(values=c("solid","longdash","dotted","dashed","dotdash","twodash")) +
  scale_x_continuous(name="Age (years)",
                     limits=c(1,7),breaks=seq(1,7,1)) +
  scale_y_continuous(name="Total length (mm)",
                     limits=c(0,250),expand=c(0,0))
PG + theme_stonecat() + theme(legend.position=c(.725,.15),
                              legend.key.height=unit(0.5,"lines"),
                              legend.key.width=unit(1.9,"lines"),
                              legend.title=element_blank(),
                              legend.text=element_text(size=8))
ggsave("doc/figures/Figure4.TIFF",width=singlewide,height=singlewide)
