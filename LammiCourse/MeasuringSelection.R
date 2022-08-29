#######################################################
########### - MEASURING NATURAL SELECTION - ###########
#######################################################



x=1:200
y=ceiling(x*rpois(200, 1))
plot(x,y)


summary(glm(y~x, poisson))

out=NULL
for(i in 1:100){
  out[i]=mean(rpois(20, 1))
}
hist(out)
mean(out)
sd(out)

setwd("C:/data/LammiCourse/")

#### - Ipomopsis - ####

# Read datafile
dat=read.table("campbellpowers.txt",header=T)
names(dat)

# Explore data
par(mfrow=c(1,3))
plot(dat$Corollalength, dat$Seeds)
plot(dat$Corollawidth, dat$Seeds)
plot(dat$Flowers, dat$Seeds)

# Define relative fitness
dat$relfit=dat$Seeds/mean(dat$Seeds,na.rm=T)

# Estimate variance-standardized univariate selection gradients
summary(lm(relfit~scale(Flowers), na=na.exclude, data=dat))

summary(lm(relfit~scale(Corollalength), na=na.exclude, data=dat))
summary(lm(relfit~scale(Corollawidth), na=na.exclude, data=dat))

# Estimate mean-standardized univariate selection gradients
summary(lm(relfit~scale(Flowers, scale=F), na=na.exclude, data=dat))
0.0167463*mean(dat$Flowers, na.rm=T)

summary(lm(relfit~scale(Corollalength, scale=F), na=na.exclude, data=dat))
0.06004*mean(dat$Corollalength, na.rm=T)

summary(lm(relfit~scale(Corollawidth, scale=F), na=na.exclude, data=dat))
0.04810*mean(dat$Corollawidth, na.rm=T)

# Estimate variance-standardized univariate selection gradients
dat$Year=as.factor(dat$Year)
anova(lm(relfit~scale(Flowers)*Year, na=na.exclude, data=dat))
summary(lm(relfit~scale(Flowers)*Year, na=na.exclude, data=dat))

anova(lm(relfit~scale(Corollalength)*Year, na=na.exclude, data=dat))
anova(lm(relfit~scale(Corollawidth)*Year, na=na.exclude, data=dat))

# Estimate variance-standardized multivariate selection gradients
summary(lm(relfit ~ scale(Flowers) + scale(Corollalength) + scale(Corollawidth), na = na.exclude, data = dat))
summary(lm(relfit ~ scale(Corollalength) + scale(Corollawidth), na = na.exclude, data = dat))


cor(dat[,3:5], use="pairwise")



summary(lm(relfit~Flowers + Corollalength + Corollawidth + 
             I(Flowers^2) + I(Corollalength^2) + I(Corollawidth^2) +
             Flowers:Corollalength + Flowers:Corollawidth + Corollalength:Corollawidth,
             na=na.exclude, data=dat))

summary(lm(relfit~Corollalength + Corollawidth + 
             I(Corollalength^2) + I(Corollawidth^2) +
             Corollalength:Corollawidth,
             na=na.exclude, data=dat))



#### - Dalechampia - ####

dal=read.csv("Dalechampia.csv")
head(dal)

cor(dal[,c(2,3,5,6)], use="pairwise")

dal$relfit=dal$total_seeds/mean(dal$total_seeds, na.rm=T)
dal$pollentot=dal$pollenfem + dal$pollenmale

m=lm(relfit~scale(UBA) + scale(GA) + scale(GSD) + scale(ASD), na=na.exclude, data=dal)
summary(m)

m=glm(pollenfem~scale(UBA)+scale(GA)+scale(GSD), family="poisson", na=na.exclude, data=dal)
summary(m)

library(MASS)
m=glm(pollentot~scale(UBA)+scale(GA)+scale(GSD)+scale(ASD), family="poisson", na=na.exclude, data=dal)
summary(m)

m=glm.nb(pollentot~scale(UBA)+scale(GA)+scale(GSD)+scale(ASD), na=na.exclude, data=dal)
summary(m)

library(lme4)
str(dal)
reddat=na.omit(subset(dal, select=c("pollentot", "UBA","GA", "GSD", "ASD", "patch")))

m=glmer(pollentot~scale(UBA)+scale(GA)+scale(GSD)+scale(ASD) + (1|patch), family="poisson", data=reddat)
summary(m)
