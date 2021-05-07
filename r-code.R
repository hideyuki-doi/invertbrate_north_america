library(tidyverse)  
library(gridExtra)
library(lme4)
library(lmerTest)
require(car)
##"Invert-temp6.csv" --> without 0
co<-read.csv("Invert-ori.csv",, fileEncoding="CP932")
c<-read.csv("Invert-temp5.csv")
cor<-read.csv("Invert-temp5or.csv")

##"Invert-temp6.csv" --> without 0
#ccd<-read.csv("MASS.csv")
#inv<-read.csv(file("Inv.csv",encoding='cp932'))

#cc$range2<-cc$max-cc$min --------------

c<-subset(c, c$temp>0)
hist(c$temp)

#latitude
geom_abline(  
  intercept = 42.7712,
  slope = -0.4855,
  lwd = 1.5,
  col = "blue"
)+
  
  ## LMM
  log10(Study_elevation_max+1)

c$Order2<-as.numeric(as.factor(fct_shuffle(c$Order)))
c$Study_elevation_max2<-c$Study_elevation_max/1000
c$FFG2<-as.numeric(as.factor(fct_shuffle(c$FFG)))
c$Volt2<-as.numeric(as.factor(fct_shuffle(c$Volt)))

#Max
summary(lm(Max_temp_reported~lat+(Study_elevation_min),c))
summary(lm(Min_temp_reported~lat+(Study_elevation_min),c))
summary(lm(temp~lat+(Study_elevation_min),c))



step(lm(Max_temp_reported~lat+(Study_elevation_min)+log10(Measured_length+1),c))
summary(lm(Max_temp_reported~lat+(Study_elevation_min),c))
AIC(lm(Max_temp_reported~lat+(Study_elevation_min)+log10(Measured_length+1),c))
AIC(lm(Max_temp_reported~lat,c))


#Min
y<-lm(Min_temp_reported~lat+(Study_elevation_min),c)
stepAIC(y)
summary(lm(Min_temp_reported~log10(Measured_length+1),c))


#temp
step(lm(temp~lat+(Study_elevation_min)+log10(Measured_length+1),c))
summary(lm(temp~Study_elevation_min,c))
AIC(lqm(temp~Study_elevation_min,c,tau=0.95))
AIC(lqm(temp~lat+(Study_elevation_min)+log10(Measured_length+1),c,tau=0.95))
summary(lqm(temp~log10(Measured_length+1),c,tau=0.95))

library(quantreg)

summary(lm(Max_temp_reported~lat+(Study_elevation_max)+log10(Measured_length+1),c))

summary(lm(Min_temp_reported~lat+(Study_elevation_max)+log10(Measured_length+1),c))
summary(lqm(Min_temp_reported~lat+(Study_elevation_max)+log10(Measured_length+1),c,tau=0.95))

summary(lm(temp~lat+(Study_elevation_max)+log10(Measured_length+1),c))
summary(lqm(temp~lat+(Study_elevation_max)+log10(Measured_length+1),c,tau=0.99))

summary(aov(temp~Volt*FFG*Order,c))
summary(aov(Min_temp_reported~Volt*FFG*Order,c))
summary(a2<-aov(Max_temp_reported~Volt*FFG*Order,c))
summary(a2<-aov(Max_temp_reported~FFG,c))

TukeyHSD(a2)

library(lqmm)
summary(lqm(temp~(Study_elevation_max)+log10(Measured_length+1),data=c,tau = 0.95))
summary(lqm(temp~(Study_elevation_max)+log10(Measured_length+1),data=c,tau = 0.95))

plot(Max_temp_reported~Order2,c)

y<-step(lm(temp~lat+(Study_elevation_max+1)+log10(Measured_length+1)+FFG,c))
summary(lm(temp~log10(Study_elevation_max+1),c))
summary(lmer(Max_temp_reported~lat+(1|Order)+(1|Study_location_state),c))
summary(lmer(Min_temp_reported~lat+(1|Order)+(1|Study_location_state),c))

shapiro.test(c$Study_elevation_max)



g1<-ggplot(c,aes(y=Max_temp_reported, x=lat))+
  geom_point()+
  ylim(0,32)+
  labs(x="Latitude N", y="Max. temperature (°C)")
g2<-ggplot(c,aes(y=Min_temp_reported, x=lat))+
  geom_point()+
  ylim(0,32)+
  labs(x="Latitude N", y="Min. temperature (°C)")
g3<-ggplot(c,aes(y=temp, x=lat))+
  geom_point()+
  ylim(0,32)+
  labs(x="Latitude N", y="Optimal temperature range (°C)")
grid.arrange(g1,g2,g3,nrow=1)


## LMM

summary(lmer(temp~lat+(1|Order)+(1|Study_location_state),c))
summary(lmer(Max_temp_reported~lat+(1|Order)+(1|Study_location_state),c))
summary(lmer(Min_temp_reported~lat+(1|Order)+(1|Study_location_state),c))

cor.test(co$Study_elevation_max,co$Study_elevation_min)



#Max elevation
g1<-ggplot(c,aes(y=Max_temp_reported, x=Study_elevation_max))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Max. elevation recorded (m)", y="Max. temperature (°C)")
g2<-ggplot(c,aes(y=Min_temp_reported, x=Study_elevation_max))+
  geom_point()+
  ylim(0,32)+
  labs(x="Max. elevation recorded (m)", y="Min. temperature (°C)")
g3<-ggplot(c,aes(y=temp, x=Study_elevation_max))+
  geom_point()+
  ylim(0,32)+
  labs(x="Max. elevation recorded (m)", y="Optimal temperature range (°C)")
grid.arrange(g1,g2,g3,nrow=1)


#Min elevation
g1<-ggplot(c,aes(y=Max_temp_reported, x=Study_elevation_min))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Min. elevation recorded (m)", y="Max. temperature (°C)")
g2<-ggplot(c,aes(y=Min_temp_reported, x=Study_elevation_min))+
  geom_point()+
  ylim(0,32)+
  labs(x="Min. elevation recorded (m)", y="Min. temperature (°C)")
g3<-ggplot(c,aes(y=temp, x=Study_elevation_min))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Min. elevation recorded (m)", y="Optimal temperature range (°C)")
grid.arrange(g1,g2,g3,nrow=1)



#Max elevation-latitude

g1<-ggplot(c,aes(y=Max_temp_reported, x=Study_elevation_max, col=lat))+
  geom_point()+ scale_colour_gradientn(colors = c("red","blue"))+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Max. elevation recorded (m)", y="Max. temperature (°C)")
g2<-ggplot(c,aes(y=Min_temp_reported, x=Study_elevation_max, col=lat))+
  geom_point()+scale_colour_gradientn(colors = c("red","blue"))+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Max. elevation recorded (m)", y="Min. temperature (°C)")
g3<-ggplot(c,aes(y=temp, x=Study_elevation_max, col=lat))+
  geom_point()+scale_colour_gradientn(colors = c("red","blue"))+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Max. elevation recorded (m)", y="Optimal temperature range (°C)")
grid.arrange(g1,g2,g3,nrow=1)


##Stat

#Min elevation
g1<-ggplot(c,aes(y=Max_temp_reported, x=Study_elevation_min))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Min. elevation recorded (m)", y="Max. temperature (°C)")
g2<-ggplot(c,aes(y=Min_temp_reported, x=Study_elevation_min))+
  geom_point()+
  ylim(0,32)+
  labs(x="Min. elevation recorded (m)", y="Min. temperature (°C)")
g3<-ggplot(c,aes(y=temp, x=Study_elevation_min))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Min. elevation recorded (m)", y="Optimal temperature range (°C)")
grid.arrange(g1,g2,g3,nrow=1)


## LMM

summary(lmer(temp~Study_elevation_min+(1|Order)+(1|Study_location_state),c))
summary(lmer(Max_temp_reported~Study_elevation_min+(1|Order)+(1|Study_location_state),c))
summary(lmer(Min_temp_reported~Study_elevation_min+(1|Order)+(1|Study_location_state),c))


#Body size

g1<-ggplot(c,aes(y=Max_temp_reported, x=log10(Measured_length)))+
  geom_point()+
  ylim(0,32)+
  labs(x=expression(Log[10]~Body~length~(mm)), y="Max. temperature (°C)")
g2<-ggplot(c,aes(y=Min_temp_reported, x=log10(Measured_length)))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x=expression(Log[10]~Body~length~(mm)), y="Min. temperature (°C)")
g3<-ggplot(c,aes(y=temp, x=log10(Measured_length)))+
  geom_point()+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x=expression(Log[10]~Body~length~(mm)), y="Optimal temperature range (°C)")
grid.arrange(g1,g2,g3,nrow=1)


## LMM

summary(lmer(temp~log10(Measured_length)+(1|Order)+(1|Study_location_state),c))
summary(lmer(Max_temp_reported~log10(Measured_length)+(1|Order)+(1|Study_location_state),c))
summary(lmer(Min_temp_reported~log10(Measured_length)+(1|Order)+(1|Study_location_state),c))

#FFG
c2<-subset(c,c$Volt!="NA")
g1<-ggplot(c2,aes(y=Max_temp_reported, x=Volt))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="", y="Max. temperature (°C)")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
g2<-ggplot(c2,aes(y=Min_temp_reported, x=Volt))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Voltinism", y="Min. temperature (°C)")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
g3<-ggplot(c2,aes(y=temp, x=Volt))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="", y="Optimal temperature range (°C)")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
grid.arrange(g1,g2,g3,nrow=1)


#FFG
g1<-ggplot(c,aes(y=Max_temp_reported, x=FFG))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="", y="Max. temperature (°C)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
g2<-ggplot(c,aes(y=Min_temp_reported, x=FFG))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Functional Feeding Group", y="Min. temperature (°C)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
g3<-ggplot(c,aes(y=temp, x=FFG))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="", y="Optimal temperature range (°C)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
grid.arrange(g1,g2,g3,nrow=1)

## LMM

summary(lmer(temp~FFG+(1|Order)+(1|Study_location_state),c))
summary(lmer(Max_temp_reported~FFG+(1|Order)+(1|Study_location_state),c))
summary(lmer(Min_temp_reported~FFG+(1|Order)+(1|Study_location_state),c))



#　order
g1<-ggplot(c,aes(y=Max_temp_reported, x=Order))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="", y="Max. temperature (°C)")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
g2<-ggplot(c,aes(y=Min_temp_reported, x=Order))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="Species Order", y="Min. temperature (°C)")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
g3<-ggplot(c,aes(y=temp, x=Order))+
  geom_boxplot(outlier.size=0)+
  geom_jitter(aes(col=Max_body_size),width = 0.2,show.legend = F)+
  stat_smooth(method="lm")+
  ylim(0,32)+
  labs(x="", y="Optimal temperature range (°C)")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
grid.arrange(g1,g2,g3,nrow=1)



## LMM

summary(lmer(temp~(Order)+(1|Study_location_state),c))
summary(lmer(Max_temp_reported~(Order)+(1|Study_location_state),c))
summary(lmer(Min_temp_reported~(Order)+(1|Study_location_state),c))









###App
g4<-ggplot(c,aes(y=Study_elevation_max, x=Order))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(colour=Max_body_size))+
  stat_smooth(method="lm")+
  ylim(0,25)+
  labs(x="Functional Feeding Group", y="Optimal temperature range (°C)")

ggplot(c,aes(y=Study_elevation_max, x=Order))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(colour=Max_body_size))+
  labs(x="Functional Feeding Group", y="Optimal temperature range (°C)")

ggplot(c,aes(y=Study_elevation_max, x=FFG))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(colour=Max_body_size))+
  labs(x="Functional Feeding Group", y="Optimal temperature range (°C)")

ggplot(c,aes(y=Study_elevation_max, x=Order))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(colour=Max_body_size))+
  labs(x="Functional Feeding Group", y="Optimal temperature range (°C)")

ggplot(c,aes(y=lat, x=Order))+
  geom_boxplot()+
  geom_jitter(width = 0.2, aes(colour=Max_body_size))+
  labs(x="Functional Feeding Group", y="Optimal temperature range (°C)")


##Appendix
par(mfrow=c(2,2), cex=0.6, ps=12)
plot(temp~Volt, c, ylab="Optimal temperature range (°C)", xlab="Voltinism")
plot(temp~Order, cor, ylab="Optimal temperature range (°C)", xlab="Taxon groups")
plot(temp~Max_body_size, c, ylab="Optimal temperature range (°C)", xlab="Max. body length groups")

plot(temp~Study_Citation, cc)






library(bestglm)

bestglmobject<-bestglm(wholematrix.scaled, IC="BIC", family = binomial) # Fit model
bestglmobject

###Leaps Model
library(leaps)
regsubsets.out <-
  regsubsets(Extirpation ~ RiverKm.Impounded + Km.from.mouth + snr + sigma.hf +
               sigma.lf + LENGMAT + AspectRatio + AGEMAT + LONG + FECUND + EggSize + 
               Rich.Nonnative + Nonnative + Fluvial.dep + DamIsolated + (1| huc8) + (1| Genus_species),
             data = wholematrix.scaled,
             nbest = 10,  # 10 best model for each number of predictors
             nvmax = 17,    # NULL for no limit on number of variables
             force.in = NULL, force.out = NULL,
             method = "exhaustive")
regsubsets.out
summary(regsubsets.out)
bestglmobject<-bestglm(wholematrix.scaled, IC="BIC", family = binomial)
bestglmobject

###glmmLasso Model
library(glmmLasso)
alabama<-glmmLasso(Extirpation ~ RiverKm.Impounded + Km.from.mouth + snr + sigma.hf +
                     sigma.lf + LENGMAT + AspectRatio + AGEMAT + LONG + FECUND + EggSize + 
                     Rich.Nonnative + as.factor(Nonnative) + as.factor(Fluvial.dep) + as.factor(DamIsolated), rnd = list(huc8=~1 + Genus_species), lambda=10,
                   data = wholematrix.scaled, switch.NR=FALSE, final.re=FALSE, control = list())

summary(alabama)


#BIC agrees with backward stepwise approach
out<-glm(Extirpation~., data=wholematrix.scaled, family=binomial)
step(out, k=log(nrow(wholematrix.scaled)))




