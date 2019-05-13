library(sn)
library(univOutl)

# Location-scale outlier detection

# generate data from normal distr.
set.seed(123456)
r <- rsn(n=200, xi=50, omega=5, alpha=0)
hist(r, xlim=c(30,70))
mc(r) # medCouple skewnes measure

a1 <- LocScaleB(x=r, method = "IQR")
a2 <- LocScaleB(x=r, method = "MAD")
a3 <- LocScaleB(x=r, method = "Sn")
a4 <- LocScaleB(x=r, method = "Qn")
a5 <- LocScaleB(x=r, method = "scaletau2")

mpars <- rbind(IQR=a1$pars, 
               MAD=a2$pars,
               Sn=a3$pars,
               Qn=a4$pars,
               scaleTau2=a5$pars)

mbounds <- rbind(IQR=a1$bounds, 
                MAD=a2$bounds,
                Sn=a3$bounds,
                Qn=a4$bounds,
                scaleTau2=a5$bounds)
mpars
mbounds
abline(v=a1$bounds, col=2, lwd=2, lty=1)
abline(v=a2$bounds, col=3, lwd=2, lty=2)
abline(v=a3$bounds, col=4, lwd=2, lty=3)
abline(v=a4$bounds, col=5, lwd=2, lty=4)
abline(v=a5$bounds, col=6, lwd=2, lty=5)


# generate data positive skewed normal distr.
set.seed(432123)
r <- rsn(n=200, xi=50, omega=5, alpha=4)
hist(r, xlim = c(40,70))
mc(r) # medCouple skewnes measure

a1 <- LocScaleB(x=r, method = "IQR")

a2 <- LocScaleB(x=r, method = "dq")
a1$pars
a2$pars

a1$bounds
a2$bounds
abline(v=median(r), col=2, lwd=2, lty=1)
abline(v=a1$bounds, col=2, lwd=2, lty=1)
abline(v=a2$bounds, col=3, lwd=2, lty=2)


mc(r)
mc(log(r))

LocScaleB(x=r, method = "MAD", logt = TRUE)


# boxplot-based outlier detection

set.seed(11122)
r <- rsn(n=200, xi=50, omega=5, alpha=5)
hist(r, xlim = c(40, 70))

mc(r) # medCouple skewnes measure


a1 <- boxB(x=r, k=1.5, method='resistant')
a2 <- boxB(x=r, k=1.5, method='asymmetric')
a3 <- boxB(x=r, k=1.5, method='adjbox')

mfenc <- rbind(std=a1$fences, 
               asym=a2$fences, 
               adjb=a3$fences)
mfenc
abline(v = median(r), col=2, lwd=2, lty=1)
abline(v = a1$fences, col=2, lwd=2, lty=3)
abline(v = a2$fences, col=3, lwd=2, lty=2)
abline(v = a3$fences, col=4, lwd=2, lty=1)

#####
# Hidiroglou-Berthelot ratios
# rice <- read_delim("AAA-Work/progetti-GdL/univOutl/2017-11 - uRos 2017 Bucharest/rice.csv", ";", escape_double = FALSE, trim_ws = TRUE)
# 
# rice1 <- subset(rice, Year==2014)
# rice2 <- subset(rice, Year==2015)
# 
# lab <- c("geographicAreaM49", "GeographicArea", "Value")
# mm <- merge(rice1[lab], rice2[lab], 
#             by=c("geographicAreaM49", "GeographicArea"),
#             all=TRUE, suffixes = c('2014', '2015'))
# RiceProd <- mm
# save(RiceProd, file='RiceProd.rda')

load(file='rice.rda')

outlRice <- HBmethod(yt1 = rice$Prod2014, 
                     yt2 = rice$Prod2015,
                     return.dataframe = TRUE, 
                     C=15)

outlRice$quartiles.E
outlRice$bounds.E

hist(outlRice$data$Escore, xlim=c(-2000, 1500))
abline(v = outlRice$quartiles.E['50%'], col=2, lwd=2, lty=1)
abline(v = outlRice$bounds.E, col=2, lwd=2, lty=3)

head(outlRice$excluded, 3)
head(outlRice$data, 3)

outl.HB <- outlRice$data$id[outlRice$data$outliers==1]


# ratioSize
oo <- ratioSize(numerator = rice$Prod2015,
                denominator = rice$Prod2014, 
                return.dataframe = T)
oo$median.r          
oo$bounds 
hist(oo$data$c.ratio)
abline(v = median(oo$data$c.ratio), col=2, lwd=2, lty=1)
abline(v = oo$bounds, col=2, lwd=2, lty=3)

head(oo$data, 3)

oo <- ratioSize(numerator = rice$Prod2015,
                denominator = rice$Prod2014, 
                return.dataframe = T, size.th = 1000)
head(oo$data)
 
# adjusted boxplot on E-scores
outlRice <- HBmethod(yt1 = rice$Prod2014, 
                     yt2 = rice$Prod2015,
                     return.dataframe = TRUE, 
                     C=5.4)

oo <- boxB(x=outlRice$data$Escore, 
           method = 'adjbox')
outlRice$bounds.E
oo$fences
outlRice.adj <- oo$outliers

intersect(outl.HB, outlRice.adj)
outl.HB
outlRice.adj

hist(outlRice$data$Escore, xlim=c(-2000, 1500))
abline(v = outlRice$quartiles.E['50%'], col=2, lwd=2, lty=1)
abline(v = outlRice$bounds.E, col=2, lwd=2, lty=3)

abline(v = oo$fences, col=3, lwd=2, lty=4)

####
## bivariate outlier detection
# 

load(file="rice.rda")

par(mfrow=c(2,2))
par(pty='s')
hist(rice$Area2014, main='Rice-Prod_2014', xlab='')
hist(rice$Prod2015, main='Rice-Prod_2015', xlab='')
hist(log(rice$Prod2014+1), main='Log(Rice-Prod_2014+1)', xlab='')
hist(log(rice$Prod2015+1), main='Log(Rice-Prod_2015+1)', xlab='')

library("splines")
library("ggplot2")
library("MASS")


# scatterplot con linear regr
rice$logProd2014 <- log(rice$Prod2014+1)
rice$logProd2015 <- log(rice$Prod2015+1)

ggplot(rice, aes(x=logProd2014, y=logProd2015)) +
    geom_point(shape=1) +    # Use hollow circles
    geom_smooth(method=lm) +
    labs(x='Log(Rice-Prod_2014+1)', y='Log(Rice-Prod_2015+1)') 

# # scatterplot with robust regression
# ggplot(rice, aes(x=z14, y=z15)) +
#     geom_point(shape=1)+
#     stat_smooth(method="rlm",fullrange=TRUE)
# 
# # scatterpot with splines
# ggplot(rice, aes(x=z14, y=z15)) +
#     geom_point() +
#     stat_smooth(method = "lm", formula = y ~ ns(x, 3), level=0.9999999) +
#     

library("mvoutlier")
par(mfrow=c(1,2))
corr.plot(rice$logProd2014, rice$logProd2015, alpha=0.01)
dd <- dd.plot(rice[,c("logProd2014", "logProd2015")], alpha=0.01)
#sp <- symbol.plot(cbind(rice$z14, rice$z15), alpha=0.01)
# cp <- color.plot(cbind(rice$z14, rice$z15), alpha=0.01)

# vcz.mcd <- covMcd(x=rice[, c('z14','z15')])
# plot(vcz.mcd, which="tolEllipsePlot", classic=TRUE, cutoff=0.01)
# plot(vcz.mcd, which="dd", cutoff=0.01)

sum(dd$outliers)
outl <- dd$outliers
head(rice[outl,], 3)

par(mfrow=c(1,1))
plot(x=rice$logProd2014[!outl], y=rice$logProd2015[!outl], 
     col=3, lwd=1, xlim=c(0,18), ylim=c(0,18),
     xlab='Log(Rice-Prod_2014+1)', ylab='Log(Rice-Prod_2015+1)')
abline(0,1)
points(x=rice$logProd2014[outl], y=rice$logProd2015[outl], 
       pch='+', col=2, lwd=3)
text(x=rice$logProd2014[outl], y=rice$logProd2015[outl], 
     labels=rice$Geographic.Area[outl], cex=0.8, pos=4)


# SeleMix, joint model (no predictor variable without errors)
library("SeleMix")
out.sel <- ml.est(y = rice[,c("logProd2014", "logProd2015")], 
                  model="N", w=0.005, w.fix=F, t.outl=0.5)
out.sel$w # estimated proportion of contaminated data
out.sel$lambda # estimated variance inflation factor
sum(out.sel$outlier) # estimated number of contamined obs
sum(out.sel$tau==1)
sum(out.sel$tau>0.98)

toCheck <- data.frame(Geographic.Area=rice$Geographic.Area,
                      postProb=out.sel$tau,
                      rice[,c("logProd2014", "logProd2015")],
                      out.sel$ypred)

toCheck <- subset(toCheck, postProb>0.5)
toCheck <- toCheck[order(toCheck$postProb, decreasing = T), ]
head(toCheck)
outl.mix <- as.logical(out.sel$outlier)

### 
par(mfrow=c(1,1))
plot(x=rice$logProd2014[!outl], y=rice$logProd2015[!outl], 
     col=3, lwd=1, xlim=c(0,18), ylim=c(0,18),
     xlab='Log(Rice-Prod_2014+1)', ylab='Log(Rice-Prod_2015+1)')
abline(0,1)
points(x=rice$logProd2014[outl], y=rice$logProd2015[outl], 
       pch='x', col=2, lwd=3)
points(x=rice$logProd2014[outl.mix], y=rice$logProd2015[outl.mix], 
       pch='+', col=4, lwd=3)
text(x=rice$logProd2014[outl], y=rice$logProd2015[outl], 
     labels=rice$Geographic.Area[outl], cex=0.8, pos=4)


