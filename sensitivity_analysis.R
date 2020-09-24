# load data
infectee = read.csv('KoreaCoVID19.csv', stringsAsFactors=F)

####### Estimate incubation period

#install.packages('coarseDataTools')
library(coarseDataTools)

# choose data whose symptom onset date is known & exposure interval is given

infectee.incub = infectee[infectee$S1 != 'Asymptomatic',]
infectee.incub = infectee.incub[infectee.incub$CL !='Unknown',]

# For estimate incubation period, we only need SL, EL, ER
infectee.incub = infectee.incub[c('S1','CL','CR')]

## Change colnames for function dic.fit
# dic.fit: from coarseDataTools library, based on 
# "Estimating incubation period distributions with coarse data."
# Reich, Nicholas G., et al. Statistics in medicine 28.22 (2009): 2769-2784.

colnames(infectee.incub) = c('SL','EL','ER')
infectee.incub$SR = infectee.incub$SL
infectee.incub = infectee.incub[c(2,3,1,4)]
for(i in 1:4){
  infectee.incub[i] = as.Date(infectee.incub[,i])
}
for(i in 1:4){
  infectee.incub[i] = as.integer(infectee.incub[,i])
}

# type: 
# If symptom onset date is known and exposure duration is given by interval, 1.
# If symptom onset date is known and exposure date is known, 2 
# We don't have a case that symptom onset date is not given.

type = c()
for(i in 1:nrow(infectee.incub)){
  if(infectee.incub$EL[i] == infectee.incub$ER[i]){
    type[i]=2
  }else{
    type[i]=1
  }
}
infectee.incub$type = type

set.seed(502)
dic.fit(infectee.incub,dist="L", n.boots=1000,
        ptiles = c(0.025, 0.25, 0.5, 0.75, 0.975))

x = seq(0,20, by = 0.01)
plot(x, dlnorm(x, meanlog = 1.055, sdlog = 0.493), type = 'l', xlab = 'Days', ylab = 'Density', xlim = c(0,15), ylim = c(0,0.5), main = 'incubation period')
lines(c(-3,25), c(0,0))

####### Estimate Serial interval distribution

# choose data whose symptom onset date of both infectee and infector is known

infectee.si = infectee[infectee$S1!='Asymptomatic',]

si = as.Date(infectee.si$S1) - as.Date(infectee.si$S2)

si= as.numeric(si)
table(si)
median(si)

# minimun of serial interval: -3
# assume mu + si ~ gamma or log-normal or weibull
# set mu = 3.1 ( make mu + si >0 )
mu = 3.1
z = si + mu

## estimate Gamma distribution (MLE)
library(rGammaGamma)
set.seed(0429)
alpha = gammaMLE(z)[1]
beta = gammaMLE(z)[2]

# compute median
med.gamma = median(rgamma(5000, shape = alpha, scale = beta)) - mu

# compute bootstrap C.I.
set.seed(0429)
B = 1000
n = length(si)
med.B.gamma = c()
para1.B.gamma = c()
para2.B.gamma = c()
for(b in 1:B){
  para = gammaMLE(mu+sample(si, n, replace = T))
  para1.B.gamma[b] = para[1]
  para2.B.gamma[b] = para[2]
  med.B.gamma[b] = median(rgamma(5000, shape = para[1], scale = para[2]))
}

CI.B.gamma = round(quantile(med.B.gamma, c(0.025, 0.975)),2) - mu
round(quantile(para1.B.gamma, c(0.025, 0.975)),2)
round(quantile(para2.B.gamma, c(0.025, 0.975)),2)

#log-likelihood of estimated distribution

LL.gamma = sum(dgamma(z, shape=alpha, scale = beta, log = T))


## estimate log-Normal distribution (MLE)
med.lnorm = exp(mean(log(z))) - mu
mu.lnorm = mean(log(z))
sd.lnorm = sd(log(z))

set.seed(0429)
B = 1000
n= length(si)
med.B.lnorm = c()
para1.B.lnorm = c()
para2.B.lnorm = c()
for(b in 1:B){
  z.lnorm = mu+sample(si, n, replace = T)
  para1.B.lnorm[b] = mean(log(z.lnorm))
  para2.B.lnorm[b] = sd(log(z.lnorm))
  med.B.lnorm[b] = exp(mean(log(z.lnorm))) - mu
}

CI.B.lnorm = quantile(med.B.lnorm, c(0.025, 0.975))
round(quantile(para1.B.lnorm, c(0.025, 0.975)),2)
round(quantile(para2.B.lnorm, c(0.025, 0.975)),2)

#log-likelihood of estimated distribution
LL.lnorm = sum(dnorm(log(z), mean(log(z)), sd.lnorm, log = T))



## estimate Weibull distribution (MLE)
library(weibullness)

set.seed(0429)
para = weibull.mle(si, threshold = -mu)
med.weibull = para$scale * (log(2) ** (1/para$shape)) + para$threshold

B = 1000
n= length(si)
med.B.weibull = c()
para1.B.weibull = c()
para2.B.weibull = c()

for(b in 1:B){
  sample.B = sample(si, n, replace = T)
  para.B = weibull.mle(sample.B, threshold = -mu)
  para1.B.weibull[b] = para.B$shape
  para2.B.weibull[b] = para.B$scale
  med.B.weibull[b] = para.B$scale * (log(2) ** (1/para.B$shape)) + para.B$threshold
}

CI.B.weibull = quantile(med.B.weibull, c(0.025, 0.975))
round(quantile(para1.B.weibull, c(0.025, 0.975)),2)
round(quantile(para2.B.weibull, c(0.025, 0.975)),2)

#log-likelihood of estimated distribution
  LL.weibull = sum(dweibull(z, shape = para$shape, scale = para$scale, log = T))
