###################################################################################################
#
#   Korea COVID-19 Dataset
#
#   72 infectees & 40 infectors
#
#   Variables
#   ID: (Unidentified) infectee ID
#   Sex: Sex of infectee
#   Age: Age of infectee
#   S1: Symptom onset date of infectee
#   Infector_ID: (Unidentified) infector ID
#   CL: Start date of the duration of exposure for infectee. If exposure is continuous, marked 'Unknown'
#   CR: End date of the duration of exposure for infectee. If exposure is continuous, same with S1.
#   S2: Symptom onset date of infector
#
###################################################################################################

infectee = read.csv('KoreaCoVID19.csv', stringsAsFactors=F)

###################################################################################################
#
#   Set prior of paramters
#
#   Model for transimission onset distribution (W = P - S)
#
#   W + mu ~  pi * gamma(phi1, phi2) + (1-pi) * U[-D,0]
#   (D = 14: fixed)
#
#   mu ~ U[0,10]
#   phi1 ~ Gamma(phi1_s, phi1_r), phi2 ~ Gamma(phi2_s, phi2_r)
#   (Set phi1_s = phi1_r = phi2_s = phi2_r = 1e-3 : very small)
#   pi ~ Beta(pi1,pi2)
#   (Set pi1 = pi2 = 1 : uniform prior)
#
###################################################################################################

set.seed(9255)

D = 14

# pi ~ Beta(pi1, pi2)
pi1 = pi2 = 1
pi = rbeta(1,pi1, pi2)

# mu ~ U[0,10]: only integer
mu = sample(0:10,1)
# if not estimate mu
# mu = 4

#phi1 ~ Gamma(phi1_s, phi1_r), phi2 ~ Gamma(phi2_s, phi2_r)
phi1_s = phi1_r = phi2_s = phi2_r = 1e-3
phi1 = dgamma(1, phi1_s, phi1_r) ; phi2 = dgamma(1, phi2_s, phi2_r)


###################################################################################################
#
#   Useful functions
#
#   - rw : generate W with parameter pi, phi1, phi2, mu (if pi=1, generate from gamma distribution)
#   - dW : density of W with parameter pi, phi1, phi2, mu
#   - rlnorm_cond : generate from conditional log-normal distribution 
#       (For impute transmission date from (current and previously estimated) incubation period distribution  
#
###################################################################################################

rW = function(pi, ph1, phi2, mu){
  coin = sample(0:1, 1,prob=c(1-pi,pi))
  if(coin){
  out = -mu + rgamma(1, phi1, phi2)
  }else{
    out = -mu + sample(-(0:D),1)
  }
  return(out)
}

dW = function(x, pi, phi1, phi2, mu){
  if(x+mu <= 0){
    d = (1-pi)/D
  }else{
    d = pi * dgamma(x+mu, phi1, phi2)
  }
  return(d)
}


rlnorm_cond = function(n, mu, sigma, lower.bound, upper.bound, replace = F){
  prob = c()
  length = upper.bound - lower.bound
  l.bound = max(lower.bound - 0.5, 0) ; u.bound = lower.bound + 0.5
  for(j in 0:length){
    prob[j+1] = plnorm(u.bound, mu, sigma) - plnorm(l.bound, mu, sigma)
    l.bound = u.bound ; u.bound = u.bound + 1
  }
  return(sample(lower.bound:upper.bound, n, replace = replace, prob = prob))
}


## infector: data frame for infector's information / current data frame is for infectee's (row: infectee's info)

infector = data.frame(ID=unique(infectee$Infector_ID))
infector$S2 = infectee$S2[match(infector$ID, infectee$Infector_ID)]

infector_L = c()
infector_U = c()
infector_P = c()
infector_n_infectee = c()
for(j in 1:nrow(infector)){
  infectee_sub = infectee[infectee$Infector_ID == infector$ID[j],]
  infector_n_infectee[j] = nrow(infectee_sub)
  if(is.element('Unknown',infectee_sub$CL)){
    infector_L[j] = -D
  }else{
    infector_L[j] = as.integer(min(as.Date(infectee_sub$CL))-as.Date(infector$S2[j]))
  }
  infector_U[j] = as.integer(max(as.Date(infectee_sub$CR))-as.Date(infector$S2[j]))
  
  if(mean(infectee_sub$CL==infectee_sub$CR)==1){
    infector_P[j] = as.character(min(as.Date(infectee_sub$CL)))
  }else{
    infector_P[j] = NA
  }
}
infector$n_infectee = infector_n_infectee
infector$L = infector_L
infector$U = infector_U
infector$P = infector_P

impute_id = which(is.na(infector$P))



### Estimated incubation period distribution: log-normal(mu.incubation, sigma.incubation)
#
# "The incubation period of coronavirus disease 2019 (COVID-19) from publicly reported confirmed cases: estimation and application." 
# Lauer, Stephen A., et al. Annals of internal medicine (2020).

mu.incubation = 1.621 ; sigma.incubation = 0.418

# Estimated from our dataset
# mu.incubation = 1.055 ; sigma.incubation = 0.493



###################################################################################################
#
#     MCMC algorithm for posterior sampling
#
#     t: MCMC step (data imputation & posterior parameter sampling). t = 1,...M(=1000)
#     print.period: print step.
#
###################################################################################################


M = 1000
print.step = 100
set.seed(425)
for(t in 1:M){
  if( t %% print.step == 0){
    print(paste(t,'th imputation & posterior sampling'))}

  ### Data imputation
  # compute acceptance rate
  acp.rate = c()
  for(j in 1:length(impute_id)){
    L = infector$L[impute_id[j]] + mu[t] ; U = infector$U[impute_id[j]] + mu[t]
    if(phi1[t] <= 1){
      if( U <= 0){
        acp.rate[j] = dW(U, pi[t], phi1[t], phi2[t], 0)
      }else if(L >0){
        acp.rate[j] = dW(L, pi[t], phi1[t], phi2[t], 0)
      }else{
        acp.rate[j] = max(dW(L, pi[t], phi1[t], phi2[t], 0), dW(1, pi[t], phi1[t], phi2[t], 0))
      }
    }else{
      mode.gamma = (phi1[t] - 1) / phi2[t]
      if( U <= 0){
        acp.rate[j] = dW(U, pi[t], phi1[t], phi2[t], 0)
      }else if( U < mode.gamma ){
        if( L <= 0 ){
          acp.rate[j] = max(dW(L, pi[t], phi1[t], phi2[t], 0), dW(U, pi[t], phi1[t], phi2[t], 0))
        }else{
          acp.rate[j] = dW(U, pi[t], phi1[t], phi2[t], 0)
        }
      }else{
        if( L < mode.gamma ){
          acp.rate[j] = dW(mode.gamma, pi[t], phi1[t], phi2[t], 0)
        }else{
          acp.rate[j] = dW(L, pi[t], phi1[t], phi2[t], 0)
        }
      } 
    }
  }


  #Data imputation
   for(j in 1:length(impute_id)){
    infectee_sub = infectee[infectee$Infector_ID == infector$ID[impute_id[j]],]
    accept = 0
    while(accept == 0){
      I = c()
      for(i in 1:nrow(infectee_sub)){
        if(infectee_sub$CL[i] == infectee_sub$CR[i]){
          # Case (i). know exact date of an exposure
          I[i] = infectee_sub$CL[i]
        }else if( infectee_sub$CL[i] != 'Unknown' ){
          # Case (ii). know duration of exposure
          I_L = as.integer(as.Date(infectee_sub$S1[i]) - as.Date(infectee_sub$CR[i]))
          I_R = as.integer(as.Date(infectee_sub$S1[i]) - as.Date(infectee_sub$CL[i]))
          I[i] = as.character(as.Date(infectee_sub$S1[i]) - rlnorm_cond(1, mu.incubation, sigma.incubation, I_L, I_R))
        }else{
          # Case (iii). continuous exposure
          I[i] = as.character(as.Date(infectee_sub$S1[i]) - round(rlnorm(1, mu.incubation, sigma.incubation)))
        }
      }
      P = min(as.Date(I))
      W = as.integer(P - as.Date(infector$S2[impute_id[j]]))
      prob.accept = dW(W, pi[t], phi1[t], phi2[t], mu[t])/acp.rate[j]
      accept = sample(0:1, 1, prob = c(1-prob.accept,prob.accept))
    }
    infector$P[impute_id[j]] = as.character(P)
  }
  
  ## Target variable(W): difference between transmisson onset date and symptom onset date of infector
  infector$W = as.integer(as.Date(infector$P) - as.Date(infector$S2))
  

  ## inference for parameters
  # initialization
  mu.mcmc = mu[t] ; phi1.mcmc = phi1[t] ; phi2.mcmc = phi2[t]

  # k : MCMC step index (k=1,...,K = 3000)
  K = 3000
  
  #sd.phi1.mh: window for metropolis- Hastings algorithm for phi1
  sd.phi1.mh = 1

  for(k in 2:K){
    #sample mu
    logp.mu = c()
    for(j in 0:10){
      logp.mu[j+1] = log(1/D) * sum(infector$W <= -j) + sum(dgamma((infector$W + j)[infector$W > -j], phi1.mcmc[k-1], phi2.mcmc[k-1], log=T) )
    }
    p.mu = exp(logp.mu - max(logp.mu))
    mu.mcmc[k] = sample(0:10, 1, prob = p.mu)
    
    # If not estimate mu,
    # mu.mcmc[k] = 4
    
    W.tmp = infector$W + mu.mcmc[k]
    ind.W.tmp = (infector$W > -mu.mcmc[k])
    
    #sample phi1
    phi1.cand = rlnorm(1, meanlog = log(phi1.mcmc[k-1]), sdlog = sd.phi1.mh)
    log.r = sum(dgamma(W.tmp[ind.W.tmp], phi1.cand,phi2.mcmc[k-1], log=T)) + dgamma(phi1.cand, phi1_s, phi1_r, log = T) - 
      sum(dgamma(W.tmp[ind.W.tmp], phi1.mcmc[k-1],phi2.mcmc[k-1], log=T)) - dgamma(phi1.mcmc[k-1], phi1_s, phi1_r, log = T)
    r = exp(log.r)
    phi1.mcmc[k] = sample(c(phi1.mcmc[k-1],phi1.cand), 1, prob = c(1-min(1,r), min(1,r)) )
    
    #sample phi2
    phi2.mcmc[k] = rgamma(1,phi2_s + phi1.mcmc[k] * sum(ind.W.tmp), phi2_r + sum(W.tmp[ind.W.tmp]))
  }
  
  #posterior sampling
  burn.in = 500
  sample.mcmc = sample(burn.in:K,1)
  mu[t+1] = mu.mcmc[sample.mcmc]
  phi1[t+1] = phi1.mcmc[sample.mcmc]
  phi2[t+1] = phi2.mcmc[sample.mcmc]
  pi[t+1] = rbeta(1,pi1 + sum(infector$W >= -mu[t+1]), pi2 + sum(infector$W < -mu[t+1]))
  
}

# check MCMC chain of mu, phi1, phi2

plot(mu, type = 'l', main = 'MCMC chain of mu')
plot(phi1, type = 'l', main = 'MCMC chain of phi1')
plot(phi2, type = 'l', main = 'MCMC chain of phi2')

# Autogorrelation plot of mu, phi1, phi2
library(coda)
autocorr.plot(mu, main = 'Autocorrelation plot of mu')
autocorr.plot(phi1, main = 'Autocorrelation plot of phi1')
autocorr.plot(phi2, main = 'Autocorrelation plot of phi2')


## prediction
# For each parameter mu, phi1, phi2, generate 5,000 samples from  -mu + gamma(phi1, phi2)
# throw away first 100 parameters (burn-in)

set.seed(420)
sample.w = c()
n.sample = 5000
for(t in 101:M){
  sample.w = rbind(sample.w, rgamma(n.sample, phi1[t], phi2[t]) - mu[t])
}

## inference : mean, sd, median, pre-symptomatic transmission probability
m=mean(sample.w)
s=sd(sample.w)
median(sample.w)
mean(sample.w<0)

mean.w = c()
for(i in 1:nrow(sample.w)){
  mean.w[i] = mean(sample.w[i,])
}
quantile(mean.w, c(0.025, 0.975))
mean(mean.w)

median.w = c()
for(i in 1:nrow(sample.w)){
  median.w[i] = median(sample.w[i,])
}
quantile(median.w, c(0.025, 0.975))
mean(median.w)

pre.sympt.w = c()
for(i in 1:nrow(sample.w)){
  pre.sympt.w[i] = mean(sample.w[i,]<0)
}
quantile(pre.sympt.w, c(0.025, 0.975))
mean(pre.sympt.w)

plot(density(sample.w[,]), xlim = c(-5,12), main = 'Predictive density of W', xlab = 'Days after symptom onset')
legend('topright', paste0('mean =',round(m,2), '\n sd=',round(s,2)), bty='n')
