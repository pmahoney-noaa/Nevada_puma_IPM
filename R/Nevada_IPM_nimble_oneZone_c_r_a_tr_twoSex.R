##################################################
# Nevada Cougar 'Integrated' Population Model
#
# Nimble Model Definition File
#
# Code developed by: Peter Mahoney, PhD
# Model based on: Tyre et al. Unpublished
##################################################

# Base Model
NevadaIPM <- nimbleCode({
  ##-------------------------------------------------
  ## 0. Year 1 initialization
  ##-------------------------------------------------
  
  # Initial populations size
  mean_fN ~ dunif(1800, 2000)  # Female starting population size, choose carefully
  mean_mN ~ dunif(1800, 2000)  # Male starting population size, choose carefully
  
  # Assign ages
  for (a in 1:nages) {
    fN[a, 1] ~ T(dpois(mean_fN * initAgeDist[a]), max(fX[a, 1], 1), )
    mN[a, 1] ~ T(dpois(mean_mN * initAgeDist[a]), max(mX[a, 1], 1), ) 
  }
  
  ##-------------------------------------------------
  ## 1. Define the priors for the parameters
  ##-------------------------------------------------
  
  ## Base non-harvest rate
  lNonHarvestK ~ dnorm(logit(0.01), sd = 0.50)  # For kittens
  lNonHarvestA ~ dnorm(logit(0.03), sd = 0.25)  # For 1+
  
  sigma_ra ~ T(dnorm(0, 0.2), 0, )
  for (t in 1:nyears) {
    epsilon_ra[t] ~ dnorm(0, sd = sigma_ra) #L: 0.1
  }
  
  ## Base harvest rate
  lHarvestK ~ dnorm(logit(0.01), sd = 0.5)                       # For kittens 
  lHarvestFA ~ T(dnorm(logit(0.08), sd = 0.25), logit(0.015), )  # For adult females (1+)
  lHarvestMA ~ T(dnorm(logit(0.08), sd = 0.25), logit(0.015), )  # For adult males (1+)
  
  sigma_pa ~ T(dnorm(0, 0.2), 0, )
  for (t in 1:nyears) {
    epsilon_pa[t] ~ dnorm(0, sd = sigma_pa) #L: 0.1
  }
  
  ## Unobserved Survival (logit scale), 1+ yo
  lSurvA ~ T(dnorm(logit(0.80), sd = 0.5), logit(0.60), logit(0.90))
  
  sigma_sa ~ T(dnorm(0, 0.2), 0, )
  for (t in 1:nyears) {
    epsilon_sa[t] ~ dnorm(0, sd = sigma_sa) #L: 0.1
  }
  
  
  ## Fecundity (log scale), adults only
  lFecA ~ dunif(0, 0.99)                    # Adult fecundity
  lFB ~ dunif(logit(0.60), logit(0.70))     # expit probability of female breeding
  
  ## Sex at birth (static/deterministic)
  #pSexKitten = 0.5  # Ashman et al. ~ 1 Female : 1.255 Male (pSexKitten = 0.44)
  pSK  ~ dunif(logit(0.40), logit(0.5))
  
  #-------------------------------------------------
  # 2. Derived parameters
  #-------------------------------------------------
  
  # Per-capita non-harvest mortality (r[a])
  logit(r_k) <- lNonHarvestK
  for (t in 1:nyears) {
    logit(r_a[t]) <- lNonHarvestA + epsilon_ra[t]
  }

  # Harvest rate (p[a,t])
  logit(p_k) <- lHarvestK
  for (t in 1:nyears) {
    logit(p_fa[t]) <- lHarvestFA + epsilon_pa[t]
    logit(p_ma[t]) <- lHarvestMA + epsilon_pa[t]
  }
  
  # Unobserved survival rate (s[t])
  #logit(s_a) <- lSurvA
  for (t in 1:nyears) {
    logit(s_a[t]) <- lSurvA + epsilon_sa[t]
  }

  
  # Fecundity
  log(fecA) <- lFecA
  logit(propFemaleBreeding) <- lFB
  
  # Sex of kittens
  logit(pSexKitten) <- pSK
  
  # Population size
  for (t in 1:(nyears + 1)) {
    Nfem[t] <- sum(fN[1:nages, t])
    Nmale[t] <- sum(mN[1:nages, t])
    Ntot[t] <- Nfem[t] + Nmale[t]
  }
  
  # Population growth rate
  for (t in 1:nyears) {
    l[t] <- Ntot[t+1] / (Ntot[t] + 0.00001)
  }
  
  
  #-------------------------------------------------
  # 3. The likelihoods of the single data sets
  #-------------------------------------------------
  
  for (t in 1:nyears) {
    # Likelihood for population population count data (state-space model)
    
    ####
    ## Observation Process
    ####
    
    # Non-harvest
    fR[1, t] ~ dbin(prob = r_k, size = fN[1, t])
    mR[1, t] ~ dbin(prob = r_k, size = mN[1, t])
    
    # Harvest
    fC[1, t] ~ dbin(prob = p_k, size = fN[1, t])
    mC[1, t] ~ dbin(prob = p_k, size = mN[1, t])
    
    for (a in 2:nages) {
      # Non-harvest
      fR[a, t] ~ dbin(prob = r_a[t], size = fN[a, t])
      mR[a, t] ~ dbin(prob = r_a[t], size = mN[a, t])
      
      # Harvest
      fC[a, t] ~ dbin(prob = p_fa[t], size = fN[a, t])
      mC[a, t] ~ dbin(prob = p_ma[t], size = mN[a, t])
    }
    
    
    ####
    ## System Process
    ####
    
    # Unobserved survival
    for (a in 1:(nages - 2)) {
      fN[a+1, t+1] ~ T(dbin(prob = s_a[t], size = fN[a, t] - fC[a, t] - fR[a, t]), fX[a, t], )
      mN[a+1, t+1] ~ T(dbin(prob = s_a[t], size = mN[a, t] - mC[a, t] - mR[a, t]), mX[a, t], )
    }
    
    oldAgeClass_f[t] <- sum(fN[(nages-1):nages, t])
    fN[nages, t+1] ~ T(dbin(prob = s_a[t], size = oldAgeClass_f[t] - fC[nages, t] - fR[nages, t]), fX[nages, t], )
    
    oldAgeClass_m[t] <- sum(mN[(nages-1):nages, t])
    mN[nages, t+1] ~ T(dbin(prob = s_a[t], size = oldAgeClass_m[t] - mC[nages, t] - mR[nages, t]), mX[nages, t], )
    
    # Reproduction
    nf_adults[t] <- sum(fN[3:nages, t]) - sum(fC[3:nages, t]) - sum(fR[3:nages, t])
    rho[t] <- (nf_adults[t] * propFemaleBreeding) * fecA # pSexKitten = 0.5 for parity
    
    # Next round alive
    nKitts[t] ~ T(dpois(rho[t]), fX[1, t] + mX[1, t], )
    fN[1, t+1] ~ T(dbin(prob = pSexKitten, size = nKitts[t]), fX[1, t], )
    mN[1, t+1] <- nKitts[t] - fN[1, t+1]
  }
})