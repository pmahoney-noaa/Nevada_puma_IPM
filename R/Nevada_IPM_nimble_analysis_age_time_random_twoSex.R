##################################################
# Nevada Cougar 'Integrated' Population Model
#
# Model runtime code
#
# Code developed by: Peter Mahoney, PhD
# Model based on: Tyre et al. Unpublished
##################################################

# Load required packages
pkgs <- c('boot', 'nimble', 'tidyverse')
sapply(pkgs, require, character = T)


#########
###
### Actual data
###
#########

#### ONE ZONE, N constrained w/2 sexes
#
## Model parameters
#
nyear <- 28  # Number of years
#nzone <- 3  # Number of management zones/areas
nsim <- 5    # Number of iterations (chains)
m <- 2500    # Number of simulations
thin <- 5    # Chain thinning value

#
## Initial param values
#
inits <- list(
  # lFecA = log(2.23),
  # lHarvestFK = logit(0.02),
  # lNonHarvestFK = logit(0.02),
  # lHarvestFA = logit(0.20),
  # lNonHarvestFA = logit(0.10),
  # lSurvFa = logit(0.80),
  # lSurvMa = logit(0.80),
  # mean_fN = 1500 
)

#
## Parameters stored during the run
#
param <- c('r_a', 'r_k',  'p_k',  
           'p_fa', 'p_ma', 
           's_a', #'s_fa', 's_ma',
           'sigma_ra', 'sigma_pa', 'sigma_sa',
           'epsilon_ra', 'epsilon_pa', 'epsilon_sa',
           'fecA', 'rho', 'nKitts', 
           'propFemaleBreeding', 
           'pSexKitten',
           'Ntot', 'Nfem', 'Nmale', 'l',
           'lHarvestK', 'lNonHarvestK',
           'lNonHarvestA',
           'lHarvestFA', 'lHarvestMA', 'lSurvA',
           'mean_fN', 'mean_mN', 'fN', 'mN', 
           'oldAgeClass_f', 'oldAgeClass_m')

#
## constants of the models
#

# Combined the oldest age group due to small sample sizes
initAgeDist <- c(0.20, 0.15, 0.134, 0.160, 0.116, 0.108, 0.04, 0.032,
                 #0.015, 0.011, 0.002, 0.001, 0.001, 0.005, 0.005)
                 sum(0.025, 0.015, 0.011, 0.002, 0.001, 0.001, 0.005, 0.005)) 

consts <- list(nyears = nyear, 
               nages = length(initAgeDist), 
               #propFemaleBreeding = 0.63, 
               #pSexKitten = 0.44,
               initAgeDist = initAgeDist)

#
## data: variables to be explained in the model
#

# Age-at-harvest data
load('./NevadaIPM_Statewide_Data_by_age.rdat')
data_IPM <- lapply(data_IPM, function(x) {
  xi <- x
  for (i in 1:ncol(xi)) {
    xi[9, i] <- sum(xi[9:nrow(xi), i])
  }
  xi <- as.matrix(xi[1:9, ])
  
  return(xi[, -c(1:18)])
})

#
## Model definition and fitting
#
source('./Nevada_IPM_nimble_oneZone_c_r_a_tr_twoSex.R')  
IPMpop <- nimbleModel(code = eval(as.name('NevadaIPM')), 
                      name = 'NevadaIPM',
                      constants = consts, 
                      data = data_IPM, 
                      inits = inits, 
                      check =T)

# Model compile
CIPMpop    <- compileNimble(IPMpop)
specIPMpop <- configureMCMC(IPMpop, monitors=param, thin=thin, 
                            control = list(maxContractionsWarning = F))
IPMpopMCMC <- buildMCMC(specIPMpop, enableWAIC = T)
CIPMpopMCMC <- compileNimble(IPMpopMCMC, project = IPMpop, 
                             resetFunctions = TRUE)
set.seed(0)

# Run model
mcmc.out <- runMCMC(CIPMpopMCMC, 
                    niter = 25000, 
                    nburnin = 10000, 
                    nchains = 3,
                    samples = T, setSeed = 1001,
                    summary = T) #, WAIC = T
View(mcmc.out$summary$all.chains)

# Save output
save(mcmc.out, file = './StateWide_cr_a_tr_2sex_full.dat')


#
## Exploratory figures
#
load('./StateWide_cr_a_tr_2sex_full.dat')
out <- mcmc.out$summary$all.chains

oTot <- out[grep('Ntot', row.names(out)),] %>%
  as.data.frame() %>%
  mutate(Year = 1987:2015)

oFem <- out[grep('Nfem', row.names(out)),] %>%
  as.data.frame() %>%
  mutate(Year = 1987:2015)

oMale <- out[grep('Nmale', row.names(out)),] %>%
  as.data.frame() %>%
  mutate(Year = 1987:2015)

ggplot(oTot, aes(x = Year, y = Mean)) +
  geom_errorbar(aes(ymin = `95%CI_low`, ymax = `95%CI_upp`), alpha = 0.5) +
  geom_point() +
  geom_line()


ggplot() +
  geom_errorbar(data = oFem, aes(x = Year, y = Mean, ymin = `95%CI_low`, ymax = `95%CI_upp`), 
                color = 'red', alpha = 0.5) +
  geom_line(data = oFem, aes(x = Year, y = Mean), color = 'red') +
  geom_point(data = oFem, aes(x = Year, y = Mean), color = 'red') +
  geom_errorbar(data = oMale, aes(x = Year, y = Mean, ymin = `95%CI_low`, ymax = `95%CI_upp`),
                color = 'black', alpha = 0.5) +
  geom_line(data = oMale, aes(x = Year, y = Mean), color = 'black') +
  geom_point(data = oMale, aes(x = Year, y = Mean), color = 'black')
