library(foreach)
library(pscl)
library(truncnorm)

# monte carlo EM algorithm to fit hierarchical model of logAge
MHFB = function(starAges, steps, chainDepth) {
  minLogAge = 9.5
  maxLogAge = 10.176
  
  nStars = ncol(starAges)
  
  sampledAges     = matrix(0, steps, nStars)
  sampledAges[1,] = starAges[1,] # sampledAges starts with a random draw from the ages. Row 1 is randomly chosen.

  tauSquared.scale = sum((sampledAges[1,] - mean(sampledAges[1,])) ^ 2)
  
  # Mean of inverse gamma is b / (a - 1) for a > 1
  tauSquared = tauSquared.scale / (nStars - 1)
  gamma      = rtruncnorm(1, minLogAge, maxLogAge, mean(sampledAges[1,]), sqrt(tauSquared[1] / nStars)) # optimal gamma
  
  for(step in 2:steps) {
    last.sampledAge = sampledAges[step - 1,]
    last.gamma      = gamma[step - 1]
    last.tau        = sqrt(tauSquared[step - 1])
    
    for(iter in 1:chainDepth) {
      proposedAge = NULL

      # Randomly draw one age for each star
      for(star in 1:nStars) {
        proposedAge[star] = sample(starAges[,star], 1)
      }
      
      probabilityRatio = dnorm(proposedAge,     last.gamma, last.tau) /
                       # --------------------------------------------
                         dnorm(last.sampledAge, last.gamma, last.tau)
      
      uniformDraws = runif(nStars)
      accept       = uniformDraws <= probabilityRatio
      
      last.sampledAge = ifelse(accept, proposedAge, last.sampledAge)
    }
    
    sampledAges[step,] = last.sampledAge

    tauSquared.scale = sum((sampledAges[step,] - last.gamma) ^ 2)
    
    tauSquared[step] = rigamma(1, (nStars - 1) / 2, tauSquared.scale / 2)
    gamma[step] = rtruncnorm(1, minLogAge, maxLogAge, mean(sampledAges[step,]), sqrt(tauSquared[step] / nStars))
  }
  
  list("sampledAges" = sampledAges,
       "gamma"       = gamma,
       "tauSquared"  = tauSquared)
}

main = function() {
  resultFiles = c("wd.sim.01.mcmc.res",
                  "wd.sim.02.mcmc.res",
                  "wd.sim.03.mcmc.res",
                  "wd.sim.04.mcmc.res",
                  "wd.sim.05.mcmc.res",
                  "wd.sim.06.mcmc.res",
                  "wd.sim.07.mcmc.res",
                  "wd.sim.08.mcmc.res",
                  "wd.sim.08.mcmc.res",
                  "wd.sim.10.mcmc.res",
                  "wd.sim.11.mcmc.res",
                  "wd.sim.12.mcmc.res",
                  "wd.sim.13.mcmc.res",
                  "wd.sim.14.mcmc.res",
                  "wd.sim.15.mcmc.res",
                  "wd.sim.16.mcmc.res",
                  "wd.sim.17.mcmc.res",
                  "wd.sim.18.mcmc.res")
  
  nStars = length(resultFiles)

  starAges = foreach(i = 1:nStars, .combine = "cbind") %do% {
    singlePopMcmcAges = read.table(resultFiles[i], header = T)$logAge
  }
  
  # chainDepth the chain
  fit = MHFB(starAges, steps = 10000, chainDepth = 400)
  
  fb.age   = fit$sampledAges
  fb.gamma = fit$gamma
  fb.tauSquared = fit$tauSquared
  
  save(fb.age, fb.gamma, fb.tauSquared, file = "ApproxFBsimulation1.Rdata")
}

main()