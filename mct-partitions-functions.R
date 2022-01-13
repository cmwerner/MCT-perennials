# Simple version of MCT partitioning functions to make running on Teton easier

### Species interactions one time step ------------

# Functions for one step of species interaction calculations. 
# Takes in a parameter df for all the species for one environmental condition
# and a list of species population sizes. 
# General to whether it's both species present or just one running to equilibrium
# if it's one, the other population is just at 0. 

pop_interactions <- function(N0, param){
  
  # seedbanking annual - nonseedbanking annual
  if(pair == 'as.an'){
    pars.as <- param %>% filter(species == 'as')
    pars.an <- param %>% filter(species == 'an')
    
    n.as <- N0$as * (1 - pars.as$seed.germ) * pars.as$seed.surv + 
      N0$as * pars.as$seed.germ *pars.as$lambda /
      (1 + pars.as$a.an * pars.an$seed.germ * N0$an +
         pars.as$a.as * pars.as$seed.germ * N0$as)
    
    n.an <- N0$an * pars.an$seed.germ *pars.an$lambda /
      (1 + pars.an$a.an * pars.an$seed.germ * N0$an +
         pars.an$a.as * pars.as$seed.germ * N0$as)
    
    n.all <- tibble(an = n.an,
                    as = n.as)
  }
  
  # perennial - nonseedbanking annual
  if(pair == 'p.an'){
    pars.s <- param %>% filter(species == 's')
    pars.p <- param %>% filter(species == 'p')
    pars.an <- param %>% filter(species == 'an')
    
    n.s <- N0$p*pars.p$lambda /
      (1 + pars.p$a.an * pars.an$seed.germ * N0$an +
         pars.p$a.s * pars.s$seed.germ * N0$s +
         pars.p$a.p * N0$p)
    
    n.p <- N0$p * pars.p$stem.surv +
      (N0$s * pars.s$seed.germ * pars.s$stem.surv) /
      (1 + pars.s$a.an * pars.an$seed.germ * N0$an +
         pars.s$a.s * pars.s$seed.germ * N0$s +
         pars.s$a.p * N0$p)
    
    n.an <- N0$an * pars.an$seed.germ *pars.an$lambda /
      (1+ pars.an$a.an * pars.an$seed.germ * N0$an +
         pars.an$a.s * pars.s$seed.germ * N0$s +
         pars.an$a.p * N0$p)
    
    n.all <- tibble(an = n.an,
                    s = n.s,
                    p = n.p)
  }
  return(n.all)
}


### Running an invasion sequence ------------

# Runs a single species (resident) for a number of initial time steps, and then invades the other species
# for the remaining time steps. The invading species comes in at a population of 1 at each time step
# and is reset for the next one with the growth rate saved for each invasion. 
# Resident species populations are tracked as both their population growth with the invader included
# to get the growth rate, and without the invader to determine the starting population for the next time step

# sp.invader and sp.resident parameters can be either a single value (ex: 'a') for annual species
# or multiple values (ex: c('s','p')) for perennial species. With perennial species as the invader,
# at each time step the perennial species is first run to equilibrium to determine
# what the stage ratio should be, and then this stage distribution is used for the invasion step

#For partitions where lambda is set to not vary, also setting germination and survival terms to not vary

run_invasion <- function(sp.invader, sp.resident, N0.r, parameters, env.cond,
                         vary.lambda = TRUE, vary.alpha = TRUE, covary = TRUE,
                         init.time.steps = 50, total.time.steps = 100) {
  
  
  # populations with resident species not impacted by the invader
  # NOTE would like to change this to be more general 
  
  if(pair == 'as.an'){
    n.all <- tibble(an = as.numeric(rep(NA, total.time.steps)),
                    as = as.numeric(rep(NA, total.time.steps)))
  } 
  else if(pair == 'p.an'){
    n.all <- tibble(an = as.numeric(rep(NA, total.time.steps)),
                    s = as.numeric(rep(NA, total.time.steps)),
                    p = as.numeric(rep(NA, total.time.steps)))
  }
  
  # setting initial populations
  n.all[1,sp.resident] <- N0.r
  n.all[1,sp.invader] <- 0
  
  # name of column(s) for resident species pop when impacted by invading pop
  sp.res.impact <- paste(sp.resident, 'impact', sep = '.')
  
  # parameters for the average env condition
  p.ave <- parameters %>% filter(treatment == 'average')
  
  # columns for the non-varying parameters
  alpha.columns <- c(paste('a', sp.invader, sep = '.'), 
                     paste('a', sp.resident, sep = '.'))
  lambda.columns <- parameters %>% 
    select(!c(all_of(alpha.columns),'species','treatment')) %>% 
    names()
  
  # environmental condition vectors for the no covariance option
  env.lambda <- sample(env.cond, size = length(env.cond), replace = FALSE)
  env.alpha <- sample(env.cond, size = length(env.cond), replace = FALSE)
  
  for(t in 1:init.time.steps){
    # parameters for the varying env condition
    p.env <- parameters %>% filter(treatment == env.cond[t])
    
    # constant lambda
    if(vary.lambda == FALSE){
      p.env[,lambda.columns] <- p.ave[,lambda.columns]
    }
    
    # constant alpha
    if(vary.alpha == FALSE){
      p.env[,alpha.columns] <- p.ave[,alpha.columns]
    }
    
    # parameters for no covariance between alpha and lambda
    if(covary == FALSE){
      lambda.env  <- parameters %>% filter(treatment == env.lambda[t])
      alpha.env <- parameters %>% filter(treatment == env.alpha[t])
      p.env[,lambda.columns] <- lambda.env[,lambda.columns]
      p.env[,alpha.columns] <- alpha.env[,alpha.columns]
    }
    
    # running one time step with the appropriate environmental parameters
    n.all[t+1,] <- pop_interactions(N0 = n.all[t,], param = p.env)
    
  }
  
  # figuring out rough invader population distribution if the invader is perennial
  # using average parameters and resident population size during the final warm-up step
  # makes it faster for running in for loop
  if(length(sp.invader) > 1) {
    n.inv.start <- 1/length(sp.invader)
    n.old <- n.all[init.time.steps,]
    n.old[sp.invader] <- n.inv.start
    n.old.inv <- n.inv.start[1]
    change <- 1
    
    while(change > 0.1){
      n.update <- pop_interactions(N0 = n.old, param = p.ave)
      n.per.total <- sum(n.update[sp.invader])
      n.update[sp.invader] <- n.update[sp.invader]/n.per.total
      n.update[sp.resident] <- n.old[sp.resident]
      n.old <- n.update
      
      n.new.inv <- n.update[sp.invader]
      change <- abs(n.old.inv - n.new.inv[1])
      n.old.inv <- n.new.inv[1]
    }
  } 
  
  # running the resident and invader for the remaining time steps
  for (t in (init.time.steps+1):(total.time.steps-1)){
    # parameters for the varying env condition
    p.env <- parameters %>% filter(treatment == env.cond[t])
    
    # adjusting non-varying or non-covarying parameters
    # constant lambda
    if(vary.lambda == FALSE){
      p.env[,lambda.columns] <- p.ave[,lambda.columns]
    }
    
    # constant alpha
    if(vary.alpha == FALSE){
      p.env[,alpha.columns] <- p.ave[,alpha.columns]
    }
    
    # parameters for no covariance between alpha and lambda
    if(covary == FALSE){
      lambda.env  <- parameters %>% filter(treatment == env.lambda[t])
      alpha.env <- parameters %>% filter(treatment == env.alpha[t])
      p.env[,lambda.columns] <- lambda.env[,lambda.columns]
      p.env[,alpha.columns] <- alpha.env[,alpha.columns]
    }
    
    
    # running invader and resident together
    N0.together <- n.all[t,]
    
    # figuring out starting invader population distribution if the invader is perennial
    # starting from the n.old found for average conditions above
    if(length(sp.invader) > 1){
      change <- 1
      while(change > 0.01){
        n.update <- pop_interactions(N0 = n.old, param = p.env)
        n.per.total <- sum(n.update[sp.invader])
        n.update[sp.invader] <- n.update[sp.invader]/n.per.total
        n.update[sp.resident] <- n.old[sp.resident]
        n.old <- n.update
        n.new.inv <- n.update[sp.invader]
        change <- abs(n.old.inv - n.new.inv[1])
        n.old.inv <- n.new.inv[1]
      }
      N0.together[sp.invader] <- n.update[sp.invader] # using the equilibrium stage distribution from above
    } else {
      N0.together[sp.invader] <- 1 # if the invader has only one life stage
    }
    
    
    # running one time step with the appropriate environmental parameters
    n.res.impact <- pop_interactions(N0 = N0.together, param = p.env)
    
    # saving this output as the invader population and the impacted resident population
    n.all[t+1,sp.invader] <- n.res.impact[sp.invader]
    n.all[t+1,sp.res.impact] <- n.res.impact[sp.resident]
    
    # a version of this with the resident unaffected by the invader
    N0.res.only <- n.all[t,] %>% select(all_of(sp.invader), all_of(sp.resident))
    N0.res.only[sp.invader] <- 0
    n.res.only <- pop_interactions(N0 = N0.res.only, param = p.env)
    n.all[t+1,sp.resident] <- n.res.only[sp.resident]
    
    # growth rate of the resident
    sp.res.gr <- paste('gr', paste(sp.resident, collapse = '.'), sep = '.') # column name
    n.all[t+1,sp.res.gr] <- sum(n.all[t+1,sp.res.impact])/sum(n.all[t,sp.resident])
    
    # growth rate of the invader
    sp.inv.gr <- paste('gr', paste(sp.invader, collapse = '.'), sep = '.') # column name
    n.all[t+1,sp.inv.gr] <- sum(n.all[t+1,sp.invader])
    
  }
  
  n.all$time <- 1:total.time.steps
  return(n.all)
  
}


### Calculating partitions --------------

# Function to calculate partitions for a given invader and resident,
# returns a df with invader growth rate, resident growth rate, and the invader-resident delta
# partitioned into the full growth rate, no env variation, variation in lambda, variation in alpha, 
# interaction between variation in lambda and variation in alpha. This last term is then 
# further partitioned into the storage effect (covariation) and a remnant term without covariation

partition_epsilons <- function(sp.inv, sp.res, pars, 
                               env.condition, time.warm.up, time.full){
  res.name <- paste(sp.res, collapse = '.')
  inv.name <- paste(sp.inv, collapse = '.')
  
  part <- tibble(sp.invader = inv.name, sp.resident = res.name)
  
  # no variation
  invade.none <-  run_invasion(sp.invader = sp.inv, sp.resident = sp.res,
                               N0.r = 100, parameters = pars, env.cond = env.condition,
                               init.time.steps = time.warm.up, 
                               total.time.steps = time.full,
                               vary.lambda = FALSE, vary.alpha = FALSE) %>%
    filter(time > time.warm.up + 1) %>% 
    pivot_longer(cols = starts_with('gr'),
                 names_to = 'species', names_prefix = 'gr.',
                 values_to = 'growth.rate',
                 values_drop_na = TRUE) %>% 
    mutate(gr.log = log(growth.rate))
  
  part$inv.e0 <- invade.none %>% 
    filter(species == inv.name) %>%
    select(gr.log) %>% colMeans()
  part$res.e0 <- invade.none %>% 
    filter(species == res.name) %>%
    select(gr.log) %>% colMeans()
  
  # lambda variation only
  invade.l <-  run_invasion(sp.invader = sp.inv, sp.resident = sp.res,
                            N0.r = 100, parameters = pars, env.cond = env.condition,
                            init.time.steps = time.warm.up, total.time.steps = time.full,
                            vary.lambda = TRUE, vary.alpha = FALSE) %>%
    filter(time > time.warm.up + 1) %>% 
    pivot_longer(cols = starts_with('gr'),
                 names_to = 'species', names_prefix = 'gr.',
                 values_to = 'growth.rate',
                 values_drop_na = TRUE) %>% 
    mutate(gr.log = log(growth.rate))
  
  l.inv.mean <- invade.l %>% 
    filter(species == inv.name) %>%
    select(gr.log) %>% colMeans()
  l.res.mean <- invade.l %>% 
    filter(species == res.name) %>%
    select(gr.log) %>% colMeans()
  
  part$inv.el <- l.inv.mean - part$inv.e0 
  part$res.el <- l.res.mean- part$res.e0
  
  # alpha variation only
  invade.a <-  run_invasion(sp.invader = sp.inv, sp.resident = sp.res,
                            N0.r = 100, parameters = pars, env.cond = env.condition,
                            init.time.steps = time.warm.up, total.time.steps = time.full,
                            vary.lambda = FALSE, vary.alpha = TRUE) %>%
    filter(time > time.warm.up + 1) %>% 
    pivot_longer(cols = starts_with('gr'),
                 names_to = 'species', names_prefix = 'gr.',
                 values_to = 'growth.rate',
                 values_drop_na = TRUE) %>% 
    mutate(gr.log = log(growth.rate))
  
  a.inv.mean <- invade.a %>% 
    filter(species == inv.name) %>%
    select(gr.log) %>% colMeans()
  a.res.mean <- invade.a %>% 
    filter(species == res.name) %>%
    select(gr.log) %>% colMeans()
  
  part$inv.ea <- a.inv.mean - part$inv.e0 
  part$res.ea <- a.res.mean- part$res.e0
  
  # full variation
  invade.full <-  run_invasion(sp.invader = sp.inv, sp.resident = sp.res,
                               N0.r = 100, parameters = pars, env.cond = env.condition,
                               init.time.steps = time.warm.up, total.time.steps = time.full,
                               vary.lambda = TRUE, vary.alpha = TRUE) %>%
    filter(time > time.warm.up + 1) %>% 
    pivot_longer(cols = starts_with('gr'),
                 names_to = 'species', names_prefix = 'gr.',
                 values_to = 'growth.rate',
                 values_drop_na = TRUE) %>% 
    mutate(gr.log = log(growth.rate))
  
  part$inv.full <- invade.full %>% 
    filter(species == inv.name) %>%
    select(gr.log) %>% colMeans()
  part$res.full <- invade.full %>% 
    filter(species == res.name) %>%
    select(gr.log) %>% colMeans()
  
  part$inv.eint <- part$inv.full - (part$inv.e0 + part$inv.el + part$inv.ea)
  part$res.eint <- part$res.full - (part$res.e0 + part$res.el + part$res.ea)
  
  # no correlation
  invade.nocov <- run_invasion(sp.invader = sp.inv, sp.resident = sp.res,
                               N0.r = 100, parameters = pars, env.cond = env.condition,
                               init.time.steps = time.warm.up, total.time.steps = time.full,
                               vary.lambda = TRUE, vary.alpha = TRUE,
                               covary = FALSE) %>%
    filter(time > time.warm.up + 1) %>% 
    pivot_longer(cols = starts_with('gr'),
                 names_to = 'species', names_prefix = 'gr.',
                 values_to = 'growth.rate',
                 values_drop_na = TRUE) %>% 
    mutate(gr.log = log(growth.rate))
  
  nocov.inv <- invade.nocov %>% 
    filter(species == inv.name) %>%
    select(gr.log) %>% colMeans()
  nocov.res <- invade.nocov %>% 
    filter(species == res.name) %>%
    select(gr.log) %>% colMeans()
  
  part$inv.nocov <- nocov.inv - (part$inv.e0 + part$inv.el + part$inv.ea)
  part$res.nocov <- nocov.res - (part$res.e0 + part$res.el + part$res.ea)
  
  # storage effect
  part$inv.storage <- part$inv.eint - part$inv.nocov
  part$res.storage <- part$res.eint - part$res.nocov
  
  return(part)
}
