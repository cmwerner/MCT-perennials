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


# Function to determine the correct stage-structured distribution for perennial as invader
perennial_distribution <- function(change.limit = 0.01, n.old = N0.together, 
                                   p.env = p.env, sp.resident, sp.invader) {
  n.old.inv <- n.old[sp.invader]
  change <- 1
  
  while(change > change.limit){
    n.update <- pop_interactions(N0 = n.old, param = p.env) # run one time step
    
    n.per.total <- sum(n.update[sp.invader]) # total perennial population
    n.update[sp.invader] <- n.update[sp.invader]/n.per.total # normalize and update perennial population
    n.update[sp.resident] <- n.old[sp.resident] # keep resident population the same (do not update)
    
    # calculate change in perennial population distribution
    n.new.inv <- n.update[sp.invader]
    change <- abs(n.old.inv[1] - n.new.inv[1])
    
    # reset values for next loop
    n.old <- n.update
    n.old.inv <- n.new.inv[1]
  }
  return(n.update) # return dataframe with all population size values
}


# Function to run one time series, first with sp.resident alone for init.time.steps alone 
# starting from population size N0.r. 
# Then the invader species is added in at each time step at low density 
# from init.time.steps to total.time.steps. 
# Parameter names sp.invader and sp.resident can be either a single value 
# (ex: 'a') or multiple values (ex: c('s','p')) to accommodate stage-structured species.

# This function runs with all parameters varying. 
# All the column names are based on the actual species names.

# First run with all parameters varying. This is where stage structure distribution is calculated
# and this what is used for all the following non-varying runs
run_full <- function(sp.invader, sp.resident, N0.r, parameters, env.cond,
                     init.time.steps = 50, total.time.steps = 100) {
  
  # populations with resident species not impacted by the invader
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
  # names of column(s) for invader stage distribution
  sp.inv.dist <- paste(sp.invader, 'dist', sep = '.')
  
  for(t in 1:init.time.steps){
    # parameters all varying
    p.env <- parameters %>% filter(treatment == env.cond[t])
    
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
    
    # while-loop with updates happens in perennial_distribution function
    n.old <- perennial_distribution(change.limit = 0.1, n.old = n.old, p.env = p.env,
                                    sp.resident, sp.invader) 
  }
  
  # running the resident and invader for the remaining time steps
  for (t in (init.time.steps+1):(total.time.steps-1)){
    # parameters for the varying env condition
    p.env <- parameters %>% filter(treatment == env.cond[t])
    
    # running invader and resident together
    N0.together <- n.all[t,]
    
    # figuring out starting invader population distribution if the invader is perennial
    # starting from the n.old found for average conditions above
    if(length(sp.invader) > 1){
      n.old[sp.resident] <- n.all[t,sp.resident]
      n.update <- perennial_distribution(change.limit = 0.01, n.old = n.old, p.env = p.env,
                                         sp.resident, sp.invader) # function
      
      N0.together[sp.invader] <- n.update[sp.invader] # using the equilibrium stage distribution from above
    } else {
      N0.together[sp.invader] <- 1 # if the invader has only one life stage
    }
    n.all[t,sp.inv.dist] <- N0.together[sp.invader]
    
    
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

# This function is similar to the above but modified to look at the impact of 
# different varying parameters
# this uses the output from the full run above
run_varying <- function(sp.invader, sp.resident, parameters, env.cond,
                        vary.lambda = TRUE, vary.alpha = TRUE, covary = TRUE,
                        init.time.steps = 50, total.time.steps = 100, n.full) {
  
  
  # populations with resident species not impacted by the invader
  if(pair == 'as.an'){
    n.all <- tibble(an = as.numeric(rep(NA, total.time.steps)),
                    as = as.numeric(rep(NA, total.time.steps)))
  } 
  else if(pair == 'p.an'){
    n.all <- tibble(an = as.numeric(rep(NA, total.time.steps)),
                    s = as.numeric(rep(NA, total.time.steps)),
                    p = as.numeric(rep(NA, total.time.steps)))
  }
  
  
  # name of column(s) for resident species pop when impacted by invading pop
  sp.res.impact <- paste(sp.resident, 'impact', sep = '.')
  # names of column(s) for invader stage distribution
  sp.inv.dist <- paste(sp.invader, 'dist', sep = '.')
  
  # parameters for the average env condition
  p.ave <- parameters %>% filter(treatment == 'average')
  
  # columns used for the non-varying parameters
  alpha.columns <- c(paste('a', sp.invader, sep = '.'), 
                     paste('a', sp.resident, sep = '.'))
  lambda.columns <- parameters %>% 
    select(!c(all_of(alpha.columns),'species','treatment')) %>% 
    names()
  
  # environmental condition vectors for the no covariance option
  # randomizes lambda columns to break the lambda - alpha covariance
  env.lambda <- sample(env.cond, size = length(env.cond), replace = FALSE)
  
  
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
    # randomizes lambda columns to break the lambda - alpha covariance
    if(covary == FALSE){
      lambda.env  <- parameters %>% filter(treatment == env.lambda[t])
      p.env[,lambda.columns] <- lambda.env[,lambda.columns]
    }
    
    
    # Starting population sizes are based on the fully varying simulation only
    # Resident starting population size is un-impacted (without the invader present)
    # Invader total starting population size is 1, stage distribution is the same
    # stage distribution from the fully varying run
    
    N0.together <- n.full[t,] 
    if(length(sp.invader) > 1){
      N0.together[sp.invader] <- N0.together[sp.inv.dist] 
    } else {
      N0.together[sp.invader] <- 1 # if the invader has only one life stage
    }
    
    
    # running one time step with the appropriate environmental parameters
    n.res.impact <- pop_interactions(N0 = N0.together, param = p.env)
    
    # saving this output as the invader population and the impacted resident population
    n.all[t+1,sp.invader] <- n.res.impact[sp.invader]
    n.all[t+1,sp.res.impact] <- n.res.impact[sp.resident]
    
    # growth rate of the resident from full variation starting point
    sp.res.gr <- paste('gr', paste(sp.resident, collapse = '.'), sep = '.') # column name
    n.all[t+1,sp.res.gr] <- sum(n.all[t+1,sp.res.impact])/sum(N0.together[sp.resident])
    
    # growth rate of the invader (always from an initial total population of 1)
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
  
  # full variation
  invade.full <-  run_full(sp.invader = sp.inv, sp.resident = sp.res,
                           N0.r = 100, parameters = pars, env.cond = env.condition,
                           init.time.steps = time.warm.up, total.time.steps = time.full) 
  invade.full.2 <- invade.full %>%
    filter(time > time.warm.up + 1) %>% 
    pivot_longer(cols = starts_with('gr'),
                 names_to = 'species', names_prefix = 'gr.',
                 values_to = 'growth.rate',
                 values_drop_na = TRUE) %>% 
    mutate(gr.log = log(growth.rate))
  
  part$inv.full <- invade.full.2 %>% 
    filter(species == inv.name) %>%
    select(gr.log) %>% colMeans()
  part$res.full <- invade.full.2 %>% 
    filter(species == res.name) %>%
    select(gr.log) %>% colMeans()
  
  # no variation full run
  invade.none <-  run_varying(sp.invader = sp.inv, sp.resident = sp.res,
                              parameters = pars, env.cond = env.condition,
                              init.time.steps = time.warm.up, 
                              total.time.steps = time.full,
                              vary.lambda = FALSE, vary.alpha = FALSE, covary = TRUE,
                              n.full = invade.full) %>%
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
  invade.l <- run_varying(sp.invader = sp.inv, sp.resident = sp.res,
                          parameters = pars, env.cond = env.condition,
                          init.time.steps = time.warm.up, total.time.steps = time.full,
                          vary.lambda = TRUE, vary.alpha = FALSE, covary = TRUE,
                          n.full = invade.full) %>%
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
  invade.a <-  run_varying(sp.invader = sp.inv, sp.resident = sp.res,
                           parameters = pars, env.cond = env.condition,
                           init.time.steps = time.warm.up, total.time.steps = time.full,
                           vary.lambda = FALSE, vary.alpha = TRUE, covary = TRUE,
                           n.full = invade.full) %>%
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
  
  part$inv.eint <- part$inv.full - (part$inv.e0 + part$inv.el + part$inv.ea)
  part$res.eint <- part$res.full - (part$res.e0 + part$res.el + part$res.ea)
  
  # no correlation
  invade.nocov <- run_varying(sp.invader = sp.inv, sp.resident = sp.res,
                              parameters = pars, env.cond = env.condition,
                              init.time.steps = time.warm.up, total.time.steps = time.full,
                              vary.lambda = TRUE, vary.alpha = TRUE, covary = FALSE,
                              n.full = invade.full) %>%
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

