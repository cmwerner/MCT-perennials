# Set the number of processors and number of simulations to be run
n_processors <- 32 # number of nodes
n_sims <- 100

# Set the working directory and load necessary data and libraries
setwd("/project/coexistence/cwerner5/") 
library(parallel)
library(Rmpi)
library(dplyr)
library(tidyr)
library(tibble)

# CHOOSE: This determines which coexistence pair we're looking at
# case <- 'as'
 case <- 'p'

pair <- paste(case, 'an', sep = '.')

# varying <- 'lambda'
 varying <- 'alpha'

# favorable <- 'same'
# favorable <- 'opposite'
 favorable <- 'same-consistent'
# favorable <- 'opposite-consistent'
# favorable <- 'same-flipped'
# favorable <- 'opposite-flipped'

# file names

file.core <- paste0(case, '_', varying, '-', favorable)
main.wd <- '/project/coexistence/cwerner5/'

file.in <- paste0(main.wd, 'model-parameters/model-parameters_', file.core, '.csv')

file.out <- paste0(main.wd, 'model-data/mct-partitions_', 'autocorrelation_', file.core, '.RData')

# SET PARAMETER VALUES 
pars <- read.csv(file.in) 


# select the relevant parameter set
if(pair == 'as.an') {
  sp.1 <- 'an'
  sp.2 <- 'as'
  sp.order <- c('an','as')
pars <- pars %>% mutate(stem.surv = NULL) 
}
if(pair == 'p.an'){
  sp.1 <- 'an'
  sp.2 <- c('s','p')
  sp.order <- c('an','s','p')
pars <- pars %>% mutate(seed.surv = NULL)
}

# number of time steps for each run
time.warm.up <- 200
time.full <- 300

# WRITE A FUNCTION TO PASS TO THE WORKER NODES. The function should perform
#    a single simulation, and will then be repeated enough times to produce
#    the desired number of simulations.

Run_part_function <- function(i){
env.ratio <- 0.5
  
  partitions.env <- tibble(sp.invader = character(), sp.resident = character(),
                           partition = character(),
                           inv = numeric(), res = numeric(), delta = numeric(), 
                           env.auto = numeric())
  env.auto.seq <- seq(-0.9, 0.9, by = 0.1) # autocorrelation strengths
  
  for(env.auto in env.auto.seq){
    # adding rows for the weighted average of the parameters
    # this could be moved out of the for loop for the autocorrelation version but 
    # leaving it here in case we want it to be a nested loop later
    pars.2 <- pars %>% mutate(weight = ifelse(treatment == 'wet', env.ratio, 1 - env.ratio))
    
    pars.ave <- pars.2 %>% 
      group_by(species) %>%
      summarise_if(is.numeric,
                   ~ weighted.mean(., weight)) %>%
      mutate(weight = NA, treatment = 'average') %>% 
      arrange(factor(species, levels = sp.order))
    
    pars.2 <- pars.2 %>% full_join(pars.ave)
    
    
    # ADDED IN AUTOCORRELATION WHICH WILL WORK FOR ALL THE PARTITIONS EXCEPT THE STORAGE EFFECT ONE
    env.from.a <- env.ratio + env.auto/2 # symmetrical autocorrelation--likelihood of going to a from a
    env.from.b <- env.ratio - env.auto/2 # likelihood of going to a from b
    
    env.draw.a <- rbinom(time.full, 1, env.from.a) # env list if coming from a wet year (a)
    env.draw.b <- rbinom(time.full, 1, env.from.b) # env list of coming from a dry year (b)
    
    env.draw <- vector(length = time.full)
    env.draw[1] <- rbinom(1, 1, env.ratio)
    for(i in 2:time.full){
      env.draw[i] <- ifelse(env.draw[i-1]==1, env.draw.a[i], env.draw.b[i])
    }
    
    env.condition <- if_else(env.draw > 0, 'wet','dry')
    
    partitions <- partition_epsilons(sp.inv = sp.2, sp.res = sp.1, pars.2, 
                                     env.condition, time.warm.up, time.full) %>%
      rbind(partition_epsilons(sp.inv = sp.1, sp.res = sp.2, pars.2, 
                               env.condition, time.warm.up, time.full)) %>%
      pivot_longer(cols = inv.e0:res.storage) %>%
      separate(name, into = c('player', 'partition')) %>%
      pivot_wider(names_from = player, values_from = value) %>%
      mutate(delta = inv - res,
             partition = factor(partition, 
                                levels = c('full','e0','ea','el','eint', 'storage','nocov')))
    partitions$env.auto <- env.auto
    
    partitions.env <- full_join(partitions.env, partitions, 
                                by = c('sp.invader','sp.resident','partition',
                                       'inv','res', 'delta','env.auto'))
  }
  return(partitions.env)
}

# CREATE A CLUSTER AND RUN THE SIMULATIONS
cl <- makeCluster(n_processors - 1, type = "MPI")

# EXPORT ANY NECESSARY R OBJECTS from the preprocessing stage to the indivdidual
#    processors (e.g. Do the processors need to know the value of any of the parameters
#    set up before?). The names of the objects need to be stored together as character
#    vectors because you are basically just telling R what to look for in the environment.

ObjectsToImport <- c("pars", "time.warm.up", "time.full", "pair", "sp.1", "sp.2", "sp.order") # Check that this is everything

clusterExport(cl, ObjectsToImport) # values and vectors not functions

# Run any commands necessary on the processors before running the simulation. This 
#    is where you can set the working directory of the processors, source any
#    necessary files, or load any necessary libraries
clusterEvalQ(cl, setwd("/project/coexistence/cwerner5/")) 
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tidyr))
clusterEvalQ(cl, library(tibble))
clusterEvalQ(cl, source("mct-partitions-functions.R")) # put the functions used on each node

# Run the simulations on the cluster. Here x corresponds to the first argument
#    of the SimFunc function. The cluster will evaluate the SimFunc for each element
#    of the vector x below. I usually just use a vector of 1 to the number of simulations,
#    but you could get creative with this too, if you wanted. The output of the SimFunc
#    function is stored in a list which you can immediately save, or process within this
#    same script.
simulations <- clusterApply(cl, x = 1:n_sims, fun = Run_part_function) # run it x amount of times


# Save as an r dataframe
save(simulations, file=file.out)
