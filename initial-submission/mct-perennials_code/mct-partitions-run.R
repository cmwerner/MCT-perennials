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

# This determines which coexistence pair we're looking at
# case <- 'as'
 case <- 'p'

pair <- paste(case, 'an', sep = '.')

varying <- 'alpha'

# This determines which environmental scenario we're running
 favorable <- 'opposite-consistent'
# favorable <- 'same-consistent'
# favorable <- 'same-flipped'

 
# FILE NAMES
file.core <- paste0(case, '_', varying, '-', favorable)
main.wd <- '/project/coexistence/cwerner5/'

file.in <- paste0(main.wd, 'model-parameters/model-parameters_', file.core, '.csv')

file.out <- paste0(main.wd, 'model-data/mct-partitions_', file.core, '.RData')

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

# FUNCTION TO PASS TO THE WORKER NODES. 
# This function performs a single simulation and will be repeated

Run_part_function <- function(i){
  
  partitions.env <- tibble(sp.invader = character(), sp.resident = character(),
                           partition = character(),
                           inv = numeric(), res = numeric(), delta = numeric(), 
                           env.ratio = numeric())
  env.seq <- seq(0, 1, by = 0.05) # ratio of wet years to dry years
  
  for(env.ratio in env.seq){
    # adding rows for the weighted average of the parameters
    pars.2 <- pars %>% mutate(weight = ifelse(treatment == 'wet', env.ratio, 1 - env.ratio))
    
    pars.ave <- pars.2 %>% 
      group_by(species) %>%
      summarise_if(is.numeric,
                   ~ weighted.mean(., weight)) %>%
      mutate(weight = NA, treatment = 'average') %>% 
      arrange(factor(species, levels = sp.order))
    
    pars.2 <- pars.2 %>% full_join(pars.ave)
    
    
    env.draw <- rbinom(time.full, 1, env.ratio) # sequence of the environment
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
    partitions$env.ratio <- env.ratio
    
    partitions.env <- full_join(partitions.env, partitions, 
                                by = c('sp.invader','sp.resident','partition',
                                       'inv','res', 'delta','env.ratio'))
  }
  return(partitions.env)
}

# CREATE A CLUSTER AND RUN THE SIMULATIONS
cl <- makeCluster(n_processors - 1, type = "MPI")

# EXPORT NECESSARY R OBJECTS from the preprocessing stage to the indivdidual processors 

ObjectsToImport <- c("pars", "time.warm.up", "time.full", "pair", "sp.1", "sp.2", "sp.order") 
clusterExport(cl, ObjectsToImport) # values and vectors not functions

# Run commands necessary on the processors before running the simulation, 
# including setting the working directory of the processors, sourcing necessary files, 
# and loading necessary libraries
clusterEvalQ(cl, setwd("/project/coexistence/cwerner5/")) 
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, library(tidyr))
clusterEvalQ(cl, library(tibble))
clusterEvalQ(cl, source("mct-partitions-functions.R")) # put the functions used on each node

# Run the simulations on the cluster n_sims times. 
# The cluster will evaluate the SimFunc for each element of the vector x below. 
simulations <- clusterApply(cl, x = 1:n_sims, fun = Run_part_function) 

# Save as an r dataframe
save(simulations, file=file.out)
