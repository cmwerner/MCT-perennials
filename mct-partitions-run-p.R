# Set the number of processors and number of simulations to be run
n_processors <- 32 # number of nodes
n_sims <- 100

# Set the working directory and load necessary data and libraries
setwd("/project/coexistence/cwerner5/") ## CHANGE this to coexistence folder
library(parallel)
library(Rmpi)
library(dplyr)
library(tidyr)
library(tibble)


# SET PARAMETER VALUES

# seedbanking annual - nonseedbanking annual parameters
as <- 0 # changed to 0 to be a fully non-seedbanking species
ag <- 0.9
es <- 0.82
eg <- 0.6

pars.as.an <- read.csv('model_parameters_hallett.csv') %>%
  mutate(species = ifelse(species == 'Avena', 'an', 'as'), # generalize species names
         treatment = ifelse(treatment == 'consistentDry','dry','wet')) %>%
  rename(a.an = aiA, a.as = aiE) %>%
  select(-X) %>%
  mutate(seed.surv = ifelse(species == 'an', as, es),
         seed.germ = ifelse(species == 'an', ag, eg)) # add seed survival and germ columns

# perennial - nonseedbanking annual parameters
pars.p.an <- read.csv('model_parameters_perennials_theory.csv') %>%
  mutate(seed.surv = NULL)

# This determines which coexistence pair we're looking at
#pair <- 'as.an'
pair <- 'p.an'

# select the relevant parameter set
if(pair == 'as.an') {
  pars <- pars.as.an
  sp.1 <- 'an'
  sp.2 <- 'as'
  sp.order <- c('an','as') 
}
if(pair == 'p.an'){
  pars <- pars.p.an
  sp.1 <- 'an'
  sp.2 <- c('s','p')
  sp.order <- c('an','s','p')
}

# number of time steps for each run
time.warm.up <- 200
time.full <- 300

# WRITE A FUNCTION TO PASS TO THE WORKER NODES. The function should perform
#    a single simulation, and will then be repeated enough times to produce
#    the desired number of simulations.

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
      arrange(factor(species, levels = sp.order)) # need to keep them in the right order for later
    
    pars.2 <- pars.2 %>% full_join(pars.ave)
    
    
    env.draw <- rbinom(time.full, 1, env.ratio) # sequence of the environment (NOTE: would want this to change if doing multiple runs)
    env.condition <- if_else(env.draw > 0, 'wet','dry')
    
    partitions <- partition_epsilons(sp.inv = sp.2, sp.res = sp.1, pars.2, 
                                     env.condition, time.warm.up, time.full) %>%
      rbind(partition_epsilons(sp.inv = sp.1, sp.res = sp.2, pars.2, 
                               env.condition, time.warm.up, time.full)) %>%
      #pivot_longer(cols = inv.e0:res.eint) %>%
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

# EXPORT ANY NECESSARY R OBJECTS from the preprocessing stage to the indivdidual
#    processors (e.g. Do the processors need to know the value of any of the parameters
#    set up before?). The names of the objects need to be stored together as character
#    vectors because you are basically just telling R what to look for in the environment.

ObjectsToImport <- c("pars", "time.warm.up", "time.full", "pair", "sp.1", "sp.2", "sp.order") # Check that this is everything

clusterExport(cl, ObjectsToImport) # values and vectors not functions

# Run any commands necessary on the processors before running the simulation. This 
#    is where you can set the working directory of the processors, source any
#    necessary files, or load any necessary libraries
clusterEvalQ(cl, setwd("/project/coexistence/cwerner5/")) # CHANGE this to coexistence folder
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
save(simulations, file=('mct-partitions-p.RData'))
