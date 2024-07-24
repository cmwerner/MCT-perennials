# MCT-perennials

## Code used in manuscript, running on computer cluster
mct-partitions-functions.R includes all the functions used for the simulations and analyses
mct-partitions-run.R is the wrapper code for executing the code on the computer cluster

## Code for generating figures
MCT-perennials_partitioning-figures.Rmd is the figure plotting code using the data generated by the simulation runs.
Inputs: data generated by simulation runs, provided in the model-data folder. Simulation outputs are saved as .RData files
Dependencies: 'tidyverse' and 'here' packages
Outputs: Figure 2, Figure 3, Supplemental Figure S2, Supplemental Figure S3, Supplemental Figure S4, Supplemental Figure S5, Supplemental Figure S6, Supplemental Figure S7, 

## Code for peer analysis
MCT-perennials_partitioning-code.Rmd is a simplified version of the code used in this manuscript which can be easily run on an individual computer, using only one simulation per parameter and environmental ratio combination rather than many runs.
Inputs: model parameters, provided in the files in the model-parameters folder.
Dependencies: 'tidyverse' and 'here' packages
Outputs: Supplemental Figure S1 is generated by the last section of this code. 