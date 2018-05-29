###########################################
############# MANUSCRIPT CODE #############
###########################################


###########################################
##------------ Load Packages ------------##
###########################################

library(INLA) # Version 17.12.01

library(sp) # Version 1.2-5
library(rgdal) # Version 1.2-13
library(maptools) # Version 0.9-2
library(ggplot2) # Version 2.2.1
library(rgeos) # Version 0.3-25
library(plyr) # Version 1.8.4
library(dplyr) # Version 0.7.4
library(raster) # Version 2.5-8

library(gridExtra) # Version 2.3
library(cowplot) # Version 0.9.2

library(TMB) # Version 1.7.11
library(coda) # Version 0.19-1

###########################################
##------- Load Helper Functions ---------##
###########################################

# Kenya
source("Code/helper_fcns/make_prediction_A_matrix.R")
source("Code/helper_fcns/obtain_results.R")

# Scotland
source("Code/helper_fcns/createU_EB.R")
source("Code/helper_fcns/createU_Hybrid.R")
source("Code/helper_fcns/setup_HMC_hybrid.R")
source("Code/helper_fcns/setup_HMC_fully_bayesian.R")


###########################################
##------------- Simulation --------------##
###########################################

source("Code/setup_simulation.R")
source("Code/run_simulation.R")


###########################################
##------------ Data Analysis ------------##
###########################################

source("Code/setup_scotland.R")
run_fb <- FALSE # switch to TRUE to run fully Bayesian approach
source("Code/run_scotland.R")
