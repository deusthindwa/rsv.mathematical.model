#Gigi's model
#18/11/2022
#RSV transmission model

#====================================================================

#load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "here", "rio", "deSolve", "RColorBrewer", "reshape2")) 

#====================================================================

# T is time points. T contains both the burn-in period and evaluation period.
# N_ages is number of age groups.
#read the initial population by age group
pop1 <- readr::read_rds(here("data", "pop1.rds"))

# birth rate in each age group
# a matrix with T rows and N_ages columns
# The first column is the population birth rate
# They are born into the 0 month age group
# The rest columns are all zeros (1 months to 80+ age group)
# For U.S., you will find this information in CDC wonder
B <- readr::read_rds(here("data", "Birth_rate.rds"))

# The contact patterns in each age group.
# The patterns affect the age-specific likelihood of a susceptible individual come into contact with an infectious individual.
# The exact values does not matter because we will estimate beta_0.
# This information can be found in published literature
# Then you may consider rearrange the matrix to reflect the age groups in your model
# I wrote an r script to rearrange the contact matrix
# see contactmatrix.r under code source folder
c2 <- readr::read_rds(here("data", "c2.rds"))

