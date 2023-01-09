#Gigi's model
#18/11/2022
#RSV transmission model

#====================================================================

#load packages for analysis
pacman::p_load(char = c("lubridate", "tidyverse", "dplyr", "here", "rio", "deSolve", "RColorBrewer", "reshape2")) 

#====================================================================

# T is time points. T contains both the burn-in period and evaluation period.
# N_ages is number of age groups.
# Read the initial population by age group
Pop1 <- readr::read_rds(here("data", "pop1.rds"))

# Birth rate in each age group
# A matrix with T rows and N_ages columns
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
# See contactmatrix.r under code source folder
c2 <- readr::read_rds(here("data", "c2.rds"))

# Could replace this with vector of actual age names
N_ages <- length(Pop1)
agenames <- paste0('Agegrp', 1:N_ages)

# Define states of transmission compartmental model
StateNames <- c('M', 'S0', 'I1', 'S1', 'I2', 'S2', 'I3', 'S3', 'I4')

# N age groups x K states
yinit.matrix <- array(NA, dim = c(N_ages, length(StateNames)))

# Assign row names and column names
dimnames(yinit.matrix)[[1]] <- agenames
dimnames(yinit.matrix)[[2]] <- StateNames

# Initializes population with infants under 3 months being protected by maternal immunity and with 1 infected person per age group in other age groups
yinit.matrix[,c('S1', 'I2', 'S2', 'I3', 'S3', 'I4')] = 0
yinit.matrix[,'M'] = c(Pop1[1:3], rep(0, N_ages-3))
yinit.matrix[,'S0'] = c(rep(0,3), Pop1[4:N_ages]-rep(N_ages-3))
yinit.matrix[,'I1'] = c(rep(0,3), rep(1, N_ages-3))

#Vectorize the states for ODE input
#Vectorize the ynit matrix
yinit.vector <- as.vector(yinit.matrix)

# Create array that has the labels by age and state, and use this to name the yinit.vector
name.array <- array(NA, dim = dim(yinit.matrix))
for(i in 1:dim(name.array)[1]){ # to number of rows of yinit.matrix
  for(j in 1:dim(name.array)[2]){ # to number of columns for yinit.matrix
    name.array[i,j] <- paste(dimnames(yinit.matrix)[[1]][i], # loop through the row names from agegrp1
                             dimnames(yinit.matrix)[[2]][j]) # loop through the column names from state M
  }
}
name.vector <- as.vector(name.array) # vectorise the name.array
names(yinit.vector) <- name.vector # assign labels to yinit.vector

# Contact matrix in use
n.cols = 100
nice.cols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(n.cols)
heatmap(c2/sum(diag(c2)), Rowv = NA, Colv = NA, scale = 'none', col = nice.cols)
