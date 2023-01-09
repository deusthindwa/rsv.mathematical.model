
# Source function of the transmission dynamic model
simple_model <- function(t, y, parms, time.step = 'month'){
  
  # Read in initial states and their names
  States <- array(y, dim = dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  # Unify the time unit of parameter inputs
  if(parms$time.step == 'month'){
    period = 12
    length.step = 30.44 #days
  }else if(parms$time.step == 'week'){
    period = 52.1775
    length.step = 7 #days
  }
  
  # Here we need to convert all the rates from 1/days to 1/length.step
  
  # Waning rate of maternal immunity (by time step)
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  
  # Aging rate (by time step)
  mu = 1/parms$WidthAgeClassMonth
  if(parms$time.step == 'week'){
    mu = 1/(WidthAgeClassMonth*4.345)
  }
  
  # rate of recovery of first, second and subsequent infection(s)
  gamma1 = 1/(parms$dur.days1/length.step)
  gamma2 = 1/(parms$dur.days2/length.step)  
  gamma3 = 1/(parms$dur.days3/length.step)  
  gamma4 = gamma3
  
  # Relative risk of second, third and subsequent infection(s)
  sigma1 = parms$sigma1
  sigma2 = parms$sigma2
  sigma3 = parms$sigma3
  
  # Relative infectiousness of second and subsequent infection(s)
  rho1 = parms$rho1
  rho2 = parms$rho2
  
  # Pull out the states  for the model as vectors
  M <-  States[,'M'] # protected by maternal immunity
  S0 <-  States[,'S0'] # fully susceptible population
  I1 <-  States[,'I1'] # first time infection (or infectiousness since exposed state is ignore -> very short)
  
  S1 <-  States[,'S1'] # susceptible population with build-up immunity after first time infection
  I2 <-  States[,'I2'] # second time infection/infectiousness
  
  S2 <-  States[,'S2'] # susceptible population with lower risk of re-infection after two previous infections
  I3 <-  States[,'I3'] # third time infection/infectiousness
  
  S3 <-  States[,'S3'] # susceptible population with lowest risk of re-infection after three previous infections
  I4 <-  States[,'I4'] # subsequent time infection/infectiousness
  
  # the number of age groups
  N.ages <- length(M)
  
  #parameters related to the force of infection
  
  # per capita transmission probability
  baseline.txn.rate = parms$baseline.txn.rate
  
  # transmission probability per unit time
  b <- baseline.txn.rate / (parms$dur.days1/length.step)
  
  # q depends on transmission type (density- or frequency-0dependent)
  q = parms$q
  
  # (whether depends on population density or not)
  contact = parms$contact # c2 is the contact matrix
  
  # transmission probability per unit time in each age group
  # 100 is a scaling factor for the contact matrix we choose
  beta <-  (b/100)/(sum(yinit.matrix)^(1-q))*contact
  
  Amp = parms$Amp # seasonal amplitude
  phi = parms$phi # seasonal phase shift
  seasonal.txn <- (1+Amp*cos((2*pi*t-phi*period)/period))
 
  # Seasonal transmission probability for a frequency-dependent transmission mode
  beta_a_i <- seasonal.txn * beta/sum(States)
  infectiousN <- I1 + rho1*I2 + rho2*I3 + rho2*I4
  
  # Force of infection plus vectorize force of infection
  lambda <- infectiousN %*% beta_a_i
  lambda <- as.vector(lambda)

  # create a matrix to record the changing variables
  dy <- matrix(NA, nrow = N_ages, ncol = ncol(States))
  colnames(dy) <- colnames(States)
  
  # Get period birth rate from annual birth rate
  # See the following page for birth rate calculation
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period # B = annual birth rate
  
  # The 'um' is death rate
  um = parms$um
  
  #mu represents aging to the next class
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  # Ordinary differential equations
  dy[,'M'] <- period.birth.rate*sum(States) - (omega+(mu+um))*M + Aging.Prop*c(0,M[1:(N.ages-1)]) 
  
  dy[,'S0'] <- omega*M - lambda*S0 - (mu + um)*S0 + Aging.Prop*c(0,S0[1:(N.ages-1)]) 
  
  dy[,'I1'] <-   lambda*S0 - (gamma1 + mu + um)*I1 + Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - sigma1*lambda*S1 - (mu+um)*S1 + Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- sigma1*lambda*S1 - gamma2*I2-(mu + um)*I2 + Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - sigma2*lambda*S2 -(mu+um)*S2 + Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- sigma2*lambda*S2 -(gamma3 + mu+um)*I3 +  Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  gamma4*I4 -sigma3*lambda*S3 - (mu + um)*S3 + Aging.Prop*c(0,S3[1:(N.ages-1)]) 
  
  dy[,'I4'] <- sigma3*lambda*S3 - gamma4*I4 - (mu + um)*I4 + Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}


# NOTE: time step here is in months, you need to adjust seasonality accordingly for the time step you choose
start_time = 1 # start date (years)
tmax = nrow(B) # end_time = 25 # end date (years)
my_times <- seq(start_time, tmax, by = 1) # gives a sequence from start to end in increments of 1

# y is the initial population in each state and age
# t is the time steps of evaluation
# func specify the ordinary differential equations
# parms are all the parameter inputs
#Using the ODE function, we get the results of the transmission dynamic model5.
results <- ode(y = yinit.vector, t = my_times, func = simple_model, parms = parms)

