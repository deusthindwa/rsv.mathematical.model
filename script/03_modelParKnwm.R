
# We first assign values to all required parameters, including rates and probabilities. Most of these are described in Table 2 of Pitzer et al, PLOS Pathogens3


# Assign values of parameters from known literature
DurationMatImmunityDays = 30 # omega = 1/duration of maternal immunity in days
dur.days1 = 10 # gamma1 = 1/duration of infectiousness in days after first infection
dur.days2 = 7 # gamma2 = 1/duration of infectiousness in days after second infection
dur.days3 = 5 # gamma3 = 1/duration of infectiousness in days after thrid and subsequent infections

Amp = 0.2 # seasonal amplitude
phi = 3.327749 # seasonal phase shift or offset (depends on the starting month; here 0 = peak @ July 1)

baseline.txn.rate = 9 # per capita transmission probability (this is entire infectious period, it needs to divide the length of infectious period

# In Density Dependent transmission, the contact rate (c) depends on the population density. 
# In Frequent Dependent transmission, the contact rate (c’) does not depend on the population density. (we assumed frequency-dependent)
c2 = c2 # contact rate
q = 1 # use this to switch between density (q=0) vs frequency-dependent (q=1) transmission

# Relative infectiousness for 2nd and subsequent infections
rho1 = 0.75
rho2 = 0.51

# Relative risk of infection following 1st, 2nd, 3rd+ infections
sigma1 = 0.76
sigma2 = 0.6
sigma3 = 0.4

# Proportion of first infections that are LRI (by age); second infections that are LRI; third infections that are LRI
delta1 = c(rep(0.40,3), rep(0.39,3), rep(0.21,3), rep(0.20,3), 0.16, rep(0.14,3), rep(0.05,N_ages-16))
delta2 = 0.5*delta1
delta3 = 0.7*delta2

# Proportion of first infection that are hospitalized; second infection that are hospitalized; subsequent infection that are hospitalized
# Proportion of infections hospitalized = hosp prob given LRI * LRI prob given infection
hosp1 = c(0.18*rep(.40,3), 0.08*rep(0.39,3), 0.07*rep(0.21,3), 0.06*rep(0.20,3), 0.06*0.16, 0.05*rep(0.14,3), 0.02*rep(0.05,N_ages-16))
hosp2 = 0.4*hosp1
hosp3 = c(rep(0,N_ages-2), 0.00001, 0.00004)

# Birth rate (births/person/YEAR)
# Matrix: T rows, N_ages columns; columns 2:N_ages all 0s
PerCapitaBirthsYear = B

# Net rate of crude deaths (+) and immigration (-) You should calibrate this parameter so we can reproduce the population growth
# What is time(-1) unit? is this right? only die from last age class--this maybe helps with this "We assumed deaths occurred from 
# the last age class and adjusted the net rate of immigration/emigration and death from other age groups in order to produce a rate 
# of population growth and age structure similar to that of the US." 
# From Giny's code: um = log(0.993)/52 #net rate of crude deaths (+) and immigration (-) from all age groups (per week): can adjust 
# this to approximate population growth in state#. Calibrated this parameter so we can reproduce the population growth
um= -0.0002227 #(from all age groups)

# Aging rate = 1/width age class (months). Vector of long N_age
WidthAgeClassMonth = c(rep(1, times = 12), rep(12, times = 4), 60, 120, 240, 240, 240)

# Save parameters in a list
parms<-list(PerCapitaBirthsYear = PerCapitaBirthsYear,
            DurationMatImmunityDays = DurationMatImmunityDays,
            WidthAgeClassMonth = WidthAgeClassMonth,
            um = um, # net growth rate
            Amp = Amp, # seasonal amplitude
            phi = phi, # seasonal peak timing
            rho1 = rho1, # Relative infectiousness (2nd)
            rho2 = rho2, # Relative infectiousness (3rd+)
            dur.days1 = dur.days1,# Duration of infectiousness
            dur.days2 = dur.days2,# Duration of infectiousness
            dur.days3 = dur.days3,# Duration of infectiousness
            yinit.matrix = yinit.matrix, # initial states
            baseline.txn.rate = baseline.txn.rate,
            q = q, # Frequency or Density dependent
            contact = c2, # contact matrix
            sigma1 = sigma1, # Relative risk of infection (2nd)
            sigma2 = sigma2, # Relative risk of infection (3rd)
            sigma3 = sigma3, # Relative risk of infection (4th+)
            time.step = 'month')
