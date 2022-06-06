#SI model for vector
# Step 1: write a function that computes the system of ODEs and returns a list containing the derivatives. The input parameters of this function are time `t`, the compartments as a vector `x`, and model parameters `parameters`. The outputs of the function are a list of derivatives. We name this function as `SIR.model`.
SI.model <- function(t, x, parameters){ 
  # Get the current number of individuals in each compartment.
  S <- x[1]
  I <- x[2]
  
  # Get the model parameters.
  gamma <- parameters$gamma
  
  # Calculate the total number of individuals in the system.
  # This is just all the compartments added together.
  N <- S + I 
  
  # Calculate the derivatives (the system of ODEs).
  dS <- - gamma*S
  dI <- gamma*S
  
  # Return the derivatives inside a list.
  derivatives <- list(c(dS, dI))
  return(derivatives)
}

#Step 2: define the initial conditions, parameters and times at which you wish to compute the solution.
# The initial parameters: how many people are in each compartment at the start of the model?
init <- c(S =  , I = 2534)

# The solution times: the sequence of times for which to compute the solution.
# This will be numbers 0 to 30.
time <- seq(from = 0, to = 30, by = 1)

#Estimating model parameters
#Time between contacts and time between recoveries
#tc = 4    # time between contacts in days
#tr = 30   # recovery time in days

#Estimating gamma and beta
#beta = 1 / tc      # contact rate in per day
#gamma = 1 / tr     # recovery rate in per day

# The model parameters: what is the transmission rate, and recovery *rate*
# (this is 1 / the infectious period).
pars <- list(gamma =  0.14)
pars
#Step 3: solve the system using the 'lsoda' function from the 'deSolve' package (try help(lsoda) for help).
# Load in the deSolve package.
library(deSolve)

# lsoda is a function with inputs: ODE system (as a function), initial conditions,
# parameters and times to compute the solution. The solution for the SIR model is
# stored in 'out'.
out <- lsoda(func = SI.model, y = init, parms = pars, times = time)

# Convert the solution into a data.frame. This makes plotting easier.
out <- as.data.frame(out)

#Step 4: plot the solution for all compartments on a single plot.
# First plot time (x-axis) against the number of susceptible (y-axis).
plot(out$time, out$S,
     type = 'l', col = 'green', ylab = "Number of individuals in each compartment",
     xlab = "Time (days)", lwd = 3, ylim = c(0, 2800))
 
# Next, we can add lines for time against infecteds, and time against recovered/removed.
lines(out$time, out$I, col = 'red', lwd = 3)

# Finally, we can add a legend so we know what lines correspond to what compartment.
legend("right", legend = c("S","I"), col = c("green","red"), lwd = 2)


#Estimating R0

