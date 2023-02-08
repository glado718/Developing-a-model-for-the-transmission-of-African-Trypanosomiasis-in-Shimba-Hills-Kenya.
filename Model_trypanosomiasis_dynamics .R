library(deSolve)

SIRSImodel <- function(t, x, parms) {
  with(as.list(c(parms, x)), {
    dScdt= - alphac * beta * (Sc * It/(Sc +Ic + Rc))
    dIcdt = alphac * beta * (Sc * It/(Sc + Ic + Rc)) - delta * Ic
    dRcdt = delta * Ic
    dStdt= bt* (St + It) - dt * St - alphat * beta * (St * Ic /(Sc + Ic+ Rc))
    dItdt = alphat * beta * (St* Ic/ (Sc + Ic + Rc))- dt * It
    res <- c(dScdt,dIcdt , dRcdt, dStdt, dItdt )
    list(res)
  })
}

## Parameters

pars  <- c(beta = 0.75, alphac = 0.46 , alphat = 0.025, delta = 0.002, 
           dt = 0.03, bt = 0.05)

## vector of timesteps
times  <- seq(0, 60, by = 1)

## Start values for steady state
y = xstart <- c(Sc= 477, Ic = 29, Rc = 8, St = 1116, It = 74)


## Solving
out <-  lsoda(y = xstart, time = times, func = SIRSImodel, parms = pars,) 
out = as.data.frame(out)
head(out)

#plot the solution for all compartments on a single plot.
plot(out$time, out$Sc,
     type = 'l', col = 'green', ylab = "Number of individuals in each compartment",
     xlab = "Time (days)", lwd = 3, ylim = c(0,3000))

lines(out$time, out$Ic, col = 'red', lwd = 3)

lines(out$time, out$Rc, col = 'magenta', lwd = 3)

lines(out$time, out$St, col = 'cyan', lwd = 3)

lines(out$time, out$It, col = 'blue', lwd = 3)

legend("topleft", legend = c("Sc","Ic", "Rc", "St", "It"), col = c("green","red", "magenta","cyan","blue"), lwd = 3)

#legend(55, 2600, legend=c("Sc", "Ic", "Rc", "St", "It"), 
       fill = c("green","red","magenta","cyan","blue")
#)
