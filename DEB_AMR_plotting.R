

# TO DO:
# 1) Plots from Day et al., Growth rate of each strain (so the y-axis will have the state variables 'Pr and Ps' plotted with different colors), drug concentration (so the parameter 'dc') on the x-axis
# 2) Plots from Day et al., Growth rate of each strain (so the y-axis will have state variables 'Pr and Ps' plotted with different colors), Food level on the x-axis (so the parameters 'Food')
# 3) Current plot set up - with both strains on the y-axis (so the y-axis will have the state variables 'Pr and Ps' plotted with different colors) and time on the y-axis and then each panel is:
# drug concentration, which is the parameter 'dc' and Food (vertical panels, so the parameter 'Food')

# Let's start with those three figures.
# 5) Repeat plot # 1 but now with additional lines showing how varying the e or sigma of each strain (Ps and Pr)...that's a lot for a single plot so
# it could be a panel split by either strain (Ps and Pr) on the left/right column. 


## Key assumptions of the model
## Importantly, I think these make very important - and testable - hypotheses
# The drug-sensitive pathogen Subscript[Ps] is a better competitor (in the absence of drugs because it is): 
# (a) more efficient at converting resources (i.e., it has a higher conversion efficiency, e) or
# (b) steals fewer resources from the host (i.e., it has a lower attack rate, \[Sigma] and is thus, less virulent)
# These traits come at the cost of higher sensitivity to the drug, hence Subscript[Ps] also have higher death rates, Subscript[m, s] 

# The drug-resistant pathogen Subscript[Pr] is a worse competitor (in the absence of drugs because it is):
# (a) less efficient at converting resources (i.e., it has a lower conversion efficiency, e) or
# (b) steals more resources from the host (i.e., it has a higher attack rate, \[Sigma] and is thus, more virulent)
# These traits come at the cost of lower  sensitivity to the drug, hence Subscript[Pr]  also have lower death rates, Subscript[m, r] 









# make sure to clear the workspace too
rm(list = ls())

# load libraries
library(deSolve)
library(magrittr)
library(MASS)


# the model:

macro_deb <- function(t, y, p) {
  #tIn <- Sys.time() ## this records the time before a block of code so that it can be optimized  ##
  
  #### STATE VARIABLES ####
  C    <- y[1]                            ## colon biomass...in our case, this is what exactly? 
  Ps   <- y[2]                            ## pathogen biomass - drug sensitive
  Pr   <- y[3]                            ## pathogen biomass - drug resistant
  
  
  
  #### HOST PARAMETERS #### 
  f       <- as.numeric(p["f"])            ## ingestion rate
  rho     <- as.numeric(p["rho"])          ## resource use by the host
  

  #### PATHOGEN PARAMETERS ####
  ePs      <- as.numeric(p["ePs"])         ## biomass conversion efficiency - drug-sensitive pathogen 
  ePr      <- as.numeric(p["ePr"])         ## biomass conversion efficiency - drug-resistant pathogen   
  hc       <- as.numeric(p["hc"])          ## pathogen handling times [g]
  sigmacs  <- as.numeric(p["sigmacs"])     ## drug-sensitive pathogen attack rates [g pathogen^(-1) day^(-1)
  sigmacr  <- as.numeric(p["sigmacr"])     ## drug-resistant pathogen attack rates [g pathogen^(-1) day^(-1)
  
  m        <- as.numeric(p["m"])           ## pathogen background morality rate [day^(-1)]
  u        <- as.numeric(p["u"])           ## mutation probability from drug-sens. to drug-resis.

  InfDosePs <- as.numeric(p["InfDosPs"])   ## infectious dose - drug-sensitive pathogen
  InfDosePr <- as.numeric(p["InfDosePr"])  ## infectious dose - drug-resistant pathogen
  dc        <- as.numeric(p["dc"])         ## drug concentration
  

  #### COMPOUND PARAMETERS ####
  dr <- 1 - tanh(15*(dc - 0.30))                                       ## dose response curve
  
  
  
  #### RATE EQUATIONS ####
  dC   <-  f - rho*C - (sigmacs*C/(hc + C))*Ps - (sigmacr*C/(hc + C))*Pr                   ## dynamics of colon biomass
  dPs  <- ePs*(sigmacs*C/(hc + C))*Ps - m*Ps - u*Ps + u*Pr - dr*Ps                         ## dynamics of drug-sensitive pathogen
  dPr  <- ePr*(sigmacr*C/(hc + C))*Pr - m*Pr - u*Pr + u*Ps                                 ## dynamics of drug-resistant pathogen
  
  return(list(c(dC, dPs, dPr)))
}


#tOut <- Sys.time()                   ## record the time after a block of code
#print(tOut-tIn)                      ## prints how long the code took to run without stopping its execution





##############################################################################
########################    **SIMULATE THE ODE SYSTEM**    ###################
#############       **Pr strain: higher attack rate, sigma **    #############
##############################################################################


##############################################################################
############################## ***PARAMETERS*** ##############################
##############################################################################

##  initial conditions 
y0             <- c(C = 100, Ps = 10, Pr = 0) 
times          <- seq(0, 100)

pars <- c(
  ## HOST ENERGY BUDGET PARAMETERS
  f = 10,  rho = 0.1,
  
  ## PATHOGEN PARAMETERS
  ePs = 20, ePr = 20, hc = 1000, sigmacs = 0.8, sigmacr = 1.0, 
  u = 0.01, m = 1, dc = 150) #%>%

out <- as.data.frame(
  ode(
    func = macro_deb,
    y = y0,
    times = times,
    parms = pars
  )
)

par(mfrow = c(3,2))
par(las = 1)
with(out, plot(times, Ps, type='l'))
with(out, plot(times, Pr, type='l'))
with(out, plot(times, C, type='l'))




# Explore effects of varying drug concentration and calorie intake -----------------------------------------------------------
# Here, the drug-resistant strain has higher attack rate 
# **This means that the resistant strain is also more virulent**
# You can see from the plots, that standard practices of both higher food and higher drug concentrations
# both lead to higher density of the drug-resistant strain.
op <- par(mgp = c(2,1,0), mar = c(3,3,1,1), mfrow = c(3,3))

y0             <- c(C = 100, Ps = 10, Pr = 0) 
times          <- seq(0, 100)

dvals          <- c(0, 150, 300)               # drug concentration
fvals          <- c(2, 5, 10)                  # food availability


for (dc in dvals) {
  for (f in fvals) {
    pars <- c(
      ## HOST ENERGY BUDGET PARAMETERS
      f = f,  rho = 0.1,
      
      ## PATHOGEN PARAMETERS
      ePs = 20, ePr = 20, hc = 1000, sigmacs = 0.8, sigmacr = 1.0, 
      u = 0.01, m = 1, dc = dc)
    out <- as.data.frame(
      ode(
        func = macro_deb, 
        y = y0, 
        times = times,
        parms = pars
      )
    )
title <- bquote(list(f ==.(f), dc==.(dc)))
par(las = 1)
plot(Ps ~ time, data = out, type = "l", col = "darkcyan", main = title, ylim = c(0, 140))
points(Pr ~ time, data = out, type = "l", col = "darkmagenta", main = title)
legend("top", c("Ps","Pr"),
       lty = c(1,1),
       lwd = 2, bty = "n", 
       col = c("darkcyan","darkmagenta"),
       ncol = 2)
  }
}












##############################################################################
########################    **SIMULATE THE ODE SYSTEM**    ###################
#############       **Pr strain: lower conversion efficiency, e**      #############
##############################################################################


##############################################################################
############################## ***PARAMETERS*** ##############################
##############################################################################

##  initial conditions 
y0             <- c(C = 100, Ps = 10, Pr = 0) 
times          <- seq(0, 100)

pars2 <- c(
  ## HOST ENERGY BUDGET PARAMETERS
  f = 10,  rho = 0.1,
  
  ## PATHOGEN PARAMETERS
  ePs = 20, ePr = 18, hc = 1000, sigmacs = 0.8, sigmacr = 0.8, 
  u = 0.01, m = 1, dc = 150) #%>%

out <- as.data.frame(
  ode(
    func = macro_deb,
    y = y0,
    times = times,
    parms = pars2
  )
)

par(mfrow = c(3,2))
par(las = 1)
with(out, plot(times, Ps, type='l'))
with(out, plot(times, Pr, type='l'))
with(out, plot(times, C, type='l'))




# Explore effects of varying drug concentration and calorie intake -----------------------------------------------------------
# Here, the drug-resistant strain has LOWER conversion efficiency, e
# You can see from the plots, that standard practices of both higher food and higher drug concentrations
# both lead to higher density of the drug-resistant strain.
op <- par(mgp = c(2,1,0), mar = c(3,3,1,1), mfrow = c(3,3))

y0             <- c(C = 100, Ps = 10, Pr = 0) 
times          <- seq(0, 100)

dvals          <- c(0, 150, 300)               # drug concentration
fvals          <- c(2, 5, 10)                  # food availability


for (dc in dvals) {
  for (f in fvals) {
    pars2 <- c(
      ## HOST ENERGY BUDGET PARAMETERS
      f = f,  rho = 0.1,
      
      ## PATHOGEN PARAMETERS
      ePs = 20, ePr = 18, hc = 1000, sigmacs = 0.8, sigmacr = 0.8, 
      u = 0.01, m = 1, dc = dc)
    out <- as.data.frame(
      ode(
        func = macro_deb, 
        y = y0, 
        times = times,
        parms = pars2
      )
    )
    title <- bquote(list(f ==.(f), dc==.(dc)))
    par(las = 1)
    plot(Ps ~ time, data = out, type = "l", col = "darkcyan", main = title, ylim = c(0, 140))
    points(Pr ~ time, data = out, type = "l", col = "darkmagenta", main = title)
    legend("top", c("Ps","Pr"),
           lty = c(1,1),
           lwd = 2, bty = "n", 
           col = c("darkcyan","darkmagenta"),
           ncol = 2)
  }
}



















  ##############################################################################
  ########################    **SIMULATE THE ODE SYSTEM**    ###################
  #############       **Pr strain: lower conversion efficiency, e**         ####
  #############       **Pr strain: higher attack rate, sigma**              ####
  #############       **Pr strain: more virulent                            ####
  ##############################################################################
  
  
  ##############################################################################
  ############################## ***PARAMETERS*** ##############################
  ##############################################################################
  
  ##  initial conditions 
  y0             <- c(C = 100, Ps = 10, Pr = 0) 
  times          <- seq(0, 100)
  
  pars2 <- c(
    ## HOST ENERGY BUDGET PARAMETERS
    f = 10,  rho = 0.1,
    
    ## PATHOGEN PARAMETERS
    ePs = 20, ePr = 18, hc = 1000, sigmacs = 0.8, sigmacr = 1.0, 
    u = 0.01, m = 1, dc = 150) #%>%
  
  out <- as.data.frame(
    ode(
      func = macro_deb,
      y = y0,
      times = times,
      parms = pars2
    )
  )
  
  par(mfrow = c(3,2))
  par(las = 1)
  with(out, plot(times, Ps, type='l'))
  with(out, plot(times, Pr, type='l'))
  with(out, plot(times, C, type='l'))
  
  
  
  
  # Explore effects of varying drug concentration and calorie intake -----------------------------------------------------------
  # Here, the drug-resistant strain has LOWER conversion efficiency, e
  # You can see from the plots, that standard practices of both higher food and higher drug concentrations
  # both lead to higher density of the drug-resistant strain.
  op <- par(mgp = c(2,1,0), mar = c(3,3,1,1), mfrow = c(3,3))
  
  y0             <- c(C = 100, Ps = 10, Pr = 0) 
  times          <- seq(0, 100)
  
  dvals          <- c(0, 150, 300)               # drug concentration
  fvals          <- c(2, 5, 10)                  # food availability
  
  
  for (dc in dvals) {
    for (f in fvals) {
      pars2 <- c(
        ## HOST ENERGY BUDGET PARAMETERS
        f = f,  rho = 0.1,
        
        ## PATHOGEN PARAMETERS
        ePs = 20, ePr = 18, hc = 1000, sigmacs = 0.8, sigmacr = 1.0, 
        u = 0.01, m = 1, dc = dc)
      out <- as.data.frame(
        ode(
          func = macro_deb, 
          y = y0, 
          times = times,
          parms = pars2
        )
      )
      title <- bquote(list(f ==.(f), dc==.(dc)))
      par(las = 1)
      plot(Ps ~ time, data = out, type = "l", col = "darkcyan", main = title, ylim = c(0, 140))
      points(Pr ~ time, data = out, type = "l", col = "darkmagenta", main = title)
      legend("top", c("Ps","Pr"),
             lty = c(1,1),
             lwd = 2, bty = "n", 
             col = c("darkcyan","darkmagenta"),
             ncol = 2)
    }
  }
  
  
  
  
  
  
  
  #### This isn't quite right....
  #### I want drug concentration on the x-axis 
  #### I can't figure out how to average over the attractor
  
  Cnew <- Psnew <- Prnew <- c()                   # initialize counters for storing data
  
  
  for (dcvals in 100:300) {
    pars <- c(f = 10, rho = 0.1, ePs = 20, ePr = 18, hc = 1000, 
              sigmacs = 0.8, sigmacr = 1.0, 
              u = 0.01, m = 1, dc = dcvals)
    set.seed(1234)
    seeds <- floor(runif(24, 1, 1e5))
    mclapply(seeds,
             function(s) deb_mod(100, tstep = 0.1, pars = pars, y0 = c(C = 100, Ps = 10, Pr = 0), seed = s),
             mc.cores = 12) -> sim
    
    Cnew <- cbind(Cnew,
                      lapply(sim, function(s) s[, 2]) %>% do.call(cbind.data.frame, .) %>% apply(., 1, mean))
    Psnew <- cbind(Psnew,
                       lapply(sim, function(s) s[, 3]) %>% do.call(cbind.data.frame, .) %>% apply(., 1, mean))
    Prnew <- cbind(Prnew,
                       lapply(sim, function(s) s[, 4]) %>% do.call(cbind.data.frame, .) %>% apply(., 1, mean))
  }
  
  par(mfrow = c(3, 1), mar = c(4, 4, 0, 0), oma = rep(0.5, 4))
  plot.new()
  plot.window(xlim = c(0, 100), ylim = range(Cnew[, 1:3]))
  axis(1); axis(2); box('plot')
  for(i in 1:3) lines(seq(0, 100, 0.1), Cnew[,i], col = i, lwd = 2)
  mtext(side = 1, line = 2, "Time")
  mtext(side = 2, line = 2, "Calorie intake")
  
  plot.new()
  plot.window(xlim = c(0, 100), ylim = range(Psnew[, 1:3]))
  axis(1); axis(2); box('plot')
  for(i in 1:3) lines(seq(0, 100, 0.1), Psnew[,i], col = i, lwd = 2)
  mtext(side = 1, line = 2, "Time")
  mtext(side = 2, line = 2, "Drug-sensitive abundance")
  
  plot.new()
  plot.window(xlim = c(0, 100), ylim = range(Prnew[, 1:3]))
  axis(1); axis(2); box('plot')
  for(i in 1:3) lines(seq(0, 100, 0.1), Prnew[,i], col = i, lwd = 2)
  mtext(side = 1, line = 2, "Time")
  mtext(side = 2, line = 2, "Drug-resistant abundance")
  
  