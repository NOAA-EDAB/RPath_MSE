#Just the code
library(here); library(Rpath); library(data.table)

#Georges Bank - EMAX (modified)----
GB.groups <- c('Phytoplankton- Primary Producers', 'Bacteria', 'Microzooplankton',
               'Small copepods', 'Large Copepods', 'Gelatinous Zooplankton', 
               'Micronekton', 'Mesopelagics', 'Macrobenthos- polychaetes', 
               'Macrobenthos- crustaceans', 'Macrobenthos- molluscs', 
               'Macrobenthos- other', 'Megabenthos- filterers', 'Megabenthos- other',
               'Shrimp et al.', 'Larval-juv fish- all', 'Small Pelagics- commercial',
               'Small Pelagics- other', 'Small Pelagics- squid', 
               'Small Pelagics- anadromous', 'Medium Pelagics- (piscivores & other)',
               'Demersals- benthivores', 'Demersals- omnivores', 
               'Demersals- piscivores', 'Sharks- pelagics', 'HMS', 'Baleen Whales',
               'Odontocetes', 'Sea Birds', 'Discard', 'Detritus-POC', 'Fishery')
types <- c(1, rep(0, 28), rep(2, 2), 3)

GB.params <- create.rpath.params(GB.groups, types)

GB.params$model$Biomass <- c(25.70472, 6.517908, 5.587981, 12.98514, 6.980794,
                             1.319463, 3.805126, 0.045, 11.40272, 10.87353,
                             9.8865, 40.02257, 3.613779, 3.965064, 0.09, 
                             0.6293996, 14.97737, 1.0737, 1.262162, 0.25,
                             0.2915301, 4.576176, 3.438957, 2.244675, 0.04334066,
                             0.025, 0.4167178, 0.1127281, 0.003496863,
                             NA, NA, NA)

GB.params$model$PB <- c(166.1342, 91.24998, 72.00002, 41.66504, 54.63586, 
                        40, 14.25, 0.9503762, 2.5, 3, 2, 2, 5, 2, 2, 15, 
                        0.3452712, 0.9571092, 0.9503762, 0.4249809, 0.459,
                        0.45, 0.45, 0.486, 0.102, 0.02, 0.03802086,
                        0.04, 0.275, NA, NA, NA)

GB.params$model$QB <- c(NA, 38.02082, 24.24243, 12.775, 10.95, 14.308, 36.5,
                        1.825, 17.5, 21, 14, 17.64, 18, 18, 5, 45, 2, 2,
                        2.75, 2, 2.3814, 0.92, 0.83, 2.205567, 0.5328019,
                        2.053014, 4.5, 13.82976, 4.379231, NA, NA, NA)

GB.params$model$BioAcc <- c(rep(0, 31), NA)

GB.params$model$Unassim <- c(0, rep(0.2, 28), 0, 0, NA)

GB.params$model$Discard <- c(rep(0, 31), 1)

GB.params$model$'Detritus-POC' <- c(rep(1, 29), rep(0, 3))

GB.params$model$Fishery <- c(rep(0, 11), 0.6478018, 0.03305462, 0, 0.00022901,
                             0.2990043, 0, 0.003389353, 0.03275818, 0.01015249,
                             0.1, 0.005314156, 0.5314448, 0, 0.003406343, 
                             rep(0, 6), NA)

GB.params$model$Fishery.disc <- c(rep(0, 5), 6.36e-7, 0, 3.7e-12, rep(0, 4), .2,
                                  0, 0, 4.8e-9, 0.0299, 0, 3.39e-4, 3.28e-3, 
                                  3.05e-3, 7.92e-2, 1.59e-3, .159, 0, 
                                  0, 1.25e-8, 1.02e-8, 0, 0, 0, NA)

#Diet
GB.params$diet[, Bacteria := c(0.24, rep(NA, 29), 0.76, NA)]

GB.params$diet[, Microzooplankton := c(0.216, 0.16, 0.12, rep(NA, 27), 0.504, NA)]

GB.params$diet[, 'Small copepods' := c(0.724, NA, 0.08, 0.065, rep(NA, 26), 
                                       0.131, NA)]

GB.params$diet[, 'Large Copepods' := c(0.546, NA, 0.0439, 0.174, 0.122, 0.0531, 
                                       rep(NA, 3), 1.92e-4, NA, 1.02e-4, 
                                       rep(NA, 18), 0.0594, NA)]

GB.params$diet[, 'Gelatinous Zooplankton' := 
                 c(0.087, 0.02, 0.051, 0.335, 0.366, 0.021, rep(NA, 9), 0.01, 
                   0.005, 0.002, 0.00043, 0.000016, rep(NA, 10), 0.102, NA)]

GB.params$diet[, Micronekton := c(0.162, NA, NA, .308, .325, NA, .041, rep(NA, 23),
                                  0.163, NA)]

GB.params$diet[, Mesopelagics := c(0.028, 0.017, 0.072, 0.34, 0.524, NA, 0.014, 
                                   rep(NA, 23), 0.004, NA)]

GB.params$diet[, 'Macrobenthos- polychaetes' := 
                 c(0.128, 0.308, rep(NA, 6), 0.015, 9.8e-4, 5.79e-4, 3.54e-3, 
                   3.56e-3, 1.9e-4, rep(NA, 15), 0.00593, 0.534, NA)]

GB.params$diet[, 'Macrobenthos- crustaceans' :=
                 c(0.212, 0.187, NA, 0.019, 0.037, rep(NA, 3), 0.015, 0.008, 
                   0.008, 0.026, 0.018, 2.9e-4, rep(NA, 7), 8.5e-4, 5.1e-4, 
                   rep(NA, 6), 0.009, 0.459, NA)]

GB.params$diet[, 'Macrobenthos- molluscs' :=
                 c(0.432, 0.199, rep(NA, 8), 0.004, 0.003, 0.013, 4.3e-4, 
                   rep(NA, 15), 0.006, 0.342, NA)]

GB.params$diet[, 'Macrobenthos- other' := 
                 c(0.215, 0.231, rep(NA, 6), 0.017, 0.025, 0.016, 0.057, 0.006,
                   0.003, rep(NA, 7), 8.5e-4, 4.2e-4, 2.4e-7, rep(NA, 5), 0.01,
                   0.418, NA)]

GB.params$diet[, 'Megabenthos- filterers' := c(0.69, 0.08, rep(NA, 28), 0.23, NA)]

GB.params$diet[, 'Megabenthos- other' := c(NA, 0.176, rep(NA, 6), 0.078, 0.098,
                                           0.032, 0.294, 0.048, 0.045, rep(NA, 7),
                                           0.003, 0.002, 2.4e-7, rep(NA, 5), 0.048,
                                           0.176, NA)]

GB.params$diet[, 'Shrimp et al.' := c(0.062, 0.365, rep(NA, 4), 0.123, NA, NA, 
                                      0.009, NA, 0.013, NA, NA, 5.4e-4, 
                                      rep(NA, 14), 0.062, 0.365, NA)]

GB.params$diet[, 'Larval-juv fish- all' := c(0.062, NA, NA, 0.456, 0.264, NA, 
                                             0.072, NA, 0.01, 0.007, 0.005, 0.005,
                                             rep(NA, 3), 0.055, rep(NA, 14), 
                                             0.062, NA)]

GB.params$diet[, 'Small Pelagics- commercial' :=
                 c(0.0113, NA, NA, 0.15, 0.436, 0.0818, 0.165, NA, 0.00892, 
                   0.0247, 0.00892, 0.0113, rep(NA, 3), 0.0963, rep(NA, 5), 
                   9.06e-4, 9.06e-4, 0.00413, rep(NA, 8))]

GB.params$diet[, 'Small Pelagics- other' :=
                 c(0.158, NA, NA, 0.115, 0.573, 0.102, 0.042, NA, NA, 6.9e-4, 
                   4.8e-4, 2.4e-4, rep(NA, 3), 0.007, rep(NA, 14), 0.0007, NA)]

GB.params$diet[, 'Small Pelagics- squid' := 
                 c(rep(NA, 4), 0.129, NA, 0.456, NA, NA, 0.099, NA, 0.018, NA, 
                   NA, 0.011, 0.176, 0.016, 0.018, 0.077, 1e-4, rep(NA, 12))]

GB.params$diet[, 'Small Pelagics- anadromous' := 
                 c(0.012, NA, NA, 0.056, 0.9, NA, 0.02, NA, 4.9e-4, 0.002, 
                   rep(NA, 5), 0.008, rep(NA, 14), 0.001, NA)]

GB.params$diet[, 'Medium Pelagics- (piscivores & other)' :=
                 c(rep(NA, 5), 0.001, NA, 0.003, NA, 0.013, NA, 0.011, 0.003, 
                   0.018, 0.001, 0.002, 0.576, 0.044, 0.114, 0.01, 0.011, 0.098,
                   0.014, 0.079, rep(NA, 6), 0.001, NA)]

GB.params$diet[, 'Demersals- benthivores' :=
                 c(rep(NA, 5), 0.005, 0.001, NA, 0.111, 0.13, 0.111, 0.133, 0.104,
                   0.133, 0.006, NA, 0.11, 0.001, 0.01, NA, NA, 0.076, 0.029,
                   0.018, rep(NA, 5), 0.011, 0.011, NA)]

GB.params$diet[, 'Demersals- omnivores' :=
                 c(rep(NA, 5), 0.006, 0.028, NA, 0.132, 0.066, 0.066, 0.066, 0.072,
                   0.286, 0.01, 0.013, 0.132, 0.002, 0.022, NA, NA, 0.048, 0.004,
                   0.025, rep(NA, 5), 0.011, 0.011, NA)]

GB.params$diet[, 'Demersals- piscivores' :=
                 c(rep(NA, 5), 0.026, 0.001, 0.007, 0.013, 0.013, 0.02, 0.118,
                   0.246, 0.019, 0.012, 0.013, 0.322, 0.032, 0.015, 0.009, 0.001,
                   0.038, 0.007, 0.086, rep(NA, 6), 0.001, NA)]

GB.params$diet[, 'Sharks- pelagics' :=
                 c(rep(NA, 4), 0.031, 0.01, NA, 0.007, NA, 0.01, NA, 0.01, 
                   rep(NA, 4), 0.216, 0.082, 0.165, 0.001, 0.124, 0.051, 0.082,
                   0.072, 0.014, 0.01, 0.01, 0.021, 0.031, NA, 0.051, NA)]

GB.params$diet[, HMS := c(rep(NA, 5), 0.103, rep(NA, 10), 0.135, 0.741, 0.022, 
                          rep(NA, 13))]

GB.params$diet[, 'Baleen Whales' := c(rep(NA, 3), 0.059, 0.473, 0.001, 0.296, NA, 
                                      NA, 0.059, 0.006, 0.024, NA, 0.004, NA, NA,
                                      0.059, 0.004, 0.003, 1.3e-4, rep(NA, 10),
                                      0.012, NA)]

GB.params$diet[, Odontocetes := c(rep(NA, 5), 0.003, 0.027, rep(NA, 9), 0.379,
                                  0.205, 0.274, 0.002, 3.2e-4, NA, 0.068, 0.04,
                                  rep(NA, 3), 0.002, rep(NA, 4))]

GB.params$diet[, 'Sea Birds' := c(rep(NA, 4), 0.034, NA, 0.137, rep(NA, 7), 0.01,
                                  NA, 0.305, 0.263, 0.068, 0.001, 5.4e-4, NA, 
                                  0.041, 0.013, rep(NA, 5), 0.126, NA, NA)]

#DCs do not sum to one...mostly rounding errors will fix
#Several groups need .001 to balance - will add to most dominate prey
GB.params$diet[Group == 'Phytoplankton- Primary Producers', 
               c('Large Copepods', 'Macrobenthos- molluscs') := 
                 as.list(c(0.547, 0.433))]
GB.params$diet[Group == 'Large Copepods', 
               c('Gelatinous Zooplankton', 'Micronekton', 'Mesopelagics', 
                 'Small Pelagics- other', 'Small Pelagics- anadromous') :=  
                 as.list(c(0.367, 0.326, 0.525, 0.574, 0.901))]
GB.params$diet[Group == 'Detritus-POC', 'Macrobenthos- other' := 0.419]
GB.params$diet[Group == 'Small Pelagics- commercial', 
               c('Medium Pelagics- (piscivores & other)',
                 'Demersals- piscivores', 'Sea Birds') := 
                 as.list(c(0.577, 0.323, 0.306))]
#Larval fish and Sharks- pelagics need .002 so I added to top 2 prey
GB.params$diet[Group == 'Small copepods', 'Larval-juv fish- all' := 0.457]
GB.params$diet[Group == 'Large Copepods', 'Larval-juv fish- all' := 0.265]
GB.params$diet[Group == 'Small Pelagics- commercial', 'Sharks- pelagics' := 0.217]
GB.params$diet[Group == 'Small Pelagics- squid',      'Sharks- pelagics' := 0.166]
#HMS has .001 too much so removed from dominate prey
GB.params$diet[Group == 'Small Pelagics- other', 'HMS' := 0.74]

#Switch discard diet to import
na.diet <- GB.params$diet[Group == 'Import', 2:30]
GB.params$diet[Group == 'Import', 2:30] <- GB.params$diet[Group == 'Discard', 2:30]
GB.params$diet[Group == 'Discard', 2:30] <- na.diet

#Run model
GB <- rpath(GB.params, 'Georges Bank')

#Harvest Control Rule----
HCR <- function(bio.ref, bmsy, fmsy, trigger = 0.5, threshold = 0.25){
  #Calculate the shape of the curve (y = mx + b)
  m <- (0 - 1) / (threshold - trigger)
  b <- -1 * m * threshold
  
  #If at or above biomass trigger:fish at Fmsy
  if(bio.ref >= trigger * bmsy){
    target.f <- fmsy
  } else {
    #If between trigger and threshold biomass: Adjust F
    #y value is f/fmsy so need to multiple by fmsy to get f
    target.f <- (m * (bio.ref / bmsy) + b) * fmsy
  }
  #If below threshold: No fishing allowed
  if(bio.ref < threshold * bmsy) target.f <- 0
  
  return(target.f)
}

#Generate Data----
GenData <- function(Rsim.output, group, sim.year, Sigma = 0.3, bias = 1, freq = 1){
  # This routine generates data (unbiased generally) for the biological groups 
  #and returns the appropriate object
  
  node     <- extract.node(Rsim.output, group)
  TrueBio  <- node$AnnualBiomass[sim.year]
  Observed <- c()
  for(Iyear in seq_along(sim.year)){
    if (Iyear %% freq == 0){
      Observed[Iyear] <- TrueBio[Iyear] * exp(rnorm(1, 0, Sigma) - Sigma^2/2)
    } else Observed[Iyear] <- -1
  }
  Catch <- as.numeric(node$AnnualTotalCatch[sim.year])
  
  out <- list(ObsBio = Observed, TotCatch = Catch, Fmort = Catch / Observed)
  
  return(out)
} 

#Surplus Production model----
Prodmodel <- function(ObsBio, Catch){
  data <- data.table(ObsBio, Catch)
  data[, Bplus1  := shift(ObsBio, type = 'lead')]
  data[, Surplus := Bplus1 - ObsBio + Catch]
  data[, BB      := ObsBio ^ 2]
  data[, prod    := Bplus1 + Catch]
  
  sprod <- lm(Surplus ~ 0 + ObsBio + BB, data)
  r <- as.numeric(sprod$coeff[1])
  K <- -1 * as.numeric(r / sprod$coeff[2])
  MSY  <- r * K / 4
  Bmsy <- K / 2
  Fmsy <- r / 2
  out  <- list(MSY = MSY, Bmsy = Bmsy, Fmsy = Fmsy, r = r, K = K)
}

#IAV-----
# Interannual variation
iav <- function(effort){
  effort <- data.table(effort = effort)
  n <- nrow(effort)
  effort[, effort1 := shift(effort, type = 'lead')]
  effort[, var := (effort1 - effort)^2]
  output <- sqrt((1/(n - 1))*sum(effort[, var], na.rm = T)) / ((1/n)*sum(effort[, effort]))
}


#MSY params----
#Calculate b0
GB.scene <- rsim.scenario(GB, GB.params, 1:100)
GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 1:100, 
                           value = 0)
GB.b0 <- rsim.run(GB.scene, years = 1:100, method = 'AB')
plot.groups <- GB.groups[which(!GB.groups %in% c('Fishery', 'Discard'))]
rsim.plot(GB.b0, spname = plot.groups)

GB.omnivore <- extract.node(GB.b0, 'Demersals- omnivores')
omni.b0 <- max(GB.omnivore$AnnualBiomass)

#Calculate r and K
GB.scene <- rsim.scenario(GB, GB.params, 1:100)
GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 1:100, 
                           value = c(seq(1, 10, .5), seq(10, 0, -.5), rep(0, 20),
                                     seq(0, 10, .5), seq(10, 1, -.5)))
GB.range <- rsim.run(GB.scene, years = 1:100, method = 'AB')
rsim.plot(GB.range, spname = plot.groups)

#Omnivores
GB.omnivore <- extract.node(GB.range, 'Demersals- omnivores')
prod.data <- data.table(biomass = GB.omnivore$AnnualBiomass, 
                        catch   = GB.omnivore$AnnualTotalCatch)
prod.data[, Bplus1  := shift(biomass, type = 'lead')]
prod.data[, Surplus := Bplus1 - biomass + catch]
prod.data[, BB      := biomass ^ 2]
prod.data[, prod    := Bplus1 + catch]

sprod <- lm(Surplus ~ 0 + biomass + BB, prod.data)
prod.data$pred <- c(sprod$fitted.values, 0)
prod.data <- prod.data[1:99, ]
setkey(prod.data, biomass)

plot(prod.data$biomass, prod.data$Surplus)
lines(prod.data$biomass, prod.data$pred, col='red')

r <- as.numeric(sprod$coeff[1])
K <- -1 * as.numeric(r / sprod$coeff[2])
omni.MSY  <- r * K / 4
omni.Bmsy <- K / 2
omni.Fmsy <- r / 2

#Pelagics
GB.medpel <- extract.node(GB.range, 'Medium Pelagics- (piscivores & other)')
prod.data <- data.table(biomass = GB.medpel$AnnualBiomass, 
                        catch   = GB.medpel$AnnualTotalCatch)
prod.data[, Bplus1  := shift(biomass, type = 'lead')]
prod.data[, Surplus := Bplus1 - biomass + catch]
prod.data[, BB      := biomass ^ 2]
prod.data[, prod    := Bplus1 + catch]

sprod <- lm(Surplus ~ 0 + biomass + BB, prod.data)
prod.data$pred <- c(sprod$fitted.values, 0)
prod.data <- prod.data[1:99, ]
setkey(prod.data, biomass)

plot(prod.data$biomass, prod.data$Surplus)
lines(prod.data$biomass, prod.data$pred, col='red')

r <- as.numeric(sprod$coeff[1])
K <- -1 * as.numeric(r / sprod$coeff[2])
medpel.MSY  <- r * K / 4
medpel.Bmsy <- K / 2
medpel.Fmsy <- r / 2


#S1----
#Strategy 1 - maximize Demersals- omnivore
set.seed(123)

#Create initial scenario - burn in to genrate data to conduct initial assessment
GB.scene <- rsim.scenario(GB, GB.params, 1:100)
# GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 3:12, 
#                            value = c(rep(1.5, 2), rep(2, 2), rep(1.5, 2), 
#                                      rep(1, 2), rep(0.5, 2)))
GB.init <- rsim.run(GB.scene, years = 1, method = 'AB')

#Initial guess at Bmsy/Fmsy
Bmsy <- 3.0
Fmsy <- 0.3

#initialize outputs
S1.effort <- c()
S1.omni   <- c()
S1.medpel <- c()
S1.iav    <- c()
for(isim in 1:100){
  GB.full <- copy(GB.init)
  
  #track changes over time
  track.Bmsy <- Bmsy
  track.Fmsy <- Fmsy
  
  #Rpath standardizes Effort = 1 so Catch = qEB simplies to C = qB so q = C/B or F
  obs.data <- GenData(GB.full, 'Demersals- omnivores', 1)
  q <- obs.data$Fmort[1]
  
  #Initial guess is to double effort to achieve an F of .3
  newE <- 2
  #Set new Effort
  GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 2, 
                             value = newE)
  
  #simulate years 2 through 100 with assessment every 5
  for(Iyear in 2:100){
    GB.full <- rsim.step(GB.scene, GB.full, method = 'AB', Iyear)
    
    #Survey biomass
    new.data      <- GenData(GB.full, 'Demersals- omnivores', sim.year = Iyear)
    obs.data$ObsBio   <- c(obs.data$ObsBio,   new.data$ObsBio)
    obs.data$TotCatch <- c(obs.data$TotCatch, new.data$TotCatch)
    obs.data$Fmort    <- c(obs.data$Fmort,    new.data$Fmort)
    
    #Run assessment every 5 years
    if(Iyear %% 5 == 0){
      oldE <- newE
      assess <- Prodmodel(obs.data$ObsBio, obs.data$TotCatch)
      Bmsy <- assess$Bmsy
      Fmsy <- assess$Fmsy
      track.Bmsy <- c(track.Bmsy, Bmsy)
      track.Fmsy <- c(track.Fmsy, Fmsy)
      
      #Check HCR and calculate new Effort
      ftarget <- HCR(obs.data$ObsBio[Iyear], Bmsy, Fmsy)
      newE <- ftarget / q
    }
    
    #Set new Effort
    GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = Iyear + 1, 
                               value = newE)
  } 
  run <- copy(GB.full)
  run.effort <- data.table(effort = GB.scene$fishing$ForcedEffort[seq(1, 1200, 12), 2], run = isim)
  omnivore  <- extract.node(run, 'Demersals- omnivores')
  medpel    <- extract.node(run, 'Medium Pelagics- (piscivores & other)')
  run.omni   <- data.table(bio.obs = obs.data$ObsBio, bio.mod = omnivore$AnnualBiomass, 
                           land = omnivore$AnnualTotalCatch, run = isim)
  run.medpel <- data.table(bio.obs = NA,              bio.mod = medpel$AnnualBiomass,   
                           land = medpel$AnnualTotalCatch, run = isim)
  run.iav <- iav(omnivore$AnnualTotalCatch + medpel$AnnualTotalCatch)
  
  #Combine outputs
  S1.effort <- rbindlist(list(S1.effort, run.effort))
  S1.omni   <- rbindlist(list(S1.omni,   run.omni))
  S1.medpel <- rbindlist(list(S1.medpel, run.medpel))
  S1.iav    <- c(S1.iav, run.iav)
}

#S2----
#Strategy 2 - set effort based on either omnivores or pelagics
set.seed(123)

#Create initial scenario - burn in to genrate data to conduct initial assessment
GB.scene <- rsim.scenario(GB, GB.params, 1:100)
# GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 3:12, 
#                            value = c(rep(1.5, 2), rep(2, 2), rep(1.5, 2), 
#                                      rep(1, 2), rep(0.5, 2)))
GB.init <- rsim.run(GB.scene, years = 1, method = 'AB')

#initialize outputs
S2.effort <- c()
S2.omni   <- c()
S2.medpel <- c()
S2.iav    <- c()
for(isim in 1:100){
  GB.full <- copy(GB.init)
  
  #Initial guess at Bmsy/Fmsy
  omni.Bmsy   <- 3.0
  omni.Fmsy   <- 0.3
  medpel.Bmsy <- 0.5
  medpel.Fmsy <- 0.3
  
  #track changes over time
  track.Bmsy <- list(omni = omni.Bmsy, medpel = medpel.Bmsy)
  track.Fmsy <- list(omni = omni.Fmsy, medpel = medpel.Fmsy)
  
  #Rpath standardizes Effort = 1 so Catch = qEB simplies to C = qB so q = C/B or F
  omni.obs.data <- GenData(GB.full, 'Demersals- omnivores', 1)
  omni.q <- omni.obs.data$Fmort[1]
  medpel.obs.data <- GenData(GB.full, 'Medium Pelagics- (piscivores & other)', 1)
  medpel.q <- medpel.obs.data$Fmort[1]
  
  #Initial guess is to double effort to achieve an F of .3
  newE <- 2
  #Set new Effort
  GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 2, 
                             value = newE)
  
  #simulate years 2 through 100 with assessment every 5
  for(Iyear in 2:100){
    GB.full <- rsim.step(GB.scene, GB.full, method = 'AB', Iyear)
    
    #Survey biomass
    omni.new.data <- GenData(GB.full, 'Demersals- omnivores', sim.year = Iyear)
    omni.obs.data$ObsBio   <- c(omni.obs.data$ObsBio,   omni.new.data$ObsBio)
    omni.obs.data$TotCatch <- c(omni.obs.data$TotCatch, omni.new.data$TotCatch)
    omni.obs.data$Fmort    <- c(omni.obs.data$Fmort,    omni.new.data$Fmort)
    
    medpel.new.data <- GenData(GB.full, 'Medium Pelagics- (piscivores & other)', 
                               sim.year = Iyear)
    medpel.obs.data$ObsBio   <- c(medpel.obs.data$ObsBio,   medpel.new.data$ObsBio)
    medpel.obs.data$TotCatch <- c(medpel.obs.data$TotCatch, medpel.new.data$TotCatch)
    medpel.obs.data$Fmort    <- c(medpel.obs.data$Fmort,    medpel.new.data$Fmort)
    
    #Run assessment every 5 years
    if(Iyear %% 5 == 0){
      oldE <- newE
      omni.assess <- Prodmodel(omni.obs.data$ObsBio, omni.obs.data$TotCatch)
      omni.Bmsy <- omni.assess$Bmsy
      omni.Fmsy <- omni.assess$Fmsy
      track.Bmsy$omni <- c(track.Bmsy$omni, omni.Bmsy)
      track.Fmsy$omni <- c(track.Fmsy$omni, omni.Fmsy)
      
      medpel.assess <- Prodmodel(medpel.obs.data$ObsBio, medpel.obs.data$TotCatch)
      medpel.Bmsy <- medpel.assess$Bmsy
      medpel.Fmsy <- medpel.assess$Fmsy
      track.Bmsy$medpel <- c(track.Bmsy$medpel, medpel.Bmsy)
      track.Fmsy$medpel <- c(track.Fmsy$medpel, medpel.Fmsy)
      
      #Check HCR and calculate new Effort
      omni.ftarget <- HCR(omni.obs.data$ObsBio[Iyear], omni.Bmsy, omni.Fmsy)
      omniE <- omni.ftarget / omni.q
      
      medpel.ftarget <- HCR(medpel.obs.data$ObsBio[Iyear], medpel.Bmsy, medpel.Fmsy)
      medpelE <- medpel.ftarget / medpel.q
      
      newE <- min(omniE, medpelE)
      if(newE < 0) newE <- 0
    }
    
    #Set new Effort
    GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = Iyear + 1, 
                               value = newE)
  }
  run <- copy(GB.full)
  run.effort <- data.table(effort = GB.scene$fishing$ForcedEffort[seq(1, 1200, 12), 2], run = isim)
  omnivore  <- extract.node(run, 'Demersals- omnivores')
  medpel    <- extract.node(run, 'Medium Pelagics- (piscivores & other)')
  run.omni   <- data.table(bio.obs = omni.obs.data$ObsBio, bio.mod = omnivore$AnnualBiomass, 
                           land = omnivore$AnnualTotalCatch, run = isim)
  run.medpel <- data.table(bio.obs = medpel.obs.data$ObsBio, bio.mod = medpel$AnnualBiomass, 
                           land = medpel$AnnualTotalCatch, run = isim)
  run.iav <- iav(omnivore$AnnualTotalCatch + medpel$AnnualTotalCatch)
  
  #Combine outputs
  S2.effort <- rbindlist(list(S2.effort, run.effort))
  S2.omni   <- rbindlist(list(S2.omni,   run.omni))
  S2.medpel <- rbindlist(list(S2.medpel, run.medpel))
  S2.iav <- c(S2.iav, run.iav)
}

#S3----
#Strategy 3 - maximize omnivores with protection for pelagics
set.seed(123)

#Create initial scenario - burn in to genrate data to conduct initial assessment
GB.scene <- rsim.scenario(GB, GB.params, 1:100)
# GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 3:12, 
#                            value = c(rep(1.5, 2), rep(2, 2), rep(1.5, 2), 
#                                      rep(1, 2), rep(0.5, 2)))
GB.init <- rsim.run(GB.scene, years = 1, method = 'AB')

#initialize outputs
S3.effort <- c()
S3.omni   <- c()
S3.medpel <- c()
S3.iav    <- c()
for(isim in 1:100){
  GB.full <- copy(GB.init)
  
  #Initial guess at Bmsy/Fmsy
  Bmsy <- 3
  Fmsy <- .3
  #track changes over time
  track.Bmsy <- Bmsy
  track.Fmsy <- Fmsy
  
  #Rpath standardizes Effort = 1 so Catch = qEB simplies to C = qB so q = C/B or F
  obs.data <- GenData(GB.full, 'Demersals- omnivores', 1)
  q <- obs.data$Fmort[1]
  
  #Initial guess is to double effort to achieve an F of .3
  newE <- 2
  #Set new Effort
  GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = 2, 
                             value = newE)
  
  #Get initial Medium Pelagic biomass
  medpel.data <- GenData(GB.full, 'Medium Pelagics- (piscivores & other)', 
                         sim.year = 1)
  medpel.init.bio <- medpel.data$ObsBio
  medpel.obs <- medpel.init.bio
  
  #simulate years 2 through 100 with assessment every 5
  for(Iyear in 2:100){
    GB.full <- rsim.step(GB.scene, GB.full, method = 'AB', Iyear)
    
    #Survey biomass
    new.data      <- GenData(GB.full, 'Demersals- omnivores', sim.year = Iyear)
    obs.data$ObsBio   <- c(obs.data$ObsBio,   new.data$ObsBio)
    obs.data$TotCatch <- c(obs.data$TotCatch, new.data$TotCatch)
    obs.data$Fmort    <- c(obs.data$Fmort,    new.data$Fmort)
    medpel.obs <- c(medpel.obs, NA)
    
    #Run assessment every 5 years
    if(Iyear %% 5 == 0){
      oldE <- newE
      assess <- Prodmodel(obs.data$ObsBio, obs.data$TotCatch)
      Bmsy <- assess$Bmsy
      Fmsy <- assess$Fmsy
      track.Bmsy <- c(track.Bmsy, Bmsy)
      track.Fmsy <- c(track.Fmsy, Fmsy)
      
      #Check HCR and calculate new Effort
      ftarget <- HCR(obs.data$ObsBio[Iyear], Bmsy, Fmsy)
      newE <- ftarget / q
      
      #Add protection for Medium Pelagics - reduce effort by 40%
      medpel.data <- GenData(GB.full, 'Medium Pelagics- (piscivores & other)', 
                             sim.year = Iyear)
      medpel.obs[Iyear] <- medpel.data$ObsBio
      if(medpel.data$ObsBio >= 0.2 * medpel.init.bio & 
         medpel.data$ObsBio < 0.5 * medpel.init.bio) newE <- 0.6 * oldE
      if(medpel.data$ObsBio < 0.2 * medpel.init.bio) newE <- 0
      
      if(newE < 0) newE <- 0
    }
    
    #Set new Effort
    GB.scene <- adjust.fishing(GB.scene, 'ForcedEffort', 'Fishery', sim.year = Iyear + 1, 
                               value = newE)
  }
  run <- copy(GB.full)
  run.effort <- data.table(effort = GB.scene$fishing$ForcedEffort[seq(1, 1200, 12), 2], run = isim)
  omnivore  <- extract.node(run, 'Demersals- omnivores')
  medpel    <- extract.node(run, 'Medium Pelagics- (piscivores & other)')
  run.omni   <- data.table(bio.obs = obs.data$ObsBio, bio.mod = omnivore$AnnualBiomass, 
                           land = omnivore$AnnualTotalCatch, run = isim)
  run.medpel <- data.table(bio.obs = medpel.obs, bio.mod = medpel$AnnualBiomass,   
                           land = medpel$AnnualTotalCatch, run = isim)
  run.iav <- iav(omnivore$AnnualTotalCatch + medpel$AnnualTotalCatch)
  
  #Combine outputs
  S3.effort <- rbindlist(list(S3.effort, run.effort))
  S3.omni   <- rbindlist(list(S3.omni,   run.omni))
  S3.medpel <- rbindlist(list(S3.medpel, run.medpel))
  S3.iav <- c(S3.iav, run.iav)
}


#compile results----
#Average results
S1.effort.10 <- c()
S2.effort.10 <- c()
S3.effort.10 <- c()
for(irun in 1:10){
  effort.s1    <- matrix(S1.effort[run == irun, effort], 100, 1)
  S1.effort.10 <- cbind(S1.effort.10, effort.s1)
  effort.s2    <- matrix(S2.effort[run == irun, effort], 100, 1)
  S2.effort.10 <- cbind(S2.effort.10, effort.s2)
  effort.s3    <- matrix(S3.effort[run == irun, effort], 100, 1)
  S3.effort.10 <- cbind(S3.effort.10, effort.s3)
}

S1.biomass.omni   <- S1.omni[,   list(mean = mean(bio.mod), var = sd(bio.mod)), by = run]
S1.landing.omni   <- S1.omni[,   list(mean = mean(land),    var = sd(land)),    by = run]
S1.biomass.medpel <- S1.medpel[, list(mean = mean(bio.mod), var = sd(bio.mod)), by = run]
S1.landing.medpel <- S1.medpel[, list(mean = mean(land),    var = sd(land)),    by = run]

S2.biomass.omni   <- S2.omni[,   list(mean = mean(bio.mod), var = sd(bio.mod)), by = run]
S2.landing.omni   <- S2.omni[,   list(mean = mean(land),    var = sd(land)),    by = run]
S2.biomass.medpel <- S2.medpel[, list(mean = mean(bio.mod), var = sd(bio.mod)), by = run]
S2.landing.medpel <- S2.medpel[, list(mean = mean(land),    var = sd(land)),    by = run]

S3.biomass.omni   <- S3.omni[,   list(mean = mean(bio.mod), var = sd(bio.mod)), by = run]
S3.landing.omni   <- S3.omni[,   list(mean = mean(land),    var = sd(land)),    by = run]
S3.biomass.medpel <- S3.medpel[, list(mean = mean(bio.mod), var = sd(bio.mod)), by = run]
S3.landing.medpel <- S3.medpel[, list(mean = mean(land),    var = sd(land)),    by = run]

outplot <- data.table(omni.BBmsy      = c(mean(S1.biomass.omni[, mean]) / omni.Bmsy,
                                          mean(S2.biomass.omni[, mean]) / omni.Bmsy,
                                          mean(S3.biomass.omni[, mean]) / omni.Bmsy),
                      omni.BBmsy.sd   = c(sd(S1.biomass.omni[, mean]) / omni.Bmsy,
                                          sd(S2.biomass.omni[, mean]) / omni.Bmsy,
                                          sd(S3.biomass.omni[, mean]) / omni.Bmsy),
                      medpel.BBmsy    = c(mean(S1.biomass.medpel[, mean]) / medpel.Bmsy,
                                          mean(S2.biomass.medpel[, mean]) / medpel.Bmsy,
                                          mean(S3.biomass.medpel[, mean]) / medpel.Bmsy),
                      medpel.BBmsy.sd = c(sd(S1.biomass.medpel[, mean]) / medpel.Bmsy,
                                          sd(S2.biomass.medpel[, mean]) / medpel.Bmsy,
                                          sd(S3.biomass.medpel[, mean]) / medpel.Bmsy),
                      omni.CMSY       = c(mean(S1.landing.omni[, mean]) / omni.MSY,
                                          mean(S2.landing.omni[, mean]) / omni.MSY,
                                          mean(S3.landing.omni[, mean]) / omni.MSY),
                      omni.CMSY.sd    = c(sd(S1.landing.omni[, mean]) / omni.MSY,
                                          sd(S2.landing.omni[, mean]) / omni.MSY,
                                          sd(S3.landing.omni[, mean]) / omni.MSY),
                      medpel.CMSY     = c(mean(S1.landing.medpel[, mean]) / medpel.MSY,
                                          mean(S2.landing.medpel[, mean]) / medpel.MSY,
                                          mean(S3.landing.medpel[, mean]) / medpel.MSY),
                      medpel.CMSY.sd  = c(sd(S1.landing.medpel[, mean]) / medpel.MSY,
                                          sd(S2.landing.medpel[, mean]) / medpel.MSY,
                                          sd(S3.landing.medpel[, mean]) / medpel.MSY))

#Example output----
par(mfrow = c(3, 2), mar = c(0,0,0,0), oma = c(4, 8, 2, 8))

#S1
plot(S1.omni[run == 1, bio.mod], typ = 'l', col = 'red', ylim = c(0, 5), axes = F)
points(S1.omni[run == 1, bio.obs])
axis(2, las = T)
box()
mtext(3, text = 'Target Species')

plot(S1.medpel[run == 1, bio.mod], typ = 'l', col = 'red', ylim = c(0, 1), axes = F)
axis(4, las = T)
box()
mtext(3, text = 'Choke Species')
legend('topright', legend = c('Model', 'Observations', 'Threshold'), bty = 'n', 
       lty = c(1, 0, 2), col = c('red', 'black', 'grey'), pch = c(NA, 1, NA))

#S2
plot(S2.omni[run == 1, bio.mod], typ = 'l', col = 'red', ylim = c(0, 5), axes = F)
points(S2.omni[run == 1, bio.obs])
axis(2, las = T)
axis(3, lab = F, tick = -.02)
box()

plot(S2.medpel[run == 1, bio.mod], typ = 'l', col = 'red', ylim = c(0, 1), axes = F)
points(S2.medpel[run == 1, bio.obs])
axis(4, las = T)
axis(3, lab = F, tick = -.02)
box()

#S3
plot(S3.omni[run == 1, bio.mod], typ = 'l', col = 'red', ylim = c(0, 5), axes = F)
points(S3.omni[run == 1, bio.obs])
axis(2, las = T)
axis(3, lab = F, tick = -.02)
axis(1)
box()

plot(S3.medpel[run == 1, bio.mod], typ = 'l', col = 'red', ylim = c(0, 1), axes = F)
points(S3.medpel[run == 1, bio.obs])
abline(h = medpel.init.bio * .5, lty = 2, col = 'grey')
abline(h = medpel.init.bio * 0.2, lty = 2, col = 'grey')
axis(4, las = T)
axis(3, lab = F, tick = -.02)
axis(1)
box()

mtext(2, text = expression('Biomass t km'^-2), outer = T, line = 3.5)
mtext(1, text = 'Simulated Year', outer = T, line = 2.5)

#Effort - Figure 3--------------------------------------------------------------
s.colors <- c('#66c2a5', '#fc8d62', '#8da0cb')

par(mfrow = c(1, 3), mar = c(0, 0, 0, 0), oma = c(4, 7, 2, 2))

#Strategy 1
plot(0, 0, ylim = c(0, 4.5), xlim = c(0, 100), typ = 'n', axes = F, xlab = '', 
     ylab = '')
lines(S1.effort[run == 1, effort], lwd = 3, col = s.colors[1])
box()
x.label <- c('', axTicks(1)[2:length(axTicks(1))])
axis(1, cex.axis = 2, at = axTicks(1), labels = x.label)
axis(2, las= T, cex.axis = 2)
text(80, 4.4, labels = 'S1', cex = 2.5)

#Strategy 2
plot(0, 0, ylim = c(0, 4.5), xlim = c(0, 100), typ = 'n', axes = F, xlab = '', 
     ylab = '')
lines(S2.effort[run == 1, effort], lwd = 3, col = s.colors[2])
box()
axis(1, cex.axis = 2, at = axTicks(1), labels = x.label)
text(80, 4.4, labels = 'S2', cex = 2.5)

#Strategy 3
plot(0, 0, ylim = c(0, 4.5), xlim = c(0, 100), typ = 'n', axes = F, xlab = '', 
     ylab = '')
lines(S3.effort[run == 1, effort], lwd = 3, col = s.colors[3])
box()
axis(1, cex.axis = 2, at = axTicks(1), labels = x.label)
text(80, 4.4, labels = 'S3', cex = 2.5)

mtext(1, text = 'Simulation Year', line = 3, outer = T, cex = 1.5)
mtext(2, text = 'Relative Effort', line = 3, outer = T, cex = 1.5)

#Within run variability----
opar <- par(mfrow = c(3, 1), mar = c(0, 0, 0, 0), oma = c(4, 6, 2, 2))
boxplot(S1.effort.10, col = s.colors[1], axes = F, range = 0)
axis(2, las = T, cex.axis = 1.8)
box()
boxplot(S2.effort.10, col = s.colors[2], ylim = c(0, 6.3), axes = F, range = 0)
axis(2, las = T, cex.axis = 1.8)
box()
boxplot(S3.effort.10, col = s.colors[3], axes = F, range = 0)
axis(2, las = T, cex.axis = 1.8)
box()
mtext(2, text = 'Relative Effort', line = 3, cex = 2, outer = T)

#Average biomass----
#B/Bmsy and F/MSY
par(mfrow = c(2,2), mar = c(0, 0, 0, 0), oma = c(2, 7, 2, 2))

#B/Bmsy - Target
plot(0, 0, ylim = c(0, 2.0), xlim = c(0, 4), typ = 'n', axes = F, xlab = '', ylab = '')
abline(h = 1, lty = 4, col = 'red')
points(outplot[, omni.BBmsy], pch = 16, cex = 3, col = s.colors)
for(i in 1:3){
  arrows(i, outplot[i, omni.BBmsy - omni.BBmsy.sd], y1 = outplot[i, omni.BBmsy + omni.BBmsy.sd], 
         length = 0.2, angle = 90, code = 3, col = s.colors[i], lwd = 3)
}
axis(2, las = T, cex.axis = 1.5)
box()
mtext(2, text = expression('B/B'['MSY']), line = 4, cex = 1.5)
mtext(3, text = 'Target Species', cex = 1.5)

#B/Bmsy - Choke
plot(0, 0, ylim = c(0, 2.0), xlim = c(0, 4), typ = 'n', axes = F, xlab = '', ylab = '')
abline(h = 1, lty = 4, col = 'red')
points(outplot[, medpel.BBmsy], pch = 16, cex = 3, col = s.colors)
for(i in 1:3){
  arrows(i, outplot[i, medpel.BBmsy - medpel.BBmsy.sd], y1 = outplot[i, medpel.BBmsy + medpel.BBmsy.sd], 
         length = 0.2, angle = 90, code = 3, col = s.colors[i], lwd = 3)
}
box()
mtext(3, text = 'Choke Species', cex = 1.5)

#Catch/MSY - Target
plot(0, 0, ylim = c(0, 2.0), xlim = c(0, 4), typ = 'n', axes = F, xlab = '', ylab = '')
abline(h = 1, lty = 4, col = 'red')
points(outplot[, omni.CMSY], pch = 16, cex = 3, col = s.colors)
for(i in 1:3){
  arrows(i, outplot[i, omni.CMSY - omni.CMSY.sd], y1 = outplot[i, omni.CMSY + omni.CMSY.sd], 
         length = 0.2, angle = 90, code = 3, col = s.colors[i], lwd = 3)
}
axis(2, las = T, cex.axis = 1.5)
box()
mtext(2, text = 'Catch/MSY', line = 4, cex = 1.5)


#Catch/MSY - Choke
plot(0, 0, ylim = c(0, 2.0), xlim = c(0, 4), typ = 'n', axes = F, xlab = '', ylab = '')
abline(h = 1, lty = 4, col = 'red')
points(outplot[, medpel.CMSY], pch = 16, cex = 3, col = s.colors)
for(i in 1:3){
  arrows(i, outplot[i, medpel.CMSY - medpel.CMSY.sd], y1 = outplot[i, medpel.CMSY + medpel.CMSY.sd], 
         length = 0.2, angle = 90, code = 3, col = s.colors[i], lwd = 3)
}
box()

#IAV----
opar <- par(mfrow = c(1, 3), mar = c(0, 0, 0, 0), oma = c(4, 10, 2, 2))
boxplot(S1.iav, col = s.colors[1], ylim = c(0, 2.5), axes = F, cex = 2)
axis(2, las = T, cex.axis = 2)
box()
boxplot(S2.iav, col = s.colors[2], ylim = c(0, 2.5), axes = F, cex = 2)
axis(2, labels = F, las = T, cex.axis = 2)
box()
boxplot(S3.iav, col = s.colors[3], ylim = c(0, 2.5), axes = F, cex = 2)
axis(2, labels = F, las = T, cex.axis = 2)
box()
mtext(2, text = 'Interannual Variation (IAV)\n of the Catch', line = 3.8, cex = 2, outer = T)

