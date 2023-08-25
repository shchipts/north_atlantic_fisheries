
# Indicates whether specified process happens before or after relative event
# during transition from time t to time t + 1
#
# @param process Life-cycle process
#   (reproduction, maturation, mortality or growth)
# @param events Order of life-cycle processes
# @param relative_to  Relative life-cycle process
#   (reproduction, maturation, mortality or growth); growth, by default
# 
# @return 0 if specified process happens before relative event, and 1 if 
# specified process happens after relative event
# 
shift <- function(process, events, relative_to = "growth") {
  
  decode <- function(x) (
    switch(
      x,
      "reproduction" = "R",
      "maturation" = "D",
      "mortality" = "M",
      "growth" = "G"))
  
  label <- decode(process)
  label_relative_to <- decode(relative_to)
  
  order <- strsplit(events, "")[[1]]
  
  inc <- Reduce(
    function(seed, x) {
      
      if (seed == -1 & x == label_relative_to)
        return (1)
      
      if (seed == -1 & x == label)
        return (0)
      
      return (seed)
    },
    order,
    init = -1)
  
  return (inc)
}

# Juvenile's body length after growth
#
# @param init_length Initial body length
# @param gamma_1 Exponent of weight-length allometry
# @param gamma_2 Exponent of energy-weight allometry
# @param alpha_1 Coefficient of weight-length allometry
# @param alpha_2 Coefficient of energy-weight allometry
# 
length_juvenile <- function(
  init_length, 
  gamma_1, 
  gamma_2, 
  alpha_1, 
  alpha_2) {
  
  ((init_length ^ (gamma_1 * gamma_2) + 
       gamma_1 * alpha_1 * alpha_2 ^ (- gamma_1)) ^ (1 / (gamma_1 * gamma_2)))
}

# Somagonadic growth of adult body length
#
# @param init_body_length Body length before growth
# @param gamma_1 Exponent of weight-length allometry
# @param gamma_2 Exponent of energy-weight allometry
# @param alpha_1 Coefficient of weight-length allometry
# @param alpha_2 Coefficient of energy-weight allometry
# @param gsi Base gonadosomatic index
# 
somagonadic_growth <- function(
  init_body_length, 
  gamma_1, 
  gamma_2, 
  alpha_1, 
  alpha_2, 
  gsi) {
  
  ((init_body_length ^ (gamma_1 * gamma_2) + 
      gamma_1 * alpha_1 * alpha_2 ^ (- gamma_1)) /
      (1 + gamma_1 * gsi)) ^ (1 / (gamma_1 * gamma_2))
}

# Determines whether there is available surplus energy for adult's
# reproductive investment
#
# @param init_body_length Body length before growth
# @param gamma_1 Exponent of weight-length allometry
# @param gamma_2 Exponent of energy-weight allometry
# @param alpha_1 Coefficient of weight-length allometry
# @param alpha_2 Coefficient of energy-weight allometry
# @param gsi Base gonadosomatic index
# 
indicator_reproductive_investment <- function(
  init_body_length, 
  gamma_1, 
  gamma_2, 
  alpha_1, 
  alpha_2, 
  gsi) {
  
  l <- somagonadic_growth(
    init_body_length, 
    gamma_1, 
    gamma_2, 
    alpha_1, 
    alpha_2, 
    gsi)
  
  return (l >= init_body_length)
}
  
# Adult's body length after growth
#
# @param init_body_length Body length before growth
# @param gamma_1 Exponent of weight-length allometry
# @param gamma_2 Exponent of energy-weight allometry
# @param alpha_1 Coefficient of weight-length allometry
# @param alpha_2 Coefficient of energy-weight allometry
# @param gsi Base gonadosomatic index
# 
length_adult <- function(
  init_body_length, 
  gamma_1, 
  gamma_2, 
  alpha_1, 
  alpha_2, 
  gsi) {
  
  l <- somagonadic_growth(
    init_body_length, 
    gamma_1, 
    gamma_2, 
    alpha_1, 
    alpha_2, 
    gsi)
  
  return (if (l < init_body_length) init_body_length else l)
}

# Juvenile's maturation probability
#
# @param age Age
# @param body_length Body length
# @param steepness Inverse steepness of the logistic curve
# @param pmrn_slope Slope of probabilistic maturation reaction norm
# @param pmrn_intercept Intercept of probabilistic maturation reaction norm
#
maturation_probability <- function(
  age,
  body_length, 
  steepness,
  pmrn_slope,
  pmrn_intercept) {
  
  pmrn_midpoint <- pmrn_slope * age + pmrn_intercept
  
  return (1 / (1 + exp(- (body_length - pmrn_midpoint) / steepness)))
}

# Determines how steeply maturation probability changes around midpoint
#
# @param pmrn_width Width of probabilistic maturation reaction norm
# @param pmrn_envelope Probability at lower bound of maturation envelope
#
maturation_steepness <- function(
  pmrn_width,
  pmrn_envelope) {
  
  pmrn_width / (log((1- pmrn_envelope) / pmrn_envelope) - 
                  log(pmrn_envelope / (1 - pmrn_envelope)))
}

# Calculates gonadosomatic index based on length before and after growth
#
# @param body_length Body length after growth
# @param init_body_length Body length before growth
# @param gamma_1 Exponent of weight-length allometry
# @param gamma_2 Exponent of energy-weight allometry
# @param alpha_1 Coefficient of weight-length allometry
# @param alpha_2 Coefficient of energy-weight allometry
# 
gsi <- function(
  body_length,
  init_body_length,
  gamma_1, 
  gamma_2, 
  alpha_1, 
  alpha_2) {
  
  (init_body_length ^ (gamma_1 * gamma_2) - 
     body_length ^ (gamma_1 * gamma_2) +
     gamma_1 * alpha_1 * alpha_2 ^ (- gamma_1)) / 
    (gamma_1 * body_length ^ (gamma_1 * gamma_2))
}

# Fish body weight
#
# @param body_length Body length
# @param gamma_2 Exponent of energy-weight allometry
# @param alpha_2 Coefficient of energy-weight allometry
# @param body_length_ref Reference length
#
weight <- function(
  body_length,
  gamma_2, 
  alpha_2,
  body_length_ref) {
  
  return (alpha_2 * (body_length / body_length_ref) ^ gamma_2)
}

# Proportion of individuals which produce spawn out of the total number of 
# individuals at the spawning grounds
#
# @param p_spf Proportion of adults caught before they spawn out of the total 
# number of adults caught at the spawning ground
# @param F_spf Instantaneous fishing mortality due to spawner fisheries
#
spawner_proportion <- function(
  p_spf,
  F_spf) {
  
  return (1 - p_spf + p_spf * exp(- F_spf))
}

# Spawning stock biomass per fish
#
# @param body_weight Body weight at reproduction
# @param p_spf Proportion of adults caught before they spawn out of the total 
# number of adults caught at the spawning ground
# @param F_spf Instantaneous fishing mortality due to spawner fisheries
#
ssb_pf <- function(
  body_weight,
  p_spf,
  F_spf) {
  
  return (body_weight * spawner_proportion(p_spf, F_spf))
}

# Effective spawning stock biomass per fish
#
# @param body_weight Body weight at reproduction
# @param gsi Gonadosomatic index at reproduction
# @param gsi_base Base gonadosomatic index
#
ssb_pf_effective <- function(
  ssb_pf,
  gsi,
  gsi_base) {
  
  return (gsi / gsi_base * ssb_pf)
}

# Spawning stock
#
# @param n Number of spawners
# @param body_weight Body weight
# @param indicator_RI Indicates whether fish has invested in reproduction
# 
ssb <- function(
  n,
  body_weight,
  indicator_RI) {
  
  return (n * body_weight * indicator_RI)
}

# Spawner-juvenile production
#
# @param ssb Spawning stock biomass
# @param sjm Spawner-juvenile model specification
#
juveniles_age_1_production <- function(ssb, sjm) {
  
  model <- sjm %>% pull(spawner_recruit_model)
  a <- sjm %>% pull(a)
  Rp <- sjm %>% pull(Rp)
  
  n_juveniles_age_1 <- case_when(
    model == "independent" ~ a * ssb,
    model == "beverton_holt" ~ a * ssb / (1 + a * ssb / Rp),
    model == "ricker" ~ a * ssb * exp(-a * ssb / (Rp * exp(1))))
  
  return (n_juveniles_age_1)
}

# Survival probability for a fish
#
# @param natural_mortality Length-dependent natural mortality
# @param fishing_mortality Length-dependent fishing mortality
# rate
#
survival_probability <- function(
  natural_mortality,
  fishing_mortality) {
  
  return (exp(-(natural_mortality + fishing_mortality)))
}

# Probability of fish being harvested
#
# @param natural_mortality Length-dependent natural mortality
# @param fishing_mortality Length-dependent fishing mortality
# rate
#
catch_probability <- function(
  natural_mortality,
  fishing_mortality) {
  
  s <- survival_probability(natural_mortality, fishing_mortality)
  z <- natural_mortality + fishing_mortality
  
  catch <- fishing_mortality / z * (1 - s)
  
  return (catch)
}

# Probability of fish dying naturally
#
# @param natural_mortality Length-dependent natural mortality
# @param fishing_mortality Length-dependent fishing mortality
# rate
#
naturally_dead_probabiity <- function(
  natural_mortality,
  fishing_mortality) {
  
  s <- survival_probability(natural_mortality, fishing_mortality)
  z <- natural_mortality + fishing_mortality
  
  naturally_dead <- natural_mortality / z * (1 - s)
  
  return (naturally_dead)
}

# Adult's survival probability at the spawning ground
#
# @param F_spf Instantaneous fishing mortality due to spawner fisheries
# @param p_spm Proportion of adults which die naturally after spawning
#
survival_probability_spawning_ground <- function(
  F_spf,
  p_spm) {
  
  return (exp(-F_spf) * (1 - p_spm))
}

# Adult's catch probability at the spawning ground
#
# @param F_spf Instantaneous fishing mortality due to spawner fisheries
#
catch_probability_spawning_ground <- function(F_spf) {
  
  return (1 - exp(-F_spf))
}

# Length-dependent fishing mortality
#
# @param length Body length
# @param F0 Size-independent component (offset)
# @param F1 Maximum size-dependent component
# @param F2 Steepness of size-dependent component
# @param F3 Inflection point of size-dependent component
# @param F4 Width of the bell-shaped curve or steepness of size-dependent 
# component for the 2nd leg
# @param F5 Inflection point of size-dependent component for the 2nd leg
# @param F6 Magnitude of the decline for the 2nd leg
#
fishing_mortality <- function(
  body_length,
  F0,
  F1,
  F2,
  F3,
  F4,
  F5,
  F6) {
  
  x <- 
    if(F6 > 0)
      F0 + F1 /(1 + exp(-F2 * (body_length - F3))) - 
        F6 / (1 + exp(-F4 * (body_length - F5)))
    else
      F0 + F1 /(1 + exp(-F2 * (body_length - F3) - F4 * (body_length - F3) ^ 2))
  
  fishing_mortality <- ifelse(x > 0, x, 0)
  
  return (fishing_mortality)
}

# Length-dependent natural mortality
#
# @param length Body length
# @param gamma_3 Exponent of natural mortality-length allometry
# @param alpha_3 Coefficient of natural mortality-length allometry 
# @param body_length_ref Reference body length
#
natural_mortality <- function(
  body_length,
  gamma_3,
  alpha_3,
  body_length_ref) {
  
  x <- alpha_3 * (body_length / body_length_ref) ^ gamma_3
  
  natural_mortality <- ifelse(x > 0, x, 0)
  
  return (natural_mortality)
}