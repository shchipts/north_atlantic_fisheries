
# Scaling factors for recruitment level to obtain juvenile and adult abundances
juveniles_adults_in_equilibrium <- function(data_transitions) {
  
  juveniles <- data_transitions %>% filter(is.na(age_maturation))
  juveniles_to_adults <- data_transitions %>% filter(age == age_maturation)
  
  juveniles_coef_eq <- Reduce(
    function(seed, coef) {
      
      prev <- tail(seed, 1)
      
      return (append(seed, coef * prev))
    },
    juveniles %>% pull(coef),
    init = c(1)) %>%
    tail(-1)
  
  juveniles_to_adults_coef_eq <- (juveniles_to_adults %>% pull(coef)) *
    append(c(1), head(juveniles_coef_eq, -1)) 
  
  juveniles_eq <- juveniles %>%
    rbind(juveniles_to_adults) %>%
    select(age, age_maturation) %>% 
    cbind(coef = append(
      juveniles_coef_eq, 
      juveniles_to_adults_coef_eq))
  
  juveniles_to_adults_eq <- juveniles_eq %>% filter(!is.na(age_maturation))
  
  adults_eq <- Reduce(
    function(seed, a_m) {
      
      adults_selected <- data_transitions %>% 
        filter(age_maturation == a_m & age > age_maturation) %>%
        arrange(age)
      
      
      adults_to_adults_same_maturation_age_coef_eq <- Reduce(
        function(seed2, idx) {
          
          coef <- adults_selected %>% slice(idx) %>% pull(coef)
          prev <- ifelse(
            length(seed2) == 0,
            juveniles_to_adults_eq %>% 
              filter(age_maturation == a_m) %>%
              pull(coef),
            tail(seed2, 1))
          
          return (append(seed2, coef * prev))
        },
        1 : nrow(adults_selected),
        init = c())
      
      return (rbind(
        seed,
        adults_selected %>%
          select(age, age_maturation) %>% 
          cbind(coef = adults_to_adults_same_maturation_age_coef_eq)))
    },
    juveniles_to_adults_eq %>% pull(age_maturation),
    init = data.frame())
  
  
  return (rbind(juveniles_eq, adults_eq) %>%
            arrange(!is.na(age_maturation), age_maturation, age))
}

# Recruitment level at demographic equilibrium
recruits_nonlinear <- function(
  stock, 
  fishing, 
  base_gsi,
  life_history,
  scaler) {
  
  ssb_scaler <- scaler %>%
    segmented_SSB_pf_in_equilibrium(stock, life_history, fishing) %>%
    effective_SSB_equilibrium(stock, base_gsi)
  
  model <- stock %>% pull(spawner_recruit_model)
  a <- stock %>% pull(a) * ssb_scaler
  Rp <- stock %>% pull(Rp)
  
  if (a < 1)
    return (NA)
  else {
  
    recruits_eq <- case_when(
      model == "beverton_holt" ~ Rp * (a - 1) / a,
      model == "ricker" ~ Rp * exp(1) * log(a) / a)
    
    return (recruits_eq)
  }
}

#' Recruitment model for transitions from spawning stock biomass (SSB) to 
#' juveniles age 1
#'
#' @param stock Fish Life-history parameters
#' @param life_history Life tables
#' @param spawner_recruit_model Initial spawner-recruit model
#'
#' @return Spawner-juvenile model
#'
#' @export
juvenile_recruitment <- function(
    stock,
    life_history,
    spawner_recruit_model) {
  
  juveniles <- life_history %>% filter(is.na(age_maturation))
  
  sjm <- spawner_juvenile_model(
    spawner_recruit_model, 
    stock %>% pull(events),
    juveniles)
  
  return (data.frame(
    spawner_recruit_model = spawner_recruit_model %>% pull(model),
    a = c(sjm$a),
    Rp = c(sjm$Rp)))
}


# Abundance vector matching class correspondence at demographic equilibrium and
# given spawning stock biomass (SSB)
#
# @param stock Fish Life-history parameters
# @param life_history Life tables
# @param fishing Table-defined fishing mortality
# @param ssb Initial SSB value 
# 
# @return Abundance per class ($N)
equilibrium_SSB <- function(
    stock,
    life_history,
    fishing,
    ssb) {
  
  scaler <- juveniles_adults_transitions(stock, life_history) %>%
    juveniles_adults_in_equilibrium()
  
  ssb_scaler <- scaler %>%
    segmented_SSB_pf_in_equilibrium(stock, life_history, fishing) %>%
    effective_SSB_equilibrium(stock, stock %>% pull(gsi))
  
  recruits <- ssb / ssb_scaler
  
  abundances <- data.frame(
    age = c(1),
    age_maturation = c(NA),
    abundance = c(recruits)) %>%
    rbind(
      scaler %>%
        mutate(abundance = coef * recruits) %>%
        select(-coef))
  
  return (list(N = abundances))
}

# Abundance vector at demographic equilibrium in case of
# density-dependent recruitment
#
# @param stock Fish Life-history parameters
# @param fishing Table-defined fishing mortality
# @param life_history Life tables
# 
# @return Abundance per class ($N)
equilibrium_density_dependent_recruitment <- function(
    stock, 
    fishing,
    life_history) {
  
  scaler <- juveniles_adults_transitions(stock, life_history) %>%
    juveniles_adults_in_equilibrium()
  
  recruits <- recruits_nonlinear(
    stock, 
    fishing, 
    stock %>% pull(gsi),
    life_history,
    scaler)
  
  abundances <- data.frame(
    age = c(1),
    age_maturation = c(NA),
    abundance = c(recruits)) %>%
    rbind(
      scaler %>%
        mutate(abundance = coef * recruits) %>%
        select(-coef))
  
  return (list(N = abundances))
}

# Abundance vector at demographic equilibrium in case of
# density-independent recruitment
#
# @param stock Fish Life-history parameters
# @param fishing Table-defined fishing mortality
# @param life_history Life tables
# 
# @return Abundance per class ($N) and growth rate ($fitness)
equilibrium_leslie_matrix <- function(
    stock,
    fishing,
    life_history) {
  
  ls <- leslie_matrix_reduced(
    stock,
    life_history,
    fishing,
    stock %>% pull(gsi))
  
  ev <- eigen(ls)
  
  if (max(Re(ev$values)) == 0) {
    return (list(
      N = life_history %>% 
        select(age, age_maturation) %>%
        cbind(data.frame(abundance = NA), 
      fitness = 0)))
  }
    
  
  pmax <- which(Re(ev$values) == max(Re(ev$values)))
  lmax <-  Re(ev$values[pmax])
  w <- abs(Re(ev$vectors[, pmax]))
  stable_distr <- w / w[1]
  
  abundances <- life_history %>% 
    cbind(data.frame(abundance = stable_distr)) %>%
    select(age, age_maturation, abundance)
  
  return (list(N = abundances, fitness = lmax))
}

# Abundance vector at the end of an annual cycle
#
# @param stock Fish Life-history parameters
# @param fishing Table-defined fishing mortality
# @param life_history_with_reproduction Life tables including reproduction
# @param state Abundance vector at the start of an annual cycle
#
# @return Abundance per class
transition <- function(
    stock,
    fishing,
    life_history_with_reproduction,
    state) {
  
  abundances_init <- state %>% pull(abundance)
  
  spawners_biomass <- sum(
    get_effective_SSB(life_history_with_reproduction) *
      abundances_init)
  
  recruits <- juveniles_age_1_production(
    spawners_biomass, 
    stock %>% select(spawner_recruit_model, a, Rp))
  
  transitions <- juveniles_adults_transitions(
    stock, 
    life_history_with_reproduction)
  
  temp_prev <- transitions %>%
    select(age, age_maturation) %>%
    cbind(abundance = head(abundances_init, -1))
  prev_juveniles <- temp_prev %>% filter(is.na(age_maturation)) 
  
  prev <- temp_prev %>%
    group_by(age_maturation) %>%
    group_modify(~ {
      
      if (is.na(.y %>% pull(age_maturation)))
        return (.x)
      
      prev_n_juveniles <- prev_juveniles %>% 
        filter(age == .y %>% pull(age_maturation)) %>%
        pull(abundance)
      
      if (nrow(.x) == 1) 
        prev_groupped_adults <- data.frame()
      else prev_groupped_adults <- .x %>% slice(2 : n())
      
      return (.x %>%
                slice(1) %>%
                mutate(abundance = prev_n_juveniles) %>%
                rbind(prev_groupped_adults))    
    }) %>%
    ungroup() %>%
    as.data.frame() %>%
    select(age, everything()) %>%
    arrange(!is.na(age_maturation), age_maturation, age)
  
  n_juveniles_adults <- (transitions %>% pull(coef)) * 
    (prev %>% pull(abundance))
  
  state_new <- data.frame(
    age = c(1),
    age_maturation = c(NA),
    abundance = c(recruits)) %>%
    rbind(
      transitions %>%
        select(age, age_maturation) %>%
        cbind(abundance = n_juveniles_adults))
  
  return (state_new)
}