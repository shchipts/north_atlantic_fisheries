
# Helper function for survival probability
get_survival_prob <- function(
  natural_mortalities,
  fishing_mortalities) {
  
  sapply(
    1 : length(natural_mortalities),
    function(idx)(
      survival_probability(
        natural_mortalities[idx],
        fishing_mortalities[idx])))
}

# Helper function for catch probability
get_catch_prob <- function(
  natural_mortalities,
  fishing_mortalities) {
  
  sapply(
    1 : length(natural_mortalities),
    function(idx)(
      catch_probability(
        natural_mortalities[idx],
        fishing_mortalities[idx])))
}

# Helper function for natural death probability
get_natural_death_prob <- function(
  natural_mortalities,
  fishing_mortalities) {
  
  sapply(
    1 : length(natural_mortalities),
    function(idx)(
      naturally_dead_probabiity(
        natural_mortalities[idx],
        fishing_mortalities[idx])))
}

# Helper function for length-dependent fishing mortality
get_fishing_mortality <- function(body_lengths, stock, fishing) {

  if (is.null(fishing)) {
    
    F_func <- function(x)(return (0))
  }
  else {

    if (nrow(fishing) == 0)
      F_func <- function(x)(
        fishing_mortality(
          x,
          stock$F0,
          stock$F1,
          stock$F2,
          stock$F3,
          stock$F4,
          stock$F5,
          stock$F6))
    else {
      
      table <- fishing %>% arrange(length)
      values <- fishing %>% pull(fishing_mortality)
      
      F_func <- approxfun(
        fishing %>% pull(length),
        values,
        method = "linear",
        yleft = 0,
        yright = tail(values, 1))
    }
  }
  
  sapply(body_lengths, F_func)
}

# Helper function for length-dependent natural mortality
get_natural_mortality <- function(body_lengths, stock) {
  
  sapply(
    body_lengths,
    function(x)(
      natural_mortality(
        x,
        stock$gamma_3,
        stock$alpha_3,
        stock$length_reference)))
}

# Indicates whether one life-cycle process occurs after another or not
occurs_before <- function(process, events, relative_to) {
  
  shift(process, events, relative_to) == 0
}

growth_dependent_pars <- function(parents, stock, events) {
  
  gsi_base <- stock %>% pull(gsi)
  
  if (occurs_before("reproduction", events, "growth")) {
    data <- parents %>% 
      mutate(
        length_at_reproduction = body_length,
        weight_at_reproduction = body_weight,
        gsi_at_reproduction = effective_gsi,
        mature_at_reproduction = ifelse(
          is.na(age_maturation),
          0,
          1))
  }   
  else {
    length_new <- sapply(
      parents %>% pull(body_length),
      function(l)(
        length_adult(
          l, 
          stock %>% pull(gamma_1), 
          stock %>% pull(gamma_2), 
          stock %>% pull(alpha_1), 
          stock %>% pull(alpha_2), 
          gsi_base)))
    
    data <- parents %>% 
      mutate(length_at_reproduction = length_new) %>%
      mutate(
        weight_at_reproduction = weight(
          length_at_reproduction,
          stock %>% pull(gamma_2), 
          stock %>% pull(alpha_2),
          stock %>% pull(length_reference)),
        gsi_at_reproduction = gsi(
          length_at_reproduction,
          body_length,
          stock %>% pull(gamma_1), 
          stock %>% pull(gamma_2), 
          stock %>% pull(alpha_1), 
          stock %>% pull(alpha_2)),
        mature_at_reproduction = ifelse(
          is.na(age_maturation),
          maturation,
          1))
  }
  
  data <- data %>%
    mutate(
      effective_SSB_pf = ifelse(
        is.na(gsi_at_reproduction),
        0,
        ssb_pf_effective(
          ssb_pf(
            weight_at_reproduction,
            stock %>% pull(p_spf),
            stock %>% pull(F_spf)),
          gsi_at_reproduction,
          gsi_base)),
      real_spawners_proportion = ifelse(
        is.na(gsi_at_reproduction),
        0,
        spawner_proportion(stock %>% pull(p_spf), stock %>% pull(F_spf))))
  
  return (data)
}

mortality_dependent_pars <- function(parents, stock, events, fishing) {

  if (occurs_before("reproduction", events, "mortality")) {
    data <- parents %>%
      mutate(survivors_at_reproduction = 1)
  }
  else {
    if (occurs_before("mortality", events, "growth")) {
      data <- parents %>%
        mutate(survivors_at_reproduction = survival_length)
    }
    else {
      data <- parents %>%
        mutate(
          survivors_at_reproduction = survival_probability(
            natural_mortality(
              length_at_reproduction,
              stock %>% pull(gamma_3),
              stock %>% pull(alpha_3),
              stock %>% pull(length_reference)),
            get_fishing_mortality(
              c(length_at_reproduction),
              stock,
              fishing)[1]))
    }
  }
  
  return (data)
}

# Life tables for juveniles
get_juveniles <- function(stock, fishing) {
  
  ages <- seq(from = 1, to = stock$age_max)
  
  body_lengths <- tail(
    Reduce(
      function(s2, a) (
        append(
          s2, 
          length_juvenile(
            tail(s2, 1), 
            stock$gamma_1, 
            stock$gamma_2, 
            stock$alpha_1, 
            stock$alpha_2))),
      ages,
      init = c(0)),
    -1)
  
  maturation_probs <- sapply(
    seq(1: length(ages)),
    function(idx)(
      maturation_probability(
        ages[idx],
        body_lengths[idx], 
        maturation_steepness(
          stock$pmrn_width,
          stock$pmrn_envelope),
        stock$pmrn_slope,
        stock$pmrn_intercept)))
  
  
  body_weights <- sapply(
    body_lengths, 
    function(l)(
      weight(l, stock$gamma_2, stock$alpha_2, stock$length_reference)))
  
  fishing_mortalities <- get_fishing_mortality(body_lengths, stock, fishing)
  natural_mortalities <- get_natural_mortality(body_lengths, stock)
  
  catch_probs_spawning_ground <- sapply(
    1 : length(body_lengths),
    function(idx)(catch_probability_spawning_ground(stock$F_spf)))
  
  return (data.frame(
    age = ages,
    body_length = body_lengths,
    body_weight = body_weights,
    maturation = maturation_probs,
    fishing_mortality = fishing_mortalities,
    natural_mortality = natural_mortalities,
    survival_length = get_survival_prob(
      natural_mortalities, 
      fishing_mortalities),
    catch_length = get_catch_prob(
      natural_mortalities, 
      fishing_mortalities),
    natural_death_length = get_natural_death_prob(
      natural_mortalities, 
      fishing_mortalities),
    catch_spawning_ground = catch_probs_spawning_ground))
}

# Life tables for adults
get_adults <- function(stock, juveniles, fishing) {
  
  age_max <- stock$age_max
  ages_maturation <- seq(from = 2, to = age_max)
  
  processes <- Reduce(
    function(s, a_m) {
      
      ages <- seq(from = a_m, to = age_max)
      
      length_as_juvenile <- juveniles %>%
        filter(age == a_m - 1) %>%
        pull(body_length)
      
      body_lengths <- tail(
        Reduce(
          function(s2, a) (
            append(
              s2, 
              length_adult(
                tail(s2, 1), 
                stock$gamma_1, 
                stock$gamma_2, 
                stock$alpha_1, 
                stock$alpha_2,
                stock$gsi))),
          ages,
          init = length_as_juvenile),
        -1)
      
      indicators_RI <- tail(
        Reduce(
          function(s2, idx) {
            
            a <- ages[idx]
            l <- ifelse(a == a_m, length_as_juvenile, body_lengths[idx - 1])
            
            return (append(
              s2,
              indicator_reproductive_investment(
                l,
                stock$gamma_1, 
                stock$gamma_2, 
                stock$alpha_1, 
                stock$alpha_2,
                stock$gsi)))
          },
          1 : length(ages),
          init = c(TRUE)),
        -1)
      
      body_lengths_with_init <- append(length_as_juvenile, body_lengths)
      
      gsi <- sapply(
        2 : length(body_lengths_with_init),
        function(idx)(
          gsi(
            body_lengths_with_init[idx],
            body_lengths_with_init[idx - 1],
            stock$gamma_1, 
            stock$gamma_2, 
            stock$alpha_1, 
            stock$alpha_2)))
      
      body_weights <- sapply(
        body_lengths, 
        function(l)(
          weight(l, stock$gamma_2, stock$alpha_2, stock$length_reference)))
      
      fishing_mortalities <- get_fishing_mortality(body_lengths, stock, fishing)
      natural_mortalities <- get_natural_mortality(body_lengths, stock)
      
      survival_probs_spawning_ground <- sapply(
        1 : length(body_lengths),
        function(idx)(
          survival_probability_spawning_ground(
            stock$F_spf,
            stock$p_spm)))
      
      catch_probs_spawning_ground <- sapply(
        1 : length(body_lengths),
        function(idx)(catch_probability_spawning_ground(stock$F_spf)))
      
      d <- data.frame(
        age = ages,
        age_maturation = rep(a_m, length(ages)),
        body_length = body_lengths,
        body_weight = body_weights,
        indicator_RI = indicators_RI,
        fishing_mortality = fishing_mortalities,
        natural_mortality = natural_mortalities,
        survival_length = get_survival_prob(
          natural_mortalities, 
          fishing_mortalities),
        catch_length = get_catch_prob(
          natural_mortalities, 
          fishing_mortalities),
        natural_death_length = get_natural_death_prob(
          natural_mortalities, 
          fishing_mortalities),
        survival_spawning_ground = survival_probs_spawning_ground,
        catch_spawning_ground = catch_probs_spawning_ground,
        effective_gsi = gsi)
      
      return (rbind(s, d))
    },
    ages_maturation,
    init = data.frame())
}

#' Life tables per age and maturation age 
#'
#' @param stock Fish life-history parameters
#' @param fishing Table-defined fishing mortality
#' 
#' @export
all_individuals <- function(stock, fishing) {
  
  age_max <- stock$age_max
  ages <- seq(from = 1, to = age_max)
  ages_maturation <- seq(from = 2, to = age_max)
  
  juveniles <- get_juveniles(stock, fishing)
  adults <- get_adults(stock, juveniles, fishing)
  
  fill_missing_juvenile_vals <- function(data)(
    data %>%
      mutate(survival_spawning_ground = NA) %>%
      mutate(indicator_RI = FALSE) %>%
      mutate(effective_gsi = NA))
  
  only_juvenile <- juveniles %>%
    mutate(age_maturation = NA) %>%
    fill_missing_juvenile_vals()
  
  data_all <- rbind(
    only_juvenile,
    Reduce(
      function(s, a_m) {
        
        d <- juveniles %>%
          filter(age < a_m) %>%
          mutate(age_maturation = a_m) %>%
          fill_missing_juvenile_vals() %>%
          rbind(
            adults %>%
              filter(age_maturation == a_m) %>%
              mutate(maturation = NA))
        
        return (rbind(s, d))
        },
      ages_maturation,
      init = data.frame())) %>%
    select(
      age, 
      age_maturation,
      indicator_RI,
      body_length,
      maturation, 
      body_weight,
      natural_mortality,
      fishing_mortality,
      survival_length, 
      catch_length,
      natural_death_length,
      survival_spawning_ground,
      catch_spawning_ground,
      effective_gsi)
  
  return (data_all)
}

# Reproduction-related life table per age and maturation age 
#
# @param stock Fish life-history parameters
# @param life_history Life tables for other processes
# @param fishing Table-defined fishing mortality
# 
reproduction <- function(stock, life_history, fishing) {
  
  events <- stock %>% pull(events)
  gsi_base <- stock %>% pull(gsi)
  
  data <- life_history %>%
    select(-survival_spawning_ground) %>%
    filter(is.na(age_maturation) | age >= age_maturation) %>%
    growth_dependent_pars(stock, events) %>%
    mortality_dependent_pars(stock, events, fishing) %>%
    select(
      age, 
      age_maturation,
      real_spawners_proportion,
      effective_SSB_pf,
      mature_at_reproduction,
      survivors_at_reproduction)
  
  return (data)
}

# Returns number of spawners
#
# @param life_history Life tables
# 
get_spawners <- function(life_history) {
  
  (life_history %>% pull(real_spawners_proportion)) *
    (life_history %>% pull(mature_at_reproduction)) *
    (life_history %>% pull(survivors_at_reproduction))
}

# Returns effective spawning biomass per mature fish 
#
# @param life_history Life tables
#
get_effective_SSB <- function(life_history) {
  
  (life_history %>% pull(effective_SSB_pf)) *
    (life_history %>% pull(mature_at_reproduction)) *
    (life_history %>% pull(survivors_at_reproduction))
}