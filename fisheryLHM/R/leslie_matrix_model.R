
# Data frame with a full-size Leslie matrix
leslie_matrix_full <- function(
  stock,
  life_history,
  fishing,
  gsi_base) {
  
  age_max <- stock$age_max
  events <- stock$events
  
  shifts <- data.frame(
    mortality = c(shift("mortality", events)),
    reproduction = c(shift("reproduction", events)))
  
  a_sjm <- stock %>% pull(a)
  reproduction_table <- reproduction(stock, life_history, fishing) %>%
    mutate(
      scaling = a_sjm *
        ssb_pf_effective(effective_SSB_pf, stock %>% pull(gsi), gsi_base) *
        mature_at_reproduction *
        survivors_at_reproduction)
  
  juveniles <- life_history %>% filter(is.na(age_maturation))
  new_adults <- life_history %>%
    filter(age == age_maturation)
  
  ### juveniles_to_juveniles
  juvenile_transitions <- sapply(
    1 : (age_max - 1),   # ages by column
    function(t){
      m <- juveniles %>% 
        slice(t) %>% 
        pull(maturation)
      s <- juveniles %>% 
        slice(t + shifts %>% pull(mortality)) %>%
        pull(survival_length)
      
      return ((1 - m) * s)
    })
  
  juveniles_to_juveniles <- matrix(
    reproduction_table %>% 
      filter(is.na(age_maturation)) %>% 
      pull(scaling), 
    nrow = 1, 
    ncol = age_max) %>%
    rbind(diag(juvenile_transitions) %>% 
            cbind(matrix(0, nrow = age_max - 1, ncol = 1)))

  ### adults_to_juveniles
  all_to_juveniles <- Reduce(
    function(seed, a_m) {
      
      reproduction_selected <- reproduction_table %>%
        filter(age_maturation == a_m)
      
      cbind(
        seed,
        sapply(
          1 : age_max,   # ages by column
          function(a)(
            if (a < a_m)
              0
            else reproduction_selected %>% 
              slice(a - a_m + 1) %>% 
              pull(scaling))) %>% 
          matrix(nrow = 1) %>%
          rbind(matrix(0, nrow = age_max - 1, ncol = age_max)))
    },
    seq(from = 2, to = age_max),   # age_maturation blocks stacked by column 
    init = juveniles_to_juveniles)
  
  all_to_all_adults <- Reduce(
    function(seed, a_m) {
      
      ### juveniles_to_adults
      juveniles_to_adults <- matrix(0, nrow = age_max, ncol = age_max)
      
      t <- a_m - 1   # juvenile age
      
      m <- juveniles %>% slice(t) %>% pull(maturation)
      s <- ifelse(
        shifts %>% pull(mortality) > 0,
        new_adults %>% slice(t) %>% pull(survival_length),
        juveniles %>% slice(t) %>% pull(survival_length))
      
      if (shifts %>% pull(reproduction) == 1) {
        s_spawning_ground <- new_adults %>%
          slice(t) %>%
          pull(survival_spawning_ground)
        
        s <- s * s_spawning_ground
      }
      
      juveniles_to_adults[t + 1, t] <- m * s
      
      all_to_adults_same_maturation_age <- Reduce(
        function(s2, a_block) {
          
          if (a_block == a_m) {
            
            ### adults_to_adults_same_maturation_age
            
            adults_selected <- life_history %>% 
              filter(age_maturation == a_m & age >= age_maturation)
            
            size <- age_max - a_m
            
            adults_to_adults <- matrix(0, nrow = a_m, ncol = age_max)
            
            if (size != 0) {
              
              adult_survivals <- sapply(
                adults_selected %>% filter(age > age_maturation) %>% pull(age),
                function(a) {
                  
                  idx <- a - a_m
                  
                  s <- adults_selected %>% 
                    slice(idx + shifts %>% pull(mortality)) %>%
                    pull(survival_length)
                  
                  s_spawning_ground <- adults_selected %>%
                    slice(idx + shifts %>% pull(reproduction)) %>%
                    pull(survival_spawning_ground)
                  
                  return (s * s_spawning_ground)
                })
              
              adults_to_adults <- rbind(
                adults_to_adults,
                matrix(0, nrow = size, ncol = a_m - 1) %>%
                  cbind(diag(adult_survivals, nrow = size)) %>%
                  cbind(matrix(0, nrow = size, ncol = 1)))
            }
          }
          else (adults_to_adults <- matrix(0, nrow = age_max, ncol = age_max))
          
          return (cbind(s2, adults_to_adults))  
        },
        seq(from = 2, to = age_max),   # age_maturation blocks stacked by column 
        init = juveniles_to_adults)
      
      return (rbind(seed, all_to_adults_same_maturation_age))
    },
    seq(from = 2, to = age_max),   # age_maturation blocks stacked by row 
    init = matrix(0, nrow = 0, ncol = age_max * age_max))
  
  return (rbind(all_to_juveniles, all_to_all_adults))
}

# Transition probabilities from juveniles to juveniles
juveniles_to_juveniles_transitions <- function(life_history, shifts) {
  
  shift_mortality <- shifts %>% pull(mortality)
  
  juveniles <- life_history %>% 
    filter(is.na(age_maturation))
  
  maturation <- juveniles %>% 
    slice(1 : (nrow(juveniles) - 1)) %>%
    pull(maturation)
  
  survival <- juveniles %>%
    slice((1 + shift_mortality) : (nrow(juveniles) - 1 + shift_mortality)) %>%
    pull(survival_length)
  
  return (juveniles %>%
            filter(age > 1) %>%
            select(age, age_maturation) %>%
            cbind(coef = (1 - maturation) * survival))
}

# Transition probabilities from juveniles to adults
juveniles_to_adults_transitions <- function(data, life_history, shifts) {
  
  new_adults <- life_history %>% filter(age == age_maturation)
  juveniles <- life_history %>% filter(is.na(age_maturation))
  
  ages <- seq(
    from = (new_adults %>% slice(1) %>% pull(age_maturation)),
    to = (new_adults %>% slice(n()) %>% pull(age_maturation))) - 1
  
  maturation <- juveniles %>%
    slice(ages) %>%
    pull(maturation)
  
  if (shifts %>% pull(mortality) > 0)
    survival <- new_adults %>% slice(ages) %>% pull(survival_length)
  else
    survival <- juveniles %>% slice(ages) %>% pull(survival_length)
  
  if (shifts %>% pull(reproduction) == 1) {
    
    survival_spawning_ground <- new_adults %>%
      slice(ages) %>%
      pull(survival_spawning_ground)
    
    survival <- survival * survival_spawning_ground
  }
  
  return (data %>%
            rbind(
              new_adults %>%
                slice(ages) %>%
                select(age, age_maturation) %>%
                cbind(coef = maturation * survival)))
}

# Transition probabilities for adults in the same maturity class
adults_to_adults_same_maturation_age <- function(
    a_m,
    life_history, 
    shifts) {
  
  adults_selected <- life_history %>% 
    filter(age_maturation == a_m & age >= age_maturation)
  
  adults_to_adults_selected <- adults_selected %>% filter(age > age_maturation)
  
  ages <- adults_to_adults_selected %>% pull(age)
  
  survival <- adults_selected %>% 
    slice(ages - a_m + shifts %>% pull(mortality)) %>%
    pull(survival_length)
  
  survival_spawning_ground <- adults_selected %>%
    slice(ages - a_m + shifts %>% pull(reproduction)) %>%
    pull(survival_spawning_ground)
  
  return (data.frame(
    age = ages,
    coef = survival * survival_spawning_ground) %>%
      mutate(age_maturation = a_m))
}

# Transition probabilities from adults to adults
adults_to_adults_transitions <- function(data, life_history, shifts) {
  
  new_adults <- data %>% filter(!is.na(age_maturation))
  
  data_new <- Reduce(
    function(seed, a_m) {
      
      t <- a_m - 1
      
      d <- adults_to_adults_same_maturation_age(
        a_m,
        life_history,
        shifts)
      
      return (rbind(seed, d))
    },
    new_adults %>% pull(age_maturation),
    init = data)
  
  return (data_new)
}

# Spawner-juvenile model
#
# @param spawner_recruit_model Spawner-recruit model
# @param events Order of life-cycle processes
# @param juveniles Life tables for juveniles
#
# @return Slope and peak recruitment
spawner_juvenile_model <- function(
    spawner_recruit_model, 
    events, 
    juveniles) {
  
  a_r <- spawner_recruit_model %>% pull(a_recruitment)
  
  if (a_r == 1)
    c <- 1 
  else {
    inc <- shift("mortality", events) 
    
    c <- Reduce(
      function(seed, idx){
        
        m <- juveniles %>% slice(idx) %>% pull(maturation)
        s <- juveniles %>% slice(idx + inc) %>% pull(survival_length)
        
        return (seed * s * (1 - m))
      },
      1 : a_r,
      init = 1)
  }
  
  a_sjm <- spawner_recruit_model %>% pull(a) / c
  
  Rp_srm <- spawner_recruit_model %>% pull(Rp)
  
  model <- spawner_recruit_model %>% pull(model)
  
  Rp_sjm <- case_when(
    model == "independent" ~ NA_real_,
    model == "beverton_holt" ~ Rp_srm / c,
    model == "ricker" ~ Rp_srm / c)
  
  return (list(a = a_sjm, Rp = Rp_sjm))
}

# Leslie matrix
#
# @param stock Fish Life-history parameters
# @param life_history Life tables
# @param fishing Table-defined fishing mortality
# @param gsi_base Initial gonadosomatic index
#
# @return A data frame with elements of a Leslie matrix 
leslie_matrix_reduced <- function(
    stock,
    life_history,
    fishing,
    gsi_base) {
  
  age_max <- stock$age_max
  
  return (Reduce(
    function(s, a_m) {
      start <- (age_max + (age_max - a_m + 2)) * (a_m - 1) / 2 + 1
      end <- start + (a_m - 2)
      
      indexes <- seq(from = start, to = end)
      
      return (s[- indexes, - indexes])
    },
    2 : age_max,
    init = leslie_matrix_full(stock, life_history, fishing, gsi_base)))
}

# Transition probabilities to next time step
#
# @param stock Fish Life-history parameters
# @param life_history Life tables
#
# @return A data with transition probabilities for a population abundance vector
juveniles_adults_transitions <- function(
    stock, 
    life_history) {
  
  events <- stock %>% pull(events)
  
  shifts <- data.frame(
    mortality = c(shift("mortality", events)),
    reproduction = c(shift("reproduction", events)))
  
  data <- juveniles_to_juveniles_transitions(life_history, shifts) %>%
    juveniles_to_adults_transitions(life_history, shifts) %>%
    adults_to_adults_transitions(life_history, shifts) %>%
    arrange(!is.na(age_maturation), age_maturation, age)
  
  return (data)
}

# Effective spawning stock biomass (SSB) per class
#
# @param equilibrium_rest Juveniles and adults
# @param stock Fish Life-history parameters
# @param life_history Life tables
# @param fishing Table-defined fishing mortality
#
# @return Scaling factors for recruitment level to obtain effective SSB per class
segmented_SSB_pf_in_equilibrium <- function(
    equilibrium_rest,
    stock, 
    life_history,
    fishing) {
  
  data <- reproduction(stock, life_history, fishing) %>%
    left_join(
      equilibrium_rest %>%
        add_row(age = 1, age_maturation = NA, coef = 1),
      by = c("age", "age_maturation")) %>%
    mutate(
      coef_SSB = effective_SSB_pf *
        mature_at_reproduction *
        survivors_at_reproduction *
        coef) %>%
    select(age, age_maturation, coef_SSB)
  
  return (data)
}

# Effective spawning stock biomass (SSB)
#
# @param SSB_pf_segmented SSB per class
# @param stock Fish Life-history parameters
# @param gsi_base Initial gonadosomatic index
#
# @return Scaling factor for recruitment level to obtain effective SSB
effective_SSB_equilibrium <- function(
    SSB_pf_segmented,
    stock,
    gsi_base) {
  
  SSB_pf <- SSB_pf_segmented %>%
    pull(coef_SSB) %>%
    sum()
  
  SSB_pf_effective <- ssb_pf_effective(
    SSB_pf,
    stock %>% pull(gsi),
    gsi_base)
  
  return (SSB_pf_effective)
}