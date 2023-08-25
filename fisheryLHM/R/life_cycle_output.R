
get_weight <- function(life_history) {
  
  life_history %>% pull(body_weight)
}

get_pmrn <- function(life_history) {
  
  life_history %>% pull(maturation)
}

get_nat_dead <- function(life_history) {
  
  life_history %>% pull(natural_death_length)
}

get_catch <- function(life_history) {
  
  life_history %>% pull(catch_length)
}

get_weighted_var <- function(
    abundances, 
    life_history,
    property, 
    get_fraction = function(p)(1)) {
  
  ns <- get_fraction(life_history) * abundances
  
  return ((life_history %>% pull(!!property)) * ns / 
            sum(ns, na.rm = TRUE))
}

# Number of individuals harvested at the spawning grounds
get_catch_spawning_ground <- function(life_history) {
  
  (life_history %>% pull(catch_spawning_ground)) *
    (life_history %>% pull(mature_at_reproduction)) *
    (life_history %>% pull(survivors_at_reproduction))
}

# Biomass of fish harvested at the spawning grounds
get_catch_spawn_biomass <- function(
    abundances,
    weights,
    life_history) {
  
  return (sum(
    weights * 
      get_catch_spawning_ground(life_history) * 
      abundances))
}

# Properties of demographic equilibrium
#
# @param equilibrium Demographic equilibrium
# @param life_history_with_reproduction Life tables
#
# @return Abundance, biomass, average age and body length per group
properties <- function(
    equilibrium,
    life_history_with_reproduction) {
  
  abundances <- equilibrium %>% pull(abundance)
  weights <- get_weight(life_history_with_reproduction)
  
  get_total_catch <- function(p)(get_catch(p) + get_catch_spawning_ground(p))
  
  data_properties <- data.frame(
    abundance = sum(abundances),
    recruits = abundances[1],
    stock_biomass = sum(weights * abundances),
    age =  sum(
      get_weighted_var(
        abundances, 
        life_history_with_reproduction, 
        "age")),
    length = sum(
      get_weighted_var(
        abundances, 
        life_history_with_reproduction,
        "body_length")),
    maturing_fish_abundance = sum(
      get_pmrn(life_history_with_reproduction) * abundances, 
      na.rm = TRUE),
    maturing_fish_biomass = sum(
      weights * get_pmrn(life_history_with_reproduction) * abundances, 
      na.rm = TRUE),
    maturing_fish_age = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction,
        "age", 
        get_pmrn), 
      na.rm = TRUE),
    maturing_fish_length = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction,
        "body_length", 
        get_pmrn), 
      na.rm = TRUE),
    nat_dead_fish_abundance = sum(
      get_nat_dead(life_history_with_reproduction) * 
        abundances),
    nat_dead_fish_biomass = sum(
      weights * 
        get_nat_dead(life_history_with_reproduction) * 
        abundances),
    nat_dead_age = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction, 
        "age", 
        get_nat_dead),
      na.rm = TRUE),
    nat_dead_length = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction, 
        "body_length", 
        get_nat_dead),
      na.rm = TRUE),
    catch_n = sum(
      get_catch(life_history_with_reproduction) * 
        abundances),
    catch_biomass = sum(
      weights * 
        get_catch(life_history_with_reproduction) *
        abundances),
    catch_spawn_n = sum(
      get_catch_spawning_ground(life_history_with_reproduction) * 
        abundances),
    catch_spawn_biomass = get_catch_spawn_biomass(
      abundances,
      weights,
      life_history_with_reproduction),
    catch_age = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction, 
        "age", 
        get_total_catch)),
    catch_length = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction,
        "body_length", 
        get_total_catch)),
    spawners_abundance = sum(
      get_spawners(life_history_with_reproduction) *
        abundances),
    spawners_biomass = sum(
      get_effective_SSB(life_history_with_reproduction) * 
        abundances),
    spawners_age = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction, 
        "age", 
        get_spawners)),
    spawners_length = sum(
      get_weighted_var(
        abundances,
        life_history_with_reproduction, 
        "body_length", 
        get_spawners))) %>%
    mutate(
      catch = catch_n + catch_spawn_n,
      yield = catch_biomass + catch_spawn_biomass) %>%
    select(-c(
      catch_n, 
      catch_spawn_n, 
      catch_biomass, 
      catch_spawn_biomass))
  
  return (data_properties)
}

# Properties of demographic equilibrium per juveniles and adults
#
# @param equilibrium Demographic equilibrium
# @param life_history_with_reproduction Life tables
#
# @return Total biomass and yield for juveniles and adults
properties_by_group <- function(
    equilibrium,
    life_history_with_reproduction) {
  
  data_cat <- life_history_with_reproduction %>%
    mutate(juvenile = as.numeric(is.na(age_maturation))) %>%
    mutate(adult = as.numeric(!juvenile))
  juveniles <- data_cat %>% pull(juvenile)
  adults <- data_cat %>% pull(adult)
  
  abundances <- equilibrium %>% pull(abundance)
  weights <- get_weight(life_history_with_reproduction)
  
  data_properties <- data.frame(
    stock_biomass_juveniles = sum(weights * abundances * juveniles),
    stock_biomass_adults = sum(weights * abundances * adults),
    catch_biomass_adults = sum(
      weights * 
        get_catch(life_history_with_reproduction) * 
        abundances * 
        adults),
    catch_biomass_juveniles = sum(
      weights * 
        get_catch(life_history_with_reproduction) * 
        abundances * 
        juveniles),
    catch_spawn_biomass = get_catch_spawn_biomass(
      abundances,
      weights,
      life_history_with_reproduction)) %>%
    mutate(
      yield_adults = catch_biomass_adults + catch_spawn_biomass,
      yield_juveniles = catch_biomass_juveniles) %>%
    select(-c(
      catch_biomass_adults,
      catch_biomass_juveniles,
      catch_spawn_biomass))
  
  return (data_properties)
}