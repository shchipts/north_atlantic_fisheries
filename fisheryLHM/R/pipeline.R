
# Helper function for life tables including reproduction process
get_life_history_with_reproduction <- function(stock, fishing) {
  
  life_history <- all_individuals(stock, fishing) %>%
    filter(is.na(age_maturation) | age >= age_maturation)
  life_history_with_reproduction <- life_history %>%
    left_join(
      reproduction(stock, life_history, fishing),
      by = c("age", "age_maturation"))
  
  return (life_history_with_reproduction)
} 

# Helper function for solution at demographic equilibrium
demographic_equilibrium <- function(stock, fishing) {
  
  get_eq <- ifelse(
    (stock %>% pull(spawner_recruit_model)) == "independent",
    equilibrium_leslie_matrix,
    equilibrium_density_dependent_recruitment)
  
  life_history <- get_life_history_with_reproduction(stock, fishing)
  
  eq <- get_eq(stock, fishing, life_history)
  
  return (append(
    eq,
    list(
      life_history = life_history)))
}

#' Properties of demographic equilibrium
#'
#' @description
#' The cases of density-independent reproduction (Leslie model) and density-
#' dependent reproduction (Beverton-Holt and Ricker models) are supported.
#' 
#' @param stock Fish life-history parameters
#' @param fishing Table-defined fishing mortality
#' 
#' @return Abundance, biomass, average age and body length per group 
#' 
#' @export
properties_at_equilibrium <- function(
    stock, 
    fishing) {
  
  output <- demographic_equilibrium(stock, fishing)
  abundances <- output$N %>% pull(abundance)
  
  data_summary <- properties(output$N, output$life_history) %>%
    cbind(properties_by_group(output$N, output$life_history)) %>%
    mutate(stock = stock %>% pull(stockname)) %>%
    select(stock, everything())
  
  return (data_summary)
}

#' Fish stock dynamics
#' 
#' @description
#' Population dynamics is determined by the spawner-recruit relationship and 
#' transition probabilities between classes corresponding to elements in the 
#' Leslie matrix. Changing life-history parameters alter spawner-recruit 
#' relationship and transition probabilities. Initial population distribution is
#' calibrated to the observed spawning stock biomass (SSB).
#'  
#' 
#' @param stock Initial life-history parameters
#' @param fishing Table-defined fishing mortality
#' @param ssb Initial SSB value
#' @param stock_dynamics Trajectory of changes in life-history parameters
#' 
#' @return Trajectory of changes in abundance, biomass, average age and 
#' body length (per group)
#' 
#' @export
dynamics <- function(
    stock,
    fishing,
    ssb,
    stock_dynamics) {
  
  life_history_init <- get_life_history_with_reproduction(stock, fishing)
  output <- equilibrium_SSB(
    stock,
    life_history_init,
    fishing,
    ssb)
  
  state_init <- output$N
  properties_init <- properties(state_init, life_history_init)
  
  data_summary <- Reduce(
    function(seed, idx) {
      
      stock_new <- stock_dynamics %>% slice(idx)
      life_history <- get_life_history_with_reproduction(stock_new, fishing)
      
      state_new <- transition(stock_new, fishing, life_history, seed$state)
      summary_new <- properties(state_new, life_history)
      
      return (list(
        state = state_new, 
        summary = rbind(seed$summary, summary_new)))
    },
    1 : nrow(stock_dynamics),
    init = list(
      state = state_init,
      summary = properties_init))$summary %>%
    mutate(stock = stock %>% pull(stockname)) %>%
    select(stock, everything())
  
  return (data_summary)
}