set.seed(123)  # For reproducibility

# Parameters
beta <- 0.8   # Infection rate
gamma <- 1/7  # Recovery rate
delta <- 1/7  # Diagnosis rate
pi <- 2/3   # Probability that an individual is an app-user
p <- 2/3    # Probability of manual tracing link

# Function to initialize the first component with a type
initialize_first_component <- function(pi) {
  type <- ifelse(runif(1) < pi, 1, 2)  # Decide the type based on probability pi
  if (type == 1) {
    return(list(type = 1, k = 1, l = 0))  # Type-1 starts with one app-user
  } else {
    return(list(type = 2, k = 0, l = 1))  # Type-2 starts with one non-app-user
  }
}

# Simulate the evolution of a component
simulate_component <- function(component, total_individuals_born) {
  new_components <- list()  # Track newly generated components
  
  while (component$k > 0 || component$l > 0) {
    k <- component$k
    l <- component$l
    total_individuals <- k + l
    
    # Rates of events
    rate_app_infection <- (k * beta * pi + l * beta * pi * p)
    rate_app_recovery <- k * gamma
    rate_non_app_infection <- total_individuals * beta * (1 - pi) * p
    rate_non_app_recovery <- l * gamma
    rate_diagnosis <- total_individuals * delta
    rate_new_component_app <- l * beta * pi * (1 - p)
    rate_new_component_non_app <- (k + l) * beta * (1 - pi) * (1 - p)
    
    # Total rate
    total_rate <- rate_app_infection + rate_app_recovery + rate_non_app_infection + rate_non_app_recovery +
      rate_diagnosis + rate_new_component_app + rate_new_component_non_app
    event_prob <- runif(1) * total_rate
    
    # Execute the event
    if (event_prob < rate_app_infection) {
      component$k <- component$k + 1
      total_individuals_born <- total_individuals_born + 1
    } else if (event_prob < rate_app_infection + rate_app_recovery) {
      component$k <- max(component$k - 1, 0)
    } else if (event_prob < rate_app_infection + rate_app_recovery + rate_non_app_infection) {
      component$l <- component$l + 1
      total_individuals_born <- total_individuals_born + 1
    } else if (event_prob < rate_app_infection + rate_app_recovery + rate_non_app_infection + rate_non_app_recovery) {
      component$l <- max(component$l - 1, 0)
    } else if (event_prob < rate_app_infection + rate_app_recovery + rate_non_app_infection + rate_non_app_recovery + rate_diagnosis) {
      component$k <- 0
      component$l <- 0
    } else if (event_prob < rate_app_infection + rate_app_recovery + rate_non_app_infection + rate_non_app_recovery + rate_diagnosis + rate_new_component_app) {
      # Generate a new component with an app-user root
      new_components <- append(new_components, list(list(type = 1, k = 1, l = 0)))
      total_individuals_born <- total_individuals_born + 1
    } else {
      # Generate a new component with a non-app-user root
      new_components <- append(new_components, list(list(type = 2, k = 0, l = 1)))
      total_individuals_born <- total_individuals_born + 1
    }
  }
  
  return(list(new_components = new_components, total_born = total_individuals_born))
}

# Recursive function to handle the creation of new components
run_simulation <- function() {
  
  total_individuals_born <- 0
  
  # Initialize the first component
  component <- initialize_first_component(pi)
  components <- list(component)  # Start with one component
  
  while (total_individuals_born < 1000 && length(components) > 0) {
    current_component <- components[[1]]
    components <- components[-1]  # Remove the first component from the list
    
    # Simulate the current component
    result <- simulate_component(current_component, total_individuals_born)
    
    # Update the total number of individuals born
    total_individuals_born <- result$total_born
    
    # Add new components to the list
    components <- append(components, result$new_components)
  }
  
  return(total_individuals_born)
}




# Run the simulation 10 000 times

count_1000_or_more <- 0

for (i in 1:10000) {
  print(i)
  total_individuals_born <- run_simulation()
  if (total_individuals_born >= 1000) {
    count_1000_or_more <- count_1000_or_more + 1
  }
}
print(paste("fraction of having >=1000 individuals born:", count_1000_or_more/10000))




