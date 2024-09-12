# Function to simulate a simple branching process where the simulation stops 
# when total born individuals exceed 1000 or no one is alive

simulate_branching_process <- function(beta, gamma,delta) {

  # Initialize variables
  population <- 1
  pop_sizes <- c(1)  # Record population sizes
  total_born <- 1  # Start with initial population as already born
  
  while (population > 0 && total_born <= 1000) {
    # Total rate (birth + death)
    total_rate <- population * (beta +gamma+delta)
    
    # Time until the next event
    delta_t <- rexp(1, rate = total_rate)
    
    # Determine if it's a birth or a death event
    if (runif(1) < beta / (beta + gamma+ delta)) {
      # Birth event
      population <- population + 1
      total_born <- total_born + 1
    } else {
      # Death event
      population <- population - 1
    }
    
    pop_sizes <- c(pop_sizes, population)
  }
  
  # Return results as a data frame
  return(total_born)
}

set.seed(123)

beta <- 0.8   # Birth rate
gamma <- 1/7  # Death rate
delta <- 1/7  # Death rate

count_1000_or_more <- 0
# Run the simulation 10,000 times
for (i in 1:10000) {
  total_born <- simulate_branching_process(beta, gamma,delta)
  if (total_born >= 1000) {
    count_1000_or_more <- count_1000_or_more + 1
  }
}
print(paste("fraction of having >=1000 individuals born:", count_1000_or_more/10000))
