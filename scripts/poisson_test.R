## Code to ensure that the used packages are installed 
list_of_packages_used <- 
  c(
    "broom",
    "DHARMa", 
    "AER",
    "dplyr"
  )

new_packages                 <- list_of_packages_used[!(list_of_packages_used %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

## Load libraries ----
lapply(list_of_packages_used, library, character.only = TRUE)

# Keep environment clean
rm(list_of_packages_used, new_packages)

# Define the function
check_poisson_assumptions <- function(model) {
  # Create a DHARMa object from the model
  simulationOutput <- simulateResiduals(fittedModel = model)
  
  # Perform tests
  dispersion_test <- testDispersion(simulationOutput, plot = F)
  zero_inflation_test <- testZeroInflation(simulationOutput, plot = F)
  
  # Test for normality of residuals
  normality_test <- shapiro.test(resid(model))
  
  # Test for homoscedasticity (constant variance) of residuals
  
  # Return a data frame with the p-values
  results <- data.frame(
    Test = c("Dispersion", "Zero Inflation", "Normality of Residuals"),
    p_value = c(dispersion_test$p.value, zero_inflation_test$p.value, normality_test$p.value),
    Decision = ifelse(c(dispersion_test$p.value, zero_inflation_test$p.value, normality_test$p.value) < 0.05, "Reject Null", "Fail to Reject Null")
  )
  results <-
    results |>
    dplyr::mutate(Decision = dplyr::case_when(
      Decision == "Reject Null" & Test == "Dispersion" ~ "Overdispersion detected",
      Decision == "Fail to Reject Null" & Test == "Dispersion" ~ "No overdispersion detected",
      Decision == "Reject Null" & Test == "Zero Inflation" ~ "Zero-inflation detected",
      Decision == "Fail to Reject Null" & Test == "Zero Inflation" ~ "No Zero-inflation detected",
      Decision == "Reject Null" & Test == "Normality of Residuals" ~ "Residuals non-normal",
      Decision == "Fail to Reject Null" & Test == "Normality of Residuals" ~ "Residuals normal"
    ))
  
  return(results)
}

# Usage:
# Assume 'model' is your Poisson GLM
# model <- glm(y ~ x, data = df, family = poisson)
# check_poisson_assumptions(model)