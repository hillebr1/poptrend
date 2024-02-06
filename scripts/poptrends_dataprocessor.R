# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                     poptrends_dataprocessor ----                                   
#                                     15/01/2024 
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# This code is an adaptation from the original poptrends.R 
# Logistic regressions were added to the code (poisson and binomial) following the discussions during the workshop on Thursday 11-01-2024 

# Load and install packages ----

## Code to ensure that the used packages are installed 
list_of_packages_used <- 
  c(
    "tidyverse", # data manipulation
    "vegan", # Analysis species data
    "lubridate", # Time-date manipulations
    "broom", # tidy output models
    "gvlma", # quick test for lm assumptions
    "MASS", # negative binomial modes
    "purrr"
    )

new_packages                 <- list_of_packages_used[!(list_of_packages_used %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

## Load libraries ----
lapply(list_of_packages_used, library, character.only = TRUE)

# Keep environment clean
rm(list_of_packages_used, new_packages)

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read in data and prepare for analysis ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# As an example we continue to use the phytoplankton dataset that's included in this github project. Replace with your own data following the same format:
# StationID, year, species, abu
# Where abu are the annual means per species per StationID

# Change input file here. If neccessary perform similar clean up steps
data <- read_delim("data/PPKTcount_noHT_28092022.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)
names(data)

# STEP 1: If your variable names do not conform to our standard (StationID, year, species, abu), rename and delete empty rows 
data<-data |> dplyr::rename(species = Species, abu = biovolume_l) 
# in our case StationID is okay, year needs to be calculated
data<-data[!is.na(data$abu),]
data<-data[!is.na(data$sampledate),]

# STEP 2: If sample date but not year is given, extract the year. Library lubridate allows to first change the dates to numeric values and then extract the year
data$date.numeric<-lubridate::dmy(data$sampledate) 
# dmy is command if data are in format dd.mm.yyyy, ymd is for yyyy.mm.dd
data$year <- lubridate::year(data$date.numeric) 
# extracts the year

# STEP 3: If your data can contain the same species more than once per sample (because different size classes or lifestages were distinguished) sum them
data<-data  |>
  dplyr::group_by(StationID, sampledate, year, species) |>  
  dplyr::summarise(abu= sum(abu, na.rm = T))
# obviously sampledate and STEP 4 can be ignored if only annual data exist

# STEP 4: Create annual means and the number of annual samples 
data<-data  |>
  dplyr::group_by(StationID, year, species) |> 
  dplyr::summarise(abu= mean(abu, na.rm = T),
            N.sample=length(abu))

# STEP 5: Further reductions
# if you need to exclude some StationID because they are not part of the WaddenSea or for other reasons
data<-data |>
  dplyr::filter(StationID %in% c("BOOMKDP", "DANTZGT", "DOOVBWT","GROOTGND","BOCHTVWTM",
                          "HUIBGOT", "MARSDND", "ROTTMPT3", "TERSLG10",
                          "Nney_W_2", "Bork_W_1", "WeMu_W_1", "JaBu_W_1"))

# In other data sets it might be needed to exclude certain species or years.

# Our example is already in long format, if your data is in wide format consider the "melt" command as in mentioned here, remove inital # to run
#data <- melt(data, 
#                 id.vars=c("StationID", "year"), # all the variables to keep but not split apart
#                 measure.vars=c(x:y), # number of the variables containing species, that shall be melted into one
#                 variable.name="species", # name for new variable that contains the names of the melded variables
#                 value.name="abu") # name of the new variables with the values


length.spec <-data |> 
  group_by(StationID, species) |> 
  summarise(N.spec=length(unique(year)))

length.stat <-data |> 
  group_by(StationID) |> 
  summarise(N.stat=length(unique(year)))
length<-merge(length.spec,length.stat,by="StationID")
length$threshold<-length$N.spec/length$N.stat
data<-merge(data,length,by=c("StationID","species"))

data<-data[data$N.spec>4,]
data<-data[data$threshold>.75,]


# Start analysis script from here ----

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Create final dataframe for running through the analysis ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Make sure to apply the filtering step to only include species and station combinations where a species has been observed for at least 5 years and a threshold number of years of at least 75%

filter.threshold <- 
  data |> 
  dplyr::group_by(StationID, species) |>
  # Remove 0's
  dplyr::filter(abu != 0) |>
  dplyr::summarise(N.spec = length(unique(year))) |>
  dplyr::left_join(data |> 
                     dplyr::group_by(StationID) |>
                     dplyr::summarise(N.stat = length(unique(year))),
                   by = "StationID") |>
  dplyr::mutate(threshold = N.spec/N.stat) |>
  dplyr::filter(N.spec >= 5 & threshold > 0.75) |>
  dplyr::mutate(UCI = paste(StationID, species, sep = " "))

# filter data and add 0's
data <-
  data |>
  dplyr::mutate(UCI = paste(StationID, species, sep = " ")) |>
  dplyr::filter(UCI %in% filter.threshold$UCI) |>
  dplyr::select(StationID, year, species, abu) |>
  tidyr::pivot_wider(id_cols = c(StationID, year), values_from = abu, names_from = species, values_fill = 0) |>
  tidyr::pivot_longer(cols = !c(StationID, year), names_to = "species", values_to = "abu") |>
  # Remove species that only contain 0 for certain stations
  dplyr::group_by(StationID, species) |>
  dplyr::filter(any(abu != 0)) |>
  dplyr::ungroup()

# Finalize the data in one format that can run through the entire script without having to change transformations everywhere
# Make sure the abundance data are integers (for poisson models_) and the year is continuous
# Also make a transformed version of abu for a linear regression. Change transformation in this step if desired, or remove transformation to perform linear regression on original data
# create a binary version of the abundance to run a logistic binomial regression for probability of occurrence

data <-
  data |>
  dplyr::mutate(abu.tr = log1p(abu),
                abu.bin = dplyr::case_when(abu == 0 ~ 0, TRUE ~ 1),
                abu = round(abu),
                year = as.numeric(year)) 

# Extract some standard information (year range for stations per species)
results.year <- 
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::summarize(
    min.year=min(year),
    max.year=max(year))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Run regressions - Linear regression ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Test linear trends ----

# get slope b of abu ~ year from linear regression
# Use the transformed version of abu as created in the data loading step (abu.tr)
results.abu <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::do(
    broom::tidy(
      lm(abu.tr~year, data = .)
    )) |>
  dplyr::filter(term == "year") |> 
  dplyr::rename(b.linear = estimate,
                SE.b.linear = std.error,
                pval.lin = p.value) |>
  dplyr::select("StationID","species","b.linear","SE.b.linear",
                "pval.lin")

# Get the model statistics using the glance functions
glance.abu <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::do(
    broom::glance(
      lm(abu.tr~year, data = .))) |> 
  rename(AIC.linear = AIC,
         nobs.linear = nobs,
         r2.lin = adj.r.squared) |>
  dplyr::select("StationID","species","AIC.linear","nobs.linear","r2.lin")

# Test model assumptions. The gvlma function from the gvlma package quickly tests for Global Validation of Linear Model Assumptions and provides p-values for linearity
assum.abu <-
  data |>
  dplyr::group_by(StationID, species) |>
  # perform the gvlma function using a list output as the do function does not work here
  dplyr::summarize(
    model = list(gvlma::gvlma(lm(abu.tr ~ year, data = cur_data()))),
    .groups = "drop"
  ) |>
  # return the output in a seperate list 
  dplyr::mutate(gvlma_summary = purrr::map(model, summary)) |>
                # Get output as seperate variables
  tidyr::unnest(cols = gvlma_summary)

# Add test names so we can filter on what we want to keep
assum.abu <-
  assum.abu |>
    dplyr::mutate(Stat = rep(c("Global Stat", "Skewness", "Kurtosis", "Link Function", "Heteroscedasticity"), length.out = nrow(assum.abu))) |>
  # Select and filter the columns we want to keep
  dplyr::select(StationID, species, Decision, Stat) |>
  dplyr::filter(Stat == "Skewness") |>
    # make wide format
    tidyr::pivot_wider(id_cols = c(StationID, species), names_from = Stat, values_from = Decision) |>
    # Add Heteroscedasticity
  dplyr::left_join(assum.abu |>
                     dplyr::mutate(Stat = rep(c("Global Stat", "Skewness", "Kurtosis", "Link Function", "Heteroscedasticity"), length.out = nrow(assum.abu))) |>
                     # Select and filter the columns we want to keep
                     dplyr::select(StationID, species, Decision, Stat) |>
                     dplyr::filter(Stat == "Heteroscedasticity") |>
                     # make wide format
                     tidyr::pivot_wider(id_cols = c(StationID, species), names_from = Stat, values_from = Decision),
                   by = c("StationID", "species"))
  
# combine results together in one dataframe
results.abu <-
  results.year |>
  # Combine the linear model output
  dplyr::left_join(results.abu,
                   by = c("StationID","species")) |>
  # and with the model statistics
  dplyr::left_join(glance.abu,
                   by = c("StationID", "species")) |>
  # And model assumptions
  dplyr::left_join(assum.abu,
                   by = c("StationID", "species"))

## Test non-linear trends (polynomials) ----

# get slope b of abu ~ poly(year,2) from quadratic regression
# Use the transformed version of abu as created in the data loading step (abu.tr)
results.abu2.lin <- 
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::tidy(
      lm(abu.tr~poly(year,2), data = .)
    )) |>
  dplyr::filter(term == "poly(year, 2)1") |>
  dplyr::rename(b.poly = estimate,
                SE.b.poly = std.error,
                pval.b.poly = p.value) |>
  dplyr::select("StationID","species","b.poly","SE.b.poly",
                "pval.b.poly")

# get slope c of abu ~ poly (year,2) from quadratic regression
results.abu2.quad <- 
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::tidy(
      lm(abu.tr~poly(year,2), data = .)
    )) |>
  dplyr::filter(term == "poly(year, 2)2") |> 
  dplyr::rename(c.poly = estimate,
                SE.c.poly = std.error,
                pval.c.poly = p.value) |>
  dplyr::select("StationID","species","c.poly","SE.c.poly",
                "pval.c.poly")

#  Get the model statistics
glance.abu2 <-
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::glance(
      lm(abu.tr~poly(year,2), data = .)
    )) |>
  dplyr::rename(AIC.poly = AIC,
                nobs.poly = nobs,
                r2.poly = adj.r.squared,
                pval.poly = p.value) |>
  dplyr::select("StationID","species","AIC.poly","nobs.poly","r2.poly","pval.poly")

assum.abu2 <-
  data |>
  dplyr::group_by(StationID, species) |>
  # perform the gvlma function using a list output as the do function does not work here
  dplyr::summarize(
    model = list(gvlma::gvlma(lm(abu.tr~poly(year,2), data = cur_data()))),
    .groups = "drop"
  ) |>
  # return the output in a seperate list 
  dplyr::mutate(gvlma_summary = purrr::map(model, summary)) |>
  # Get output as seperate variables
  tidyr::unnest(cols = gvlma_summary)

# Add test names so we can filter on what we want to keep
assum.abu2 <-
  assum.abu2 |>
  dplyr::mutate(Stat = rep(c("Global Stat", "Skewness", "Kurtosis", "Link Function", "Heteroscedasticity"), length.out = nrow(assum.abu2))) |>
  # Select and filter the columns we want to keep
  dplyr::select(StationID, species, Decision, Stat) |>
  dplyr::filter(Stat == "Skewness") |>
  # make wide format
  tidyr::pivot_wider(id_cols = c(StationID, species), names_from = Stat, values_from = Decision) |>
  # Add Heteroscedasticity
  dplyr::left_join(assum.abu2 |>
                     dplyr::mutate(Stat = rep(c("Global Stat", "Skewness", "Kurtosis", "Link Function", "Heteroscedasticity"), length.out = nrow(assum.abu2))) |>
                     # Select and filter the columns we want to keep
                     dplyr::select(StationID, species, Decision, Stat) |>
                     dplyr::filter(Stat == "Heteroscedasticity") |>
                     # make wide format
                     tidyr::pivot_wider(id_cols = c(StationID, species), names_from = Stat, values_from = Decision),
                   by = c("StationID", "species")) |>
  # rename 
  dplyr::rename(Skewness_poly = Skewness,
                Heteroscedasticity_poly = Heteroscedasticity)


# Combine into one dataframe
results.abu2 <-
  # b slope
  results.abu2.lin |>
  # c slope
  dplyr::left_join(results.abu2.quad,
                   by = c("StationID", "species")) |>
  # Model statistics
  dplyr::left_join(glance.abu2,
                   by=c("StationID","species")) |>
  # And model assumptions
  dplyr::left_join(assum.abu2,
                   by = c("StationID", "species"))

# Merge linear and quadratic regression
output <-
  results.abu |>
  dplyr::left_join(results.abu2,
                   by=c("StationID","species"))

## Create the votes ----

output <-
  output |>
  dplyr::mutate(vote = NA) |>
  dplyr::mutate(vote = dplyr::case_when(
    # None of the regressions is significant
    pval.lin >= 0.1 & pval.poly >= 0.1 ~ "no",
    
    # Only linear is significant
    pval.lin < 0.1 & pval.poly >= 0.1 & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.poly >= 0.1 & b.linear > 0 ~ "positive linear",
    
    # Only quadratic is significant
    pval.lin >= 0.1 & pval.poly < 0.1 & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin >= 0.1 & pval.poly < 0.1 & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin >= 0.1 & pval.poly < 0.1 & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin >= 0.1 & pval.poly < 0.1 & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # both are significant, requiring involving AIC, first if AIC is lower for the polynomial
    pval.lin < 0.1 & pval.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin < 0.1 & pval.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin < 0.1 & pval.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin < 0.1 & pval.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # same but when AIC is lower for linear regression
    pval.lin < 0.1 & pval.poly < 0.1 & AIC.linear < AIC.poly & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.poly < 0.1 & AIC.linear < AIC.poly & b.linear > 0 ~ "positive linear",
    pval.lin < 0.1 & is.na(pval.poly) & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & is.na(pval.poly) & b.linear > 0 ~ "positive linear",
    
    # sometimes, the polynomial or even the linear do not produce results, such that the p-value is NA, make these neutral trends as well
    is.na(vote) & pval.lin >= 0.1 & is.na(pval.poly) ~ "no",
    is.na(vote) & is.na(pval.lin) ~ "no"
  )) 

## MOStest ----

## Apply the MOS test on the votes that pass criteria i and ii  
# create unique case identifier for data and for output
data <-
  data |>
  dplyr::mutate(UCI = paste(StationID, species))

output <-
  output |>
  dplyr::mutate(UCI = paste(StationID, species))

UCI <- unique(output$UCI[output$vote=="negative to positive"|
                           output$vote=="positive to negative"])

# as there might be cases when there is no such case, i do the loop as in ifelse clause
if(length(UCI)>0){
  # If cases with species fullfilling criterea i or ii:
  # Create an empty dataframe for output
  MOS.abu<-data.frame()
  for(i in 1:length(UCI)){
    # For each of the cases
    # Create temporary dataframe and extract the data from the data df into it
    temp <- data[data$UCI==UCI[i], ]
    if(dim(temp)[1]>2){
      # Needs more than 2 observations
      # Perform MOStest
      MOS <- vegan::MOStest(temp$year,temp$abu.tr, family = "gaussian")
      MOSout <- MOS$isBracketed
      # Combine output back into MOS.abu
      MOS.abu<- rbind(MOS.abu,data.frame(UCI=temp$UCI[1],MOSout))
      rm(temp)
    }
  }
} else{
  # Otherwise create NA
  MOS.abu<-data.frame(matrix("", ncol = 1, nrow = dim(output)[1]))
  MOS.abu$UCI<-output$UCI
  MOS.abu$MOSout<-NA}

# Inspect created dataframe
head(MOS.abu)

# Combine with original output and change votes if MOS test was FALSE
output <-
  output |>
  dplyr::full_join(MOS.abu |>
                     dplyr::select(UCI,MOSout),
                   by = "UCI") |>
  # change votes if MOS Test was FALSE
  dplyr::mutate(vote = dplyr::case_when(
    vote == "positive to negative" & MOSout == "FALSE" ~ "positive decelerating",
    vote == "negative to positive" & MOSout == "FALSE" ~ "negative decelerating",
    TRUE ~ vote
  ))

summary(as.factor(output$vote))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Run regressions - Poisson regression ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Note that generalized linear models do not output adjusted R2's. it is possible to calculate pseudo R2 but there different methods to do this.
# Do we want these measures??
# Glance function also does not output p-value as a model statistic to be used in the distribution of the votes. Instead here chosen to use the p-value of poly(year,2) instead to compare to the linear trend

## Test linear trends ----

# get slope b of abu ~ year from poisson regression
# Use the the rounded version of abu
results.abu.pois <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::do(
    broom::tidy(
      glm(abu~year, family = "poisson", data = .)
    )) |>
  dplyr::filter(term == "year") |> 
  dplyr::rename(b.linear = estimate,
                SE.b.linear = std.error,
                pval.lin = p.value) |>
  dplyr::select("StationID","species","b.linear","SE.b.linear",
                "pval.lin")

# Get the model statistics using the glance function
glance.abu.pois <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::do(
    broom::glance(
      glm(abu~year, family = "poisson", data = .))) |> 
  rename(AIC.linear = AIC,
         nobs.linear = nobs) |>
  dplyr::select("StationID","species","AIC.linear","nobs.linear")

# Model assumptions
# no simple function like gvlma for a linear model so needs a bit more coding
# source the function to do it in one go
source("scripts/poisson_test.R")
assum.abu.pois <-
  data |>
  dplyr::group_by(StationID, species) |>
  # perform the gvlma function using a list output as the do function does not work here
  dplyr::summarize(
    model = list(check_poisson_assumptions(glm(abu~year, family = "poisson", data = cur_data()))),
    .groups = "drop"
  ) |>
  tidyr::unnest(cols = model)
  
# Return to a wide output again
assum.abu.pois<-
  assum.abu.pois |>
  dplyr::select(-p_value) |>
  dplyr::filter(Test == "Dispersion") |>
  tidyr::pivot_wider(id_cols = c(StationID, species), values_from = Decision, names_from = Test) |>
  dplyr::left_join(assum.abu.pois |>
                     dplyr::select(-p_value) |>
                     dplyr::filter(Test == "Zero Inflation") |>
                     tidyr::pivot_wider(id_cols = c(StationID, species), values_from = Decision, names_from = Test)) |>
  dplyr::left_join(assum.abu.pois |>
                     dplyr::select(-p_value) |>
                     dplyr::filter(Test == "Normality of Residuals") |>
                     tidyr::pivot_wider(id_cols = c(StationID, species), values_from = Decision, names_from = Test))

# combine results together in one dataframe
results.abu.pois <-
  results.year |>
  # Combine the linear model output
  dplyr::left_join(results.abu.pois,
                   by = c("StationID","species")) |>
  # and with the model statistics
  dplyr::left_join(glance.abu.pois,
                   by = c("StationID", "species")) |>
  # Model assumptions
  dplyr::left_join(assum.abu.pois,
                   by = c("StationID", "species"))

## Test non-linear trends (polynomials) ----

# get slope b of abu ~ poly(year,2) from quadratic regression
# Use the the rounded version of abu
results.abu2.lin.pois <- 
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::tidy(
      glm(abu~poly(year,2), family = "poisson", data = .)
    )) |>
  dplyr::filter(term == "poly(year, 2)1") |>
  dplyr::rename(b.poly = estimate,
                SE.b.poly = std.error,
                pval.b.poly = p.value) |>
  dplyr::select("StationID","species","b.poly","SE.b.poly",
                "pval.b.poly")

# get slope c of abu ~ poly (year,2) from quadratic regression
results.abu2.quad.pois <- 
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::tidy(
      glm(abu~poly(year,2), family = "poisson", data = .)
    )) |>
  dplyr::filter(term == "poly(year, 2)2") |> 
  dplyr::rename(c.poly = estimate,
                SE.c.poly = std.error,
                pval.c.poly = p.value) |>
  dplyr::select("StationID","species","c.poly","SE.c.poly",
                "pval.c.poly")

#  Get the model statistics
glance.abu2.pois <-
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::glance(
      glm(abu~poly(year,2), family = "poisson", data = .)
    )) |>
  dplyr::rename(AIC.poly = AIC,
                nobs.poly = nobs) |>
  dplyr::select("StationID","species","AIC.poly","nobs.poly")

# Model assumptions
assum.abu.pois2 <-
  data |>
  dplyr::group_by(StationID, species) |>
  # perform the gvlma function using a list output as the do function does not work here
  dplyr::summarize(
    model = list(check_poisson_assumptions(glm(abu~poly(year,2), family = "poisson", data = cur_data()))),
    .groups = "drop"
  ) |>
  tidyr::unnest(cols = model)

# Return to a wide output again
assum.abu.pois2<-
  assum.abu.pois2 |>
  dplyr::select(-p_value) |>
  dplyr::filter(Test == "Dispersion") |>
  tidyr::pivot_wider(id_cols = c(StationID, species), values_from = Decision, names_from = Test) |>
  dplyr::left_join(assum.abu.pois2 |>
                     dplyr::select(-p_value) |>
                     dplyr::filter(Test == "Zero Inflation") |>
                     tidyr::pivot_wider(id_cols = c(StationID, species), values_from = Decision, names_from = Test)) |>
  dplyr::left_join(assum.abu.pois2 |>
                     dplyr::select(-p_value) |>
                     dplyr::filter(Test == "Normality of Residuals") |>
                     tidyr::pivot_wider(id_cols = c(StationID, species), values_from = Decision, names_from = Test)) |>
  # Change names
  dplyr::rename(Dispersion_poly = "Dispersion",
                `Zero Inflation_poly` = "Zero Inflation",
                `Normality of Residuals_poly` = "Normality of Residuals")

# Combine into one dataframe
results.abu2.pois <-
  # b slope
  results.abu2.lin.pois |>
  # c slope
  dplyr::left_join(results.abu2.quad.pois,
                   by = c("StationID", "species")) |>
  # Model statistics
  dplyr::left_join(glance.abu2.pois,
                   by=c("StationID","species")) |>
  # Model assumptions
  dplyr::left_join(assum.abu.pois2,
                   by = c("StationID", "species"))

# Merge linear and quadratic regression
output.pois <-
  results.abu.pois |>
  dplyr::left_join(results.abu2.pois,
                   by=c("StationID","species"))

## Create the votes ----

output.pois <-
  output.pois |>
  dplyr::mutate(vote = NA) |>
  dplyr::mutate(vote = dplyr::case_when(
    # None of the regressions is significant
    pval.lin >= 0.1 & pval.c.poly >= 0.1 ~ "no",
    
    # Only linear is significant
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear > 0 ~ "positive linear",
    
    # Only quadratic is significant
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # both are significant, requiring involving AIC, first if AIC is lower for the polynomial
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # same but when AIC is lower for linear regression
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear > 0 ~ "positive linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear > 0 ~ "positive linear",
    
    # sometimes, the polynomial or even the linear do not produce results, such that the p-value is NA, make these neutral trends as well
    is.na(vote) & pval.lin >= 0.1 & is.na(pval.c.poly) ~ "no",
    is.na(vote) & is.na(pval.lin) ~ "no"
  )) 

## MOStest ----

## Apply the MOS test on the votes that pass criteria i and ii  
# create unique case identifier for data and for output
data <-
  data |>
  dplyr::mutate(UCI = paste(StationID, species))

output.pois <-
  output.pois |>
  dplyr::mutate(UCI = paste(StationID, species))

UCI <- unique(output.pois$UCI[output.pois$vote=="negative to positive"|
                                output.pois$vote=="positive to negative"])

# as there might be cases when there is no such case, i do the loop as in ifelse clause
if(length(UCI)>0){
  # If cases with species fullfilling criterea i or ii:
  # Create an empty dataframe for output
  MOS.abu.pois<-data.frame()
  for(i in 1:length(UCI)){
    # For each of the cases
    # Create temporary dataframe and extract the data from the data df into it
    temp <- data[data$UCI==UCI[i], ]
    if(dim(temp)[1]>2){
      # Needs more than 2 observations
      # Perform MOStest with error handling
      MOS <- tryCatch({
        vegan::MOStest(temp$year,temp$abu, family = "poisson")
      }, error = function(e) {
        return(list(isBracketed = NA))
      })
      MOSout <- MOS$isBracketed
      # Combine output back into MOS.abu.pois
      MOS.abu.pois<- rbind(MOS.abu.pois,data.frame(UCI=temp$UCI[1],MOSout))
      rm(temp)
    }
  }
} else{
  # Otherwise create NA
  MOS.abu.pois<-data.frame(matrix("", ncol = 1, nrow = dim(output.pois)[1]))
  MOS.abu.pois$UCI<-output.pois$UCI
  MOS.abu.pois$MOSout<-NA}

# Inspect created dataframe
head(MOS.abu.pois)

# Combine with original output and change votes if MOS test was FALSE
output.pois <-
  output.pois |>
  dplyr::full_join(MOS.abu.pois |>
                     dplyr::select(UCI,MOSout),
                   by = "UCI") |>
  # change votes if MOS Test was FALSE
  dplyr::mutate(vote = dplyr::case_when(
    vote == "positive to negative" & MOSout == "FALSE" ~ "positive decelerating",
    vote == "negative to positive" & MOSout == "FALSE" ~ "negative decelerating",
    TRUE ~ vote
  ))

summary(as.factor(output.pois$vote))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Run regressions - Negative binomial regression ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## In case of overdisperison a negative binomial model could be used instead of a poisson model

# get slope b of abu ~ year from poisson regression
# Use the the rounded version of abu

# The negative binomial model tends to throw errors, stopping the model run. Requires the ommision of errors
# Define a safe version of your model function
safe_glm.nb <- possibly(function(df) broom::tidy(glm.nb(abu~year, data = df)), otherwise = tibble(term = "year"))

results.abu.nb <-
  data |>
  dplyr::group_by(StationID, species) |>
  # The negative binomial model tends to throw errors, stopping the model run. Requires the ommision of errors
  dplyr::summarise(model_results = list(safe_glm.nb(cur_data())), .groups = "drop") |>
  # Filter for the year term
  dplyr::mutate(model_results = purrr::map(model_results, ~dplyr::filter(., term == "year"))) |>
  tidyr::unnest(model_results) |>
  dplyr::rename(b.linear = estimate,
                SE.b.linear = std.error,
                pval.lin = p.value) |>
  dplyr::select("StationID","species","b.linear","SE.b.linear",
                "pval.lin")

# Define a safe function for glance
safe_glm.nb <- possibly(function(df) broom::glance(glm.nb(abu~year, data = df)), otherwise = tibble(null.deviance = NA, df.null = NA, logLik = NA, AIC = NA, BIC = NA, deviance = NA, df.residual = NA, nobs = NA))

# Get the model statistics using the glance function
glance.abu.nb <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::summarise(model_results = list(safe_glm.nb(cur_data())), .groups = "drop") |>
  # Remove the logLik column, it causes errors
  dplyr::mutate(model_results = purrr::map(model_results, ~ .x %>% dplyr::select(-logLik))) |>
  # Filter for the year term
  tidyr::unnest(model_results) |>
  rename(AIC.linear = AIC,
         nobs.linear = nobs) |>
  dplyr::select("StationID","species","AIC.linear","nobs.linear")


# combine results together in one dataframe
results.abu.nb <-
  results.year |>
  # Combine the linear model output
  dplyr::left_join(results.abu.nb,
                   by = c("StationID","species")) |>
  # and with the model statistics
  dplyr::left_join(glance.abu.nb,
                   by = c("StationID", "species")) 

## Test non-linear trends (polynomials) ----

# get slope b of abu ~ poly(year,2) from quadratic regression
# Use the the rounded version of abu
# Define a safe version of your model function
safe_glm.nb <- possibly(function(df) broom::tidy(glm.nb(abu~poly(year,2), data = df)), otherwise = tibble(term = "year"))


results.abu2.lin.nb <- 
  data |>
  dplyr::group_by(StationID, species) |>
  # The negative binomial model tends to throw errors, stopping the model run. Requires the ommision of errors
  dplyr::summarise(model_results = list(safe_glm.nb(cur_data())), .groups = "drop") |>
  # Filter for the year term
  dplyr::mutate(model_results = purrr::map(model_results, ~dplyr::filter(., term == "poly(year, 2)1"))) |>
  tidyr::unnest(model_results) |>
  dplyr::rename(b.poly = estimate,
                SE.b.poly = std.error,
                pval.b.poly = p.value) |>
  dplyr::select("StationID","species","b.poly","SE.b.poly",
                "pval.b.poly")

# get slope c of abu ~ poly (year,2) from quadratic regression
results.abu2.quad.nb <- 
  data |>
  dplyr::group_by(StationID, species) |>
  # The negative binomial model tends to throw errors, stopping the model run. Requires the ommision of errors
  dplyr::summarise(model_results = list(safe_glm.nb(cur_data())), .groups = "drop") |>
  # Filter for the year term
  dplyr::mutate(model_results = purrr::map(model_results, ~dplyr::filter(., term == "poly(year, 2)2")))  |> 
  tidyr::unnest() |>
  dplyr::rename(c.poly = estimate,
                SE.c.poly = std.error,
                pval.c.poly = p.value) |>
  dplyr::select("StationID","species","c.poly","SE.c.poly",
                "pval.c.poly")

#  Get the model statistics
# Define a safe function for glance
safe_glm.nb <- possibly(function(df) broom::glance(glm.nb(abu~poly(year,2), data = df)), otherwise = tibble(null.deviance = NA, df.null = NA, logLik = NA, AIC = NA, BIC = NA, deviance = NA, df.residual = NA, nobs = NA))

glance.abu2.nb <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::summarise(model_results = list(safe_glm.nb(cur_data())), .groups = "drop") |>
  # Remove the logLik column, it causes errors
  dplyr::mutate(model_results = purrr::map(model_results, ~ .x %>% dplyr::select(-logLik))) |>
  tidyr::unnest(model_results) |>
  dplyr::rename(AIC.poly = AIC,
                nobs.poly = nobs) |>
  dplyr::select("StationID","species","AIC.poly","nobs.poly")

# Combine into one dataframe
results.abu2.nb <-
  # b slope
  results.abu2.lin.nb |>
  # c slope
  dplyr::left_join(results.abu2.quad.nb,
                   by = c("StationID", "species")) |>
  # Model statistics
  dplyr::left_join(glance.abu2.nb,
                   by=c("StationID","species"))

# Merge linear and quadratic regression
output.nb <-
  results.abu.nb |>
  dplyr::left_join(results.abu2.nb,
                   by=c("StationID","species"))

## Create the votes ----

output.nb <-
  output.nb |>
  dplyr::mutate(vote = NA) |>
  dplyr::mutate(vote = dplyr::case_when(
    # None of the regressions is significant
    pval.lin >= 0.1 & pval.c.poly >= 0.1 ~ "no",
    
    # Only linear is significant
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear > 0 ~ "positive linear",
    
    # Only quadratic is significant
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # both are significant, requiring involving AIC, first if AIC is lower for the polynomial
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # same but when AIC is lower for linear regression
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear > 0 ~ "positive linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear > 0 ~ "positive linear",
    
    # sometimes, the polynomial or even the linear do not produce results, such that the p-value is NA, make these neutral trends as well
    is.na(vote) & pval.lin >= 0.1 & is.na(pval.c.poly) ~ "no",
    is.na(vote) & is.na(pval.lin) ~ "no"
  )) 

## MOStest ----

## Not available for negative binomial models... 

# create unique case identifier for data and for output
data <-
  data |>
  dplyr::mutate(UCI = paste(StationID, species))

output.nb <-
  output.nb |>
  dplyr::mutate(UCI = paste(StationID, species))

UCI <- unique(output.nb$UCI[output.nb$vote=="negative to positive"|
                                output.nb$vote=="positive to negative"])

summary(as.factor(output.nb$vote))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Run regressions - Binomial regression ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Note that generalized linear models do not output adjusted R2's. it is possible to calculate pseudo R2 but there different methods to do this.
# Do we want these measures??
# Glance function also does not output p-value as a model statistic to be used in the distribution of the votes. Instead here chosen to use the p-value of poly(year,2) instead to compare to the linear trend

## Test linear trends ----

# get slope b of abu ~ year from binomial regression
# Use the the binary version of abu
results.abu.bin <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::do(
    broom::tidy(
      glm(abu.bin~year, family = "binomial", data = .)
    )) |>
  dplyr::filter(term == "year") |> 
  dplyr::rename(b.linear = estimate,
                SE.b.linear = std.error,
                pval.lin = p.value) |>
  dplyr::select("StationID","species","b.linear","SE.b.linear",
                "pval.lin")

# Get the model statistics using the glance function
glance.abu.bin <-
  data |>
  dplyr::group_by(StationID, species) |>
  dplyr::do(
    broom::glance(
      glm(abu.bin~year, family = "binomial", data = .))) |> 
  rename(AIC.linear = AIC,
         nobs.linear = nobs) |>
  dplyr::select("StationID","species","AIC.linear","nobs.linear")

# combine results together in one dataframe
results.abu.bin <-
  results.year |>
  # Combine the linear model output
  dplyr::left_join(results.abu.bin,
                   by = c("StationID","species")) |>
  # and with the model statistics
  dplyr::left_join(glance.abu.bin,
                   by = c("StationID", "species"))

## Test non-linear trends (polynomials) ----

# get slope b of abu ~ poly(year,2) from quadratic regression
# Use the transformed binary version of abu (abu.bin)
results.abu2.lin.bin <- 
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::tidy(
      glm(abu.bin~poly(year,2), family = "binomial", data = .)
    )) |>
  dplyr::filter(term == "poly(year, 2)1") |>
  dplyr::rename(b.poly = estimate,
                SE.b.poly = std.error,
                pval.b.poly = p.value) |>
  dplyr::select("StationID","species","b.poly","SE.b.poly",
                "pval.b.poly")

# get slope c of abu ~ poly (year,2) from quadratic regression
results.abu2.quad.bin <- 
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::tidy(
      glm(abu.bin~poly(year,2), family = "binomial", data = .)
    )) |>
  dplyr::filter(term == "poly(year, 2)2") |> 
  dplyr::rename(c.poly = estimate,
                SE.c.poly = std.error,
                pval.c.poly = p.value) |>
  dplyr::select("StationID","species","c.poly","SE.c.poly",
                "pval.c.poly")

#  Get the model statistics
glance.abu2.bin <-
  data |>
  dplyr::group_by(StationID,species) |>
  dplyr::do(
    broom::glance(
      glm(abu.bin~poly(year,2), family = "binomial", data = .)
    )) |>
  dplyr::rename(AIC.poly = AIC,
                nobs.poly = nobs) |>
  dplyr::select("StationID","species","AIC.poly","nobs.poly")

# Combine into one dataframe
results.abu2.bin <-
  # b slope
  results.abu2.lin.bin |>
  # c slope
  dplyr::left_join(results.abu2.quad.bin,
                   by = c("StationID", "species")) |>
  # Model statistics
  dplyr::left_join(glance.abu2.bin,
                   by=c("StationID","species"))

# Merge linear and quadratic regression
output.bin <-
  results.abu.bin |>
  dplyr::left_join(results.abu2.bin,
                   by=c("StationID","species"))

## Create the votes ----

output.bin <-
  output.bin |>
  dplyr::mutate(vote = NA) |>
  dplyr::mutate(vote = dplyr::case_when(
    # None of the regressions is significant
    pval.lin >= 0.1 & pval.c.poly >= 0.1 ~ "no",
    
    # Only linear is significant
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear > 0 ~ "positive linear",
    
    # Only quadratic is significant
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # both are significant, requiring involving AIC, first if AIC is lower for the polynomial
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # same but when AIC is lower for linear regression
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear > 0 ~ "positive linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear > 0 ~ "positive linear",
    
    # sometimes, the polynomial or even the linear do not produce results, such that the p-value is NA, make these neutral trends as well
    is.na(vote) & pval.lin >= 0.1 & is.na(pval.c.poly) ~ "no",
    is.na(vote) & is.na(pval.lin) ~ "no"
  )) 

## MOStest ----

## Apply the MOS test on the votes that pass criteria i and ii  
# create unique case identifier for data and for output
data <-
  data |>
  dplyr::mutate(UCI = paste(StationID, species))

output.bin <-
  output.bin |>
  dplyr::mutate(UCI = paste(StationID, species))

UCI <- unique(output.bin$UCI[output.bin$vote=="negative to positive"|
                               output.bin$vote=="positive to negative"])

# as there might be cases when there is no such case, i do the loop as in ifelse clause
if(length(UCI)>0){
  # If cases with species fullfilling criterea i or ii:
  # Create an empty dataframe for output
  MOS.abu.bin<-data.frame()
  for(i in 1:length(UCI)){
    # For each of the cases
    # Create temporary dataframe and extract the data from the data df into it
    temp <- data[data$UCI==UCI[i], ]
    if(dim(temp)[1]>2){
      # Needs more than 2 observations
      # Perform MOStest
      MOS <- vegan::MOStest(temp$year,temp$abu.bin, family = "binomial")
      MOSout <- MOS$isBracketed
      # Combine output back into MOS.abu.bin
      MOS.abu.bin<- rbind(MOS.abu.bin,data.frame(UCI=temp$UCI[1],MOSout))
      rm(temp)
    }
  }
} else{
  # Otherwise create NA
  MOS.abu.bin<-data.frame(matrix("", ncol = 1, nrow = dim(output.bin)[1]))
  MOS.abu.bin$UCI<-output.bin$UCI
  MOS.abu.bin$MOSout<-NA}

# Inspect created dataframe
head(MOS.abu.bin)

# Combine with original output and change votes if MOS test was FALSE
output.bin <-
  output.bin |>
  dplyr::full_join(MOS.abu.bin |>
                     dplyr::select(UCI,MOSout),
                   by = "UCI") |>
  # change votes if MOS Test was FALSE
  dplyr::mutate(vote = dplyr::case_when(
    vote == "positive to negative" & MOSout == "FALSE" ~ "positive decelerating",
    vote == "negative to positive" & MOSout == "FALSE" ~ "negative decelerating",
    TRUE ~ vote
  ))

summary(as.factor(output.bin$vote)) 

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Combine all the votes and make a decision tree for the model output to use ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# Combine output files with some renaming
all_output<- 
  output |>
  dplyr::left_join(
    output.pois |>
      dplyr::select(StationID, species, 
                    b.linear, SE.b.linear, pval.lin, AIC.linear, nobs.linear, Dispersion, `Zero Inflation`, `Normality of Residuals`,
                    b.poly, SE.b.poly, pval.b.poly, c.poly, SE.c.poly, pval.c.poly, AIC.poly, nobs.poly, Dispersion_poly, `Zero Inflation_poly`, `Normality of Residuals_poly`) |>
      dplyr::rename_with(~paste0(., "_pois")) |>
      dplyr::rename(StationID = StationID_pois,
                    species = species_pois)
  ) |>
  dplyr::left_join(
    output.nb |>
      dplyr::select(StationID, species, 
                    b.linear, SE.b.linear, pval.lin, AIC.linear, nobs.linear, 
                    b.poly, SE.b.poly, pval.b.poly, c.poly, SE.c.poly, pval.c.poly, AIC.poly, nobs.poly) |>
      dplyr::rename_with(~paste0(., "_nb")) |>
      dplyr::rename(StationID = StationID_nb,
                    species = species_nb)
  ) |>
  dplyr::ungroup()

# make a decision on which model output to use in each case

all_output <-
  all_output |>
  dplyr::mutate(
    # Model chosen
    best_linear = dplyr::case_when(
      Skewness == "Assumptions acceptable." & Heteroscedasticity == "Assumptions acceptable." & AIC.linear < AIC.linear_pois & AIC.linear < AIC.linear_nb ~ "lm",
      Skewness == "Assumptions acceptable." & Heteroscedasticity == "Assumptions acceptable." & AIC.linear > AIC.linear_pois & Dispersion_pois == "No overdispersion detected" & `Zero Inflation_pois` == "No Zero-inflation detected" & `Normality of Residuals_pois` == "Residuals normal" ~ "poisson", 
      Skewness == "Assumptions acceptable." & Heteroscedasticity == "Assumptions acceptable." & AIC.linear > AIC.linear_nb & Dispersion_pois == "Overdispersion detected" ~ "negative binomial",
      Skewness == "Assumptions acceptable." & Heteroscedasticity == "Assumptions acceptable." & AIC.linear > AIC.linear_nb & AIC.linear < AIC.linear_pois & Dispersion_pois == "No overdispersion detected" ~ "lm",
      
      Skewness == "Assumptions NOT satisfied!" & Heteroscedasticity == "Assumptions acceptable." & Dispersion_pois == "No overdispersion detected" & `Zero Inflation_pois` == "No Zero-inflation detected" & `Normality of Residuals_pois` == "Residuals normal" ~ "poisson",
      Skewness == "Assumptions acceptable." & Heteroscedasticity == "Assumptions NOT satisfied!" & Dispersion_pois == "No overdispersion detected" & `Zero Inflation_pois` == "No Zero-inflation detected" & `Normality of Residuals_pois` == "Residuals normal" ~ "poisson",
      Skewness == "Assumptions NOT satisfied!" & Heteroscedasticity == "Assumptions NOT satisfied!" & Dispersion_pois == "No overdispersion detected" & `Zero Inflation_pois` == "No Zero-inflation detected" & `Normality of Residuals_pois` == "Residuals normal" ~ "poisson",
      
      Skewness == "Assumptions NOT satisfied!" | Heteroscedasticity == "Assumptions NOT satisfied!" & Dispersion_pois == "No overdispersion detected" & `Zero Inflation_pois` == "Zero-inflation detected" | `Normality of Residuals_pois` == "Residuals non-normal" ~ "No best model",
      
      Skewness == "Assumptions acceptable." & Heteroscedasticity == "Assumptions NOT satisfied!" & Dispersion_pois == "Overdispersion detected" & AIC.linear_nb < AIC.linear_pois ~ "negative binomial", 
      Skewness == "Assumptions NOT satisfied!" & Heteroscedasticity == "Assumptions NOT satisfied!" & Dispersion_pois == "Overdispersion detected" & AIC.linear_nb < AIC.linear_pois ~ "negative binomial", 
      Skewness == "Assumptions NOT satisfied!" & Heteroscedasticity == "Assumptions acceptable." & Dispersion_pois == "Overdispersion detected" & AIC.linear_nb < AIC.linear_pois ~ "negative binomial", 
    TRUE ~ "No best model")
  ) |> 
  dplyr::mutate(
    # Model chosen
    best_poly = dplyr::case_when(
      Skewness_poly == "Assumptions acceptable." & Heteroscedasticity_poly == "Assumptions acceptable." & AIC.poly < AIC.poly_pois & AIC.poly < AIC.poly_nb ~ "lm",
      Skewness_poly == "Assumptions acceptable." & Heteroscedasticity_poly == "Assumptions acceptable." & AIC.poly > AIC.poly_pois & Dispersion_poly_pois == "No overdispersion detected" & `Zero Inflation_poly_pois` == "No Zero-inflation detected" & `Normality of Residuals_poly_pois` == "Residuals normal" ~ "poisson", 
      Skewness_poly == "Assumptions acceptable." & Heteroscedasticity_poly == "Assumptions acceptable." & AIC.poly > AIC.poly_nb & Dispersion_poly_pois == "Overdispersion detected" ~ "negative binomial",
      Skewness_poly == "Assumptions acceptable." & Heteroscedasticity_poly == "Assumptions acceptable." & AIC.poly > AIC.poly_nb & AIC.poly < AIC.poly_pois & Dispersion_poly_pois == "No overdispersion detected" ~ "lm",
      
      Skewness_poly == "Assumptions NOT satisfied!" & Heteroscedasticity_poly == "Assumptions acceptable." & Dispersion_poly_pois == "No overdispersion detected" & `Zero Inflation_poly_pois` == "No Zero-inflation detected" & `Normality of Residuals_poly_pois` == "Residuals normal" ~ "poisson",
      Skewness_poly == "Assumptions acceptable." & Heteroscedasticity_poly == "Assumptions NOT satisfied!" & Dispersion_poly_pois == "No overdispersion detected" & `Zero Inflation_poly_pois` == "No Zero-inflation detected" & `Normality of Residuals_poly_pois` == "Residuals normal" ~ "poisson",
      Skewness_poly == "Assumptions NOT satisfied!" & Heteroscedasticity_poly == "Assumptions NOT satisfied!" & Dispersion_poly_pois == "No overdispersion detected" & `Zero Inflation_poly_pois` == "No Zero-inflation detected" & `Normality of Residuals_poly_pois` == "Residuals normal" ~ "poisson",
      
      Skewness_poly == "Assumptions NOT satisfied!" | Heteroscedasticity_poly == "Assumptions NOT satisfied!" & Dispersion_poly_pois == "No overdispersion detected" & `Zero Inflation_poly_pois` == "Zero-inflation detected" | `Normality of Residuals_poly_pois` == "Residuals non-normal" ~ "No best model",
      
      Skewness_poly == "Assumptions acceptable." & Heteroscedasticity_poly == "Assumptions NOT satisfied!" & Dispersion_poly_pois == "Overdispersion detected" & AIC.poly_nb < AIC.poly_pois ~ "negative binomial", 
      Skewness_poly == "Assumptions NOT satisfied!" & Heteroscedasticity_poly == "Assumptions NOT satisfied!" & Dispersion_poly_pois == "Overdispersion detected" & AIC.poly_nb < AIC.poly_pois ~ "negative binomial", 
      Skewness_poly == "Assumptions NOT satisfied!" & Heteroscedasticity_poly == "Assumptions acceptable." & Dispersion_poly_pois == "Overdispersion detected" & AIC.poly_nb < AIC.poly_pois ~ "negative binomial", 
      TRUE ~ "No best model"
    )
  ) 
    
# Take the output of the correct model for determining trends
all_output <- 
  all_output |>
  # Linear model
  dplyr::mutate(
    b.linear = dplyr::case_when(
      best_linear == "lm" ~ b.linear,
      best_linear == "poisson" ~ b.linear_pois,
      best_linear == "negative binomial" ~ b.linear_nb,
      best_linear == "No best model" ~ NA
    ),
    SE.b.linear = dplyr::case_when(
      best_linear == "lm" ~ SE.b.linear,
      best_linear == "poisson" ~ SE.b.linear_pois,
      best_linear == "negative binomial" ~ SE.b.linear_nb,
      best_linear == "No best model" ~ NA
    ),
    pval.lin = dplyr::case_when(
      best_linear == "lm" ~ pval.lin,
      best_linear == "poisson" ~ pval.lin_pois,
      best_linear == "negative binomial" ~ pval.lin_nb,
      best_linear == "No best model" ~ NA
    ),
    AIC.linear = dplyr::case_when(
      best_linear == "lm" ~ AIC.linear,
      best_linear == "poisson" ~ AIC.linear_pois,
      best_linear == "negative binomial" ~ AIC.linear_nb,
      best_linear == "No best model" ~ NA
    ),
    nobs.linear = dplyr::case_when(
      best_linear == "lm" ~ nobs.linear,
      best_linear == "poisson" ~ nobs.linear_pois,
      best_linear == "negative binomial" ~ nobs.linear_nb,
      best_linear == "No best model" ~ NA
    ),
    r2.lin = dplyr::case_when(
      best_linear == "lm" ~ r2.lin,
      best_linear == "poisson" ~ NA,
      best_linear == "negative binomial" ~ NA,
      best_linear == "No best model" ~ NA
    )) |>
  # Polynomial
  dplyr::mutate(
    b.poly = dplyr::case_when(
      best_poly == "lm" ~ b.poly,
      best_poly == "poisson" ~ b.poly_pois,
      best_poly == "negative binomial" ~ b.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    SE.b.poly = dplyr::case_when(
      best_poly == "lm" ~ SE.b.poly,
      best_poly == "poisson" ~ SE.b.poly_pois,
      best_poly == "negative binomial" ~ SE.b.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    pval.b.poly = dplyr::case_when(
      best_poly == "lm" ~ pval.b.poly,
      best_poly == "poisson" ~ pval.b.poly_pois,
      best_poly == "negative binomial" ~ pval.b.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    c.poly = dplyr::case_when(
      best_poly == "lm" ~ c.poly,
      best_poly == "poisson" ~ c.poly_pois,
      best_poly == "negative binomial" ~ c.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    SE.c.poly = dplyr::case_when(
      best_poly == "lm" ~ SE.c.poly,
      best_poly == "poisson" ~ SE.c.poly_pois,
      best_poly == "negative binomial" ~ SE.c.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    pval.c.poly = dplyr::case_when(
      best_poly == "lm" ~ pval.c.poly,
      best_poly == "poisson" ~ pval.c.poly_pois,
      best_poly == "negative binomial" ~ pval.c.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    AIC.poly = dplyr::case_when(
      best_poly == "lm" ~ AIC.poly,
      best_poly == "poisson" ~ AIC.poly_pois,
      best_poly == "negative binomial" ~ AIC.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    nobs.poly = dplyr::case_when(
      best_poly == "lm" ~ nobs.poly,
      best_poly == "poisson" ~ nobs.poly_pois,
      best_poly == "negative binomial" ~ nobs.poly_nb,
      best_poly == "No best model" ~ NA
    ),
    r2.poly = dplyr::case_when(
      best_poly == "lm" ~ r2.poly,
      best_poly == "poisson" ~ NA,
      best_poly == "negative binomial" ~ NA,
      best_poly == "No best model" ~ NA
    )
  ) |>
  # clean up the dataframe
  dplyr::select(StationID, species, min.year, max.year,
                b.linear, SE.b.linear, pval.lin, AIC.linear, nobs.linear, r2.lin, best_linear,
                b.poly, SE.b.poly, pval.b.poly, c.poly, SE.c.poly, pval.c.poly, AIC.poly, nobs.poly, r2.poly, best_poly)

# Create votes
all_output <-
  all_output |>
  dplyr::mutate(vote = NA) |>
  dplyr::mutate(vote = dplyr::case_when(
    # None of the regressions is significant
    pval.lin >= 0.1 & pval.c.poly >= 0.1 ~ "no",
    
    # Only linear is significant
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly >= 0.1 & b.linear > 0 ~ "positive linear",
    
    # Only quadratic is significant
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin >= 0.1 & pval.c.poly < 0.1 & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # both are significant, requiring involving AIC, first if AIC is lower for the polynomial
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly < 0 ~ "negative accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly > 0 ~ "positive accelerating",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly < 0 & c.poly > 0 ~ "negative to positive",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear > AIC.poly & b.poly > 0 & c.poly < 0 ~ "positive to negative",
    
    # same but when AIC is lower for linear regression
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & pval.c.poly < 0.1 & AIC.linear < AIC.poly & b.linear > 0 ~ "positive linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear < 0 ~ "negative linear",
    pval.lin < 0.1 & is.na(pval.c.poly) & b.linear > 0 ~ "positive linear",
    
    # sometimes, the polynomial or even the linear do not produce results, such that the p-value is NA, make these neutral trends as well
    is.na(vote) & pval.lin >= 0.1 & is.na(pval.c.poly) ~ "no",
    is.na(vote) & is.na(pval.lin) ~ "no"
  )) 


## MOStest ----
# Perform MOS test. depending on the model used different distributions are possible

## Apply the MOS test on the votes that pass criteria i and ii  
# create unique case identifier for data and for output
data <-
  data |>
  dplyr::mutate(UCI = paste(StationID, species))

all_output <-
  all_output |>
  dplyr::mutate(UCI = paste(StationID, species))

UCI <- unique(all_output$UCI[all_output$vote=="negative to positive"|
                               all_output$vote=="positive to negative"])
MOD <- all_output$best_poly[all_output$vote=="negative to positive"|
                              all_output$vote=="positive to negative"]

# as there might be cases when there is no such case, i do the loop as in ifelse clause
if(length(UCI)>0){
  # If cases with species fullfilling criterea i or ii:
  # Create an empty dataframe for output
  MOS.abu<-data.frame()
  for(i in 1:length(UCI)){
    # For each of the cases
    # Create temporary dataframe and extract the data from the data df into it
    temp <- data[data$UCI==UCI[i], ]
    temp$MOD <- MOD[i]
    if(dim(temp)[1]>2){
      if(temp$MOD[i] == "lm"){
        # Needs more than 2 observations
        # Perform MOStest
        MOS <- vegan::MOStest(temp$year,temp$abu.tr, family = "gaussian")
        MOSout <- MOS$isBracketed
        # Combine output back into MOS.abu
        MOS.abu<- rbind(MOS.abu,data.frame(UCI=temp$UCI[1],MOSout))
      } else{
      if(temp$MOD[i] == "poisson"){
      # Needs more than 2 observations
      # Perform MOStest with error handling
      MOS <- tryCatch({
        vegan::MOStest(temp$year,temp$abu, family = "poisson")
      }, error = function(e) {
        return(list(isBracketed = NA))
      })
      MOSout <- MOS$isBracketed
      # Combine output back into MOS.abu.pois
      MOS.abu<- rbind(MOS.abu,data.frame(UCI=temp$UCI[1],MOSout))
      rm(temp)
      }
      }
  } else{
  # Otherwise create NA
  MOS.abu<-data.frame(matrix("", ncol = 1, nrow = dim(all_output)[1]))
  MOS.abu$UCI<-all_output$UCI
  MOS.abu$MOSout<-NA}}
}

# Inspect created dataframe
head(MOS.abu)

# Combine with original output and change votes if MOS test was FALSE
all_output <-
  all_output |>
  dplyr::full_join(MOS.abu |>
                     dplyr::select(UCI,MOSout),
                   by = "UCI") |>
  # change votes if MOS Test was FALSE
  dplyr::mutate(vote = dplyr::case_when(
    vote == "positive to negative" & MOSout == "FALSE" ~ "positive decelerating",
    vote == "negative to positive" & MOSout == "FALSE" ~ "negative decelerating",
    TRUE ~ vote
  ))

summary(as.factor(all_output$vote))


# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Output for meta analysis ----
# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# For now assuming the output of the linear models with the transformation
# Rename the output (otherwise we overwrite each other)

phyto1 <- 
  all_output |>
  # the data set needs a unique identifier (see above, will use this as random term)
  # suggest do use initials XXX and numbers that you can then allocate yourself
  dplyr::mutate(UDI = "HLH001",
                # Information on the organism group
                organism.group = "phytoplankton", # entries can be: mammals, birds, macroinvertebrates, fish, zooplankton, phytoplankton, macrophytes, bacteria
                #Information whether species-biomass or -abundance were used
                measure = "biomass",
                transformation = "log1p" # transformation applied to the data for the lm model
                ) 



# Information on the phylogeny
# In the exemplary case this can be taken from the original data
phylo <- read_delim("data/PPKTcount_noHT_28092022.csv", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)
phylo <- 
  phylo |> 
  dplyr::rename(species = Species) |>
  dplyr::select(Phylum, Class, Order, Family, Genus, species) |>
  dplyr::distinct()

phyto1 <- dplyr::left_join(phyto1,phylo,by="species")

# If the phylogeny is not available, you need to add all five variables via mydata$Phylum<-NA, same for Class, Order, Family, Genus. Instead of NA, you can if you have also add text, e.g. mydata$Phylum = "Mollusca"

# Finally, a few variables for the data origin
# where did the data come from
phyto1$origin<-"Rijkswaterstaat NL"
phyto1$origin[phyto1$StationID=="Nney_W_2"]<-"NLWKN"
phyto1$origin[phyto1$StationID=="JaBu_W_1"]<-"NLWKN"
phyto1$origin[phyto1$StationID=="Bork_W_1"]<-"NLWKN"
phyto1$origin[phyto1$StationID=="WeMu_W_1"]<-"NLWKN"
# who entered them into our synthesis
phyto1$enter<-"Helmut Hillebrand"
#is there a DOI where we can read more on this data
phyto1$DOI<-"10.1007/s12526-023-01382-9"

#We can discuss at a later point if we also want to map things, then we need to add lat and long for the StationID, currently I only add placeholders
phyto1$lat<-NA
phyto1$long<-NA

# The final data set has the following names
names(phyto1)
# and needs to be saved like this
write.csv(phyto1,file="output/phylo1.csv")