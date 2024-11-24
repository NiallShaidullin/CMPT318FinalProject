library(dplyr)
library(corrplot)
library(ggplot2)
library(lubridate)
library(zoo)
library(imputeTS)
library("depmixS4")
library(ggbiplot)
library(factoextra)
library(plotly)
set.seed(1)

### Functions ###
# add a specified number of seconds to a POSIXlt time object
addTimeToPOSIXlt <- function(time, numSecs) {
  return(as.POSIXlt(as.POSIXct(time)+numSecs))
}

# Convert POSIXlt to POSIXct
POSIXlt_to_ct <- function(datetime, format_str) {
  return(as.POSIXct(strftime(datetime, format = format_str), format = format_str, tz = timezone))
}

# Calculate the rolling mean for the 'Global_intensity' column with a defined window size
rollMeanGlobInt <- function(df, window) {
  rollmean(df[['Global_intensity']], window)
}

# Function to train and evaluate HMM models
train_evaluate_hmm <- function(train_data, test_data, states_range, train_n_times, test_n_times) {
  results <- list() # Store log-likelihoods and BIC for comparison
  
  for (num_states in states_range) {
    # Step 1: Define the HMM model for training data
    model <- depmix(
      response = list(Global_intensity ~ 1, Global_active_power ~ 1),
      data = train_data,
      nstates = num_states,
      family = list(gaussian(), gaussian()),
      ntimes = train_n_times,
      emcontrol = list(maxit = 500)
    )
    
    # Step 2: Fit the model on the training data
    fitted_model <- fit(model)

    # Step 3: Extract log-likelihood for training data and normalize
    log_likelihood_train <- logLik(fitted_model)
    normalized_train_log_likelihood <- log_likelihood_train / 107  # Normalize by total observations
    
    # Step 4: Extract BIC for the training model
    bic <- BIC(fitted_model)
    print(bic)
    # Step 5: Create a similar model for the test data using the same structure
    # This is done using the same specification as the trained model, but with test data
    new_mod <- depmix(
      response = list(Global_intensity ~ 1, Global_active_power ~ 1),
      data = test_data,
      nstates = num_states,
      family = list(gaussian(), gaussian()),
      ntimes = test_n_times,
      emcontrol = list(maxit = 500)
      )
    
    # Step 6: Set the parameters of the new model to be the same as the trained model
    model2 <- setpars(new_mod, getpars(fitted_model))  # Set trained parameters into the test model
    
    # Step 7: Evaluate the log-likelihood of the test data on the trained model
    log_likelihood_test <- logLik(model2)  # Calculate log-likelihood of the test data
    # Check if it's a valid number (not NaN or Inf)
    print(log_likelihood_test)
    
    # Step 8: Normalize the test log-likelihood
    normalized_test_log_likelihood <- log_likelihood_test / 48  # Normalize by total observations
    
    # Step 9: Store the results for this model configuration
    results[[as.character(num_states)]] <- list(
      num_states = num_states,
      train_log_likelihood = as.numeric(log_likelihood_train),
      test_log_likelihood = as.numeric(log_likelihood_test),
      bic = bic
    )
  }
  
  return(results)
}



timezone <- "UTC"
window_size <- 10 # moving time window

# Loading data 
data <- read.table("TermProjectData.txt", header = TRUE, sep = ",")

# see how many na values for each feature
colSums(is.na(data))

# go through all the data, if it is numeric then interpolate the na values (From Assignment1)
data <- as.data.frame(lapply(data, function(col) {
  if (is.numeric(col)) {
    return(na.approx(col, na.rm = FALSE)) # Apply interpolation
  } 
  else {
    return(col) # Leave non-numeric columns as is
  }
}))
# Replace the first NA in Global_active_power with the median
data$Global_active_power[1] <- median(data$Global_active_power, na.rm = TRUE)

# check that there are no na values left
colSums(is.na(data))

columns = c("Global_active_power","Global_reactive_power","Voltage","Global_intensity","Sub_metering_1", "Sub_metering_2", "Sub_metering_3" )

# Create a new dataframe to store z-scores
z_scores <- data.frame(matrix(nrow = nrow(data), ncol = length(columns)))
colnames(z_scores) <- columns

# Calculate z-scores for the selected columns
for (col in columns) {
  z_scores[[col]] <- (data[[col]] - mean(data[[col]], na.rm = TRUE)) / sd(data[[col]], na.rm = TRUE)
}

# Create a matrix to identify anomalies where z-scores are greater than 3 (or less than -3)
anomaly_matrix <- abs(z_scores) > 3

# Remove rows that have at least one anomaly
data_cleaned <- data[!apply(anomaly_matrix, 1, any), ]

# Creates a new Datetime column
data_cleaned$Datetime <- as.POSIXlt(paste(data_cleaned$Date, data_cleaned$Time), format = "%d/%m/%Y %H:%M:%S", tz = timezone)

# Converting Date column to actual Date data
data_cleaned$Date <- as.Date(data_cleaned$Date, format = "%d/%m/%Y")

# Check the date range to verify appropriate splitting
print(summary(data_cleaned$Date))

# Divide data into Test and Train 
train_data <- data_cleaned %>% filter(Date < as.Date("2009-01-01"))
test_data <- data_cleaned %>% filter(Date >= as.Date("2009-01-01"))

# Applying Standardization on Train Data #
train_scaled_data <- train_data 
num_cols <- names(train_data)[sapply(train_data, is.numeric)]
train_scaled_data[num_cols] <- scale(train_data[num_cols])

test_scaled_data <- test_data
test_num_cols <- names(test_data)[sapply(test_data, is.numeric)]
test_scaled_data[test_num_cols] <- scale(test_data[test_num_cols])






###### Part 2 ######


# Drop the time, date and timedate variables, not needed for PCA
# will be added back for training the HMM
pca_df <- train_scaled_data[, -c(1, 2, ncol(train_scaled_data))]

# Perform PCA
pca <- princomp(pca_df, scale = FALSE)

# View contributions to variance
summary_pca <- summary(pca)


######### THIS MAY BE REMOVED ########## 
# Variance explained by each principal component (standard deviations squared)
variance_explained <- pca$sdev^2
# Proportion of variance explained by each component
proportion_variance <- variance_explained / sum(variance_explained)
# Extract the proportion of variance explained by each of the first 3 components
var_exp <- proportion_variance[1:2]  # Proportion of variance explained
# Extract the loadings for the first 3 components
loadings <- pca$loadings[,1:2]
loadings
# Calculate the squared loadings for each feature and component
squared_loadings <- loadings^2
squared_loadings
# Repeat the variance explained for each feature (to match the number of rows in squared_loadings)
var_exp_repeated <- matrix(var_exp, nrow = nrow(squared_loadings), ncol = length(var_exp), byrow = TRUE)
var_exp_repeated
# Multiply squared loadings by the proportion of variance for each component
importance <- rowSums(squared_loadings * var_exp_repeated)
# Print the overall importance for each feature
importance
#########################################


# view the load of each feature for first 3 components
pca$loadings[, 1:3]

# scree plot
fviz_eig(pca, addlabels = TRUE)



fviz_cos2(pca, choice = "var", axes = 1:2)

# coloured biplot of component 1 and 2
fviz_pca_var(pca, col.var = "cos2",
             gradient.cols = c("black","red", "orange", "green","purple","blue"),
             repel = TRUE)

# another method to find biplot 
ggbiplot(pca,
         labels = pca$st ,
         circle = TRUE,
         varname.size = 4,
         varname.color = "red") 

# biplot between component 2 and 3
fviz_pca_var(pca, axes = c(2, 3), col.var = "black")
fviz_cos2(pca, choice = "var", axes = 2:3)
fviz_pca_var(pca, axes = c(2, 3), col.var = "cos2",
             gradient.cols = c("black", "red", "orange", "green", "purple", "blue"),
             repel = TRUE)

# loadings <- pca$rotation ## probably wont need, will clean up




### P3 ###
### USE THIS DF TO TRAIN THE HMM #### 
# The values in it may change later, but will still use the same name so
# code can be written for part 3
p3_data <- train_scaled_data
## use Global_intensity, Global_active_power, Sub_metering_3

# Specify response variables for HMM training (based on PCA results)
response_vars <- c("Global_intensity", "Global_active_power", "Sub_metering_3")


# columns for the weekday and week number
p3_data$Weekday <- wday(p3_data$Date, label = TRUE, abbr = TRUE)
p3_data$Week <- isoweek(p3_data$Date) # ISO week number


# Filter between 2 AM and 6 AM
mon_train_data <- filter(p3_data,
                           Weekday == "Mon" & 
                           format(Datetime, "%H:%M:%S") >= "02:00:00" & 
                           format(Datetime, "%H:%M:%S") <= "5:59:00")


# columns for the weekday and week number
test_scaled_data$Weekday <- wday(test_scaled_data$Date, label = TRUE, abbr = TRUE)
test_scaled_data$Week <- isoweek(test_scaled_data$Date) # ISO week number


# Filter between 2 AM and 6 AM
mon_test_data <- filter(test_scaled_data,
                         Weekday == "Mon" & 
                           format(Datetime, "%H:%M:%S") >= "02:00:00" & 
                           format(Datetime, "%H:%M:%S") <= "5:59:00")



weekly_counts <- group_by(mon_train_data, Date)
weekly_counts <- summarise(weekly_counts, Data_Points = n())

# Vector to be used for ntimes in the model
count_vector <- weekly_counts$Data_Points

count_vector

test_weekly_counts <- group_by(mon_test_data, Week)
test_weekly_counts <- summarise(test_weekly_counts, Data_Points = n())

# Vector to be used for ntimes in the model
test_count_vector <- test_weekly_counts$Data_Points

test_count_vector
# Subset training data with selected response variables
train_hmm_data <- mon_train_data[, response_vars]

# Prepare test data with the same response variables
test_hmm_data <- mon_test_data[, response_vars]

# Ensure no missing or invalid values in train and test datasets 
train_hmm_data <- train_hmm_data[complete.cases(train_hmm_data), ]
test_hmm_data <- test_hmm_data[complete.cases(test_hmm_data), ]

# Define the range of states to evaluate
range <- 4:10

states_range <- range[seq(1, length(range), by = 2)]

anyNA(mon_test_data)
# Train and evaluate HMM models
hmm_results <- train_evaluate_hmm(mon_train_data, mon_test_data, states_range, count_vector, test_count_vector)

# Extract results into a data frame for comparison
hmm_results_df <- do.call(rbind, lapply(hmm_results, function(x) {
  data.frame(
    num_states = x$num_states,
    train_log_likelihood = x$train_log_likelihood,
    test_log_likelihood = x$test_log_likelihood,
    bic = x$bic
  )
}))

# Visualize results: Log-Likelihood and BIC
ggplot(hmm_results_df, aes(x = num_states)) +
  geom_line(aes(y = train_log_likelihood, color = "Train Log-Likelihood")) +
  geom_line(aes(y = test_log_likelihood, color = "Test Log-Likelihood")) +
  geom_line(aes(y = bic, color = "BIC")) +
  labs(title = "HMM Evaluation Results", x = "Number of States", y = "Value") +
  scale_color_manual("", values = c("Train Log-Likelihood" = "blue", "Test Log-Likelihood" = "green", "BIC" = "red")) +
  theme_minimal()

# Select the best model based on BIC and log-likelihood
best_model <- hmm_results[[which.min(sapply(hmm_results, function(x) x$bic))]]

# Print best model information
print(paste("Best Model: ", best_model$num_states, "states"))
print(paste("Train Log-Likelihood: ", best_model$train_log_likelihood))
print(paste("Test Log-Likelihood: ", best_model$test_log_likelihood))
print(paste("BIC: ", best_model$bic))

