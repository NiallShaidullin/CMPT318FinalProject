library(dplyr)
library(corrplot)
library(ggplot2)
library(lubridate)
library(zoo)
library(imputeTS)
library("depmixS4")
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


