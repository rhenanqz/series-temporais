library(rugarch)
library(xts)
library(httr)
library(jsonlite)
library(lubridate)
library(plumber)
library(BatchGetSymbols)
library(dplyr)
library(tidyr)
library(tseries)
library(forecast)
library(quantmod)
library(fGarch)
library(rugarch)
library(rmgarch)
library(MTS)
library(ggplot2)
library(dlm)
library(timetk)
library(dynlm)
library(gridExtra)
library(grid)
library(patchwork)
library(xts)
library(zoo)
library(lubridate)
library(vars)
library(FinTS)
library(yfR)

# Function to get daily and houly data from cryptocurrencies

fetch_cryto_data <- function(fsym, 
                             tsym, 
                             timestamp_now = as.numeric(as.POSIXct("2023-08-31")),
                             frequency = "daily") {
  base_url <- ifelse(frequency == "daily",
                     "https://min-api.cryptocompare.com/data/v2/histoday",
                     "https://min-api.cryptocompare.com/data/v2/histohour")
  
  # Define the start and end times
  timestamp_now <- timestamp_now #as.numeric(Sys.time())
  timestamp_2018 <- as.numeric(as.POSIXct("2018-01-01"))
  
  params <- list(
    fsym = fsym,
    tsym = tsym,
    limit = 2000,  # Maximum allowed by the API for one call
    toTs = timestamp_now
  )
  
  # Fetch the data
  response <- GET(url = base_url, query = params)
  data <- fromJSON(content(response, as = "text"))$Data$Data
  
  # Continue fetching if more data is needed (due to the limit)
  while (data[1,1] > timestamp_2018) {
    params$toTs <- data[1,1]
    response <- GET(url = base_url, query = params)
    more_data <- fromJSON(content(response, as = "text"))$Data$Data
    data <- rbind(more_data, data)
  }
  
  data <- data %>%
    # Convert the timestamp to a datetime column
    mutate(datetime = as.POSIXct(time, origin = "1970-01-01", tz = "UTC")) %>%
    
    # Filter out rows before 2018-01-01
    filter(datetime >= as.POSIXct("2018-01-01", tz = "UTC")) %>%
    
    # Create a column for just the day
    mutate(day = as.Date(datetime)) %>%
    
    # Calculate hourly returns
    arrange(day, datetime) %>%
    mutate(returns = (close/lag(close) - 1) * 100) %>%
    
    # Select relevant columns
    dplyr::select(datetime, day, open, close, returns)
  return(data)
}

# Function extract data returns and realized vol
data_to_garch <- function(crypto = "BTC", comparative = "USD", timestamp_now = as.numeric(as.POSIXct("2023-08-31"))) {
  
  data_h <- fetch_cryto_data(fsym = crypto, 
                             tsym = comparative, 
                             timestamp_now = timestamp_now,
                             frequency = "hourly")
  data_d <- fetch_cryto_data(fsym = crypto, 
                             tsym = comparative, 
                             timestamp_now = timestamp_now,
                             frequency = "daily")
  
  daily <- data_d[!is.na(data_d$returns), ]
  daily_xts <- xts(daily$returns, order.by = as.Date(daily$day))
  
  hourly <- data_h[!is.na(data_h$returns), ]
  returns_h <- hourly$returns
  hourly$Date <- as.Date(hourly$datetime)  # Extract the date from the timestamp
  hourly_realized_var <- aggregate((returns/100)^2 ~ Date, data = hourly, sum)  # Sum of squared returns by day
  colnames(hourly_realized_var)[2] <- "realized_var"
  hourly_realized_vol <- sqrt(hourly_realized_var$realized_var)*100
  hourly_realized_vol_xts <- xts(hourly_realized_vol, order.by = as.Date(hourly_realized_var$Date))
  
  # Merge returns and realized volatility into a single xts object
  merged_data <- merge(daily_xts, hourly_realized_vol_xts, join="inner", suffixes = c("_returns", "_vol"))
  return(merged_data)
}


# Function to fit GARCH model and calculate metrics
fit_garch_model <- function(crypto, data_returns, data_vol, model_type, p, q, dist) {
  
  # This will hold the results of ugarchfit and ugarchroll
  fit <- NULL
  garchroll <- NULL
  
  # Try to fit the model and catch any errors or warnings
  tryCatch({
    if (model_type == "realGARCH") {
      rgarch.spec <- ugarchspec(mean.model = list(armaOrder= c(0,0),  include.mean = FALSE),
                                variance.model = list(model= 'realGARCH', garchOrder= c(p,q)), 
                                distribution.model = dist)
      fit <- ugarchfit(spec=rgarch.spec, data=data_returns, realizedVol=data_vol, solver = 'hybrid')
      garchroll <- ugarchroll(spec = rgarch.spec,
                              data= data_returns,
                              n.ahead = 1, 
                              forecast.length = 360,
                              refit.every = 5,
                              solver= 'hybrid',
                              realizedVol= data_vol,
                              refit.window = "recursive")
    } else {
      spec <- ugarchspec(variance.model = list(model = model_type, garchOrder = c(p, q)),
                         mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                         distribution.model = dist)
      fit <- ugarchfit(spec, data = data_returns, solver = 'hybrid')
      garchroll <- ugarchroll(spec,
                              data = data_returns, 
                              n.ahead = 1, 
                              forecast.length = 360,
                              refit.every = 5,
                              realizedVol= data_vol,
                              refit.window = "recursive",
                              solver= 'hybrid')
    }
  }, error = function(e) {
    fit <- NULL
    garchroll <- NULL
  }, warning = function(w) {
    fit <- NULL
    garchroll <- NULL
  })
  
  # Check if either 'fit' or 'garchroll' encountered an error
  if (is.null(fit) || is.null(garchroll)) {
    return(data.frame(Crypto = crypto, Model = model_type, p = p, q = q, Distribution = dist,
                      MSPE_in = NA, MAE_in = NA, MAPE_in = NA,
                      MSPE_out = NA, MAE_out = NA, MAPE_out = NA, BIC = NA, LogLikelihood = NA))
  } else {
    
    #in-sample
    realized_vol_fit = data_vol
    predicted_vol_fit <- xts(fit@fit$sigma, order.by = index(realized_vol_fit))
    
    
    #out-of-sample
    realized_vol_= tail(data_vol, 360)
    predicted_vol <- xts(garchroll@forecast$density$Sigma, order.by = index(realized_vol_))
    
    #metrics
    mspe_in <- mean((predicted_vol_fit - realized_vol_fit)^2, na.rm = TRUE)
    mae_in <- mean(abs(realized_vol_fit - predicted_vol_fit))
    mape_in <- mean(100 * abs((realized_vol_fit - predicted_vol_fit) / realized_vol_fit))
    
    mspe_out <- mean((predicted_vol - realized_vol_)^2, na.rm = TRUE)
    mae_out <- mean(abs(realized_vol_ - predicted_vol))
    mape_out <- mean(100 * abs((realized_vol_ - predicted_vol) / realized_vol_))
    bic <- infocriteria(fit)[2]
    llh <- fit@fit$LLH
    
    return(data.frame(Crypto = crypto, Model = model_type, p = p, q = q, Distribution = dist, 
                      MSPE_in = mspe_in, MAE_in = mae_in, MAPE_in = mape_in,
                      MSPE_out = mspe_out, MAE_out = mae_out, MAPE_out = mape_out, BIC = bic, LogLikelihood = llh))
  }
}


######### BITCOIN ############

merged_btc_data <- data_to_garch(crypto = "BTC", 
                                 comparative = "USD", 
                                 timestamp_now = as.numeric(as.POSIXct("2023-08-31")))
# Initialize an empty dataframe to store results

n = 1
#"sGARCH", "realGARCH", "gjrGARCH", 'fiGARCH', 'csGARCH',
# Loop through different model types, p, q, and dist
for (model in c("sGARCH", "realGARCH", "gjrGARCH",  'fiGARCH', 'csGARCH',  "eGARCH")) {
  final_results <- data.frame()
  for (p in 1:4) {
    for (q in 1:4) {
      for (dist in c('norm', 'std', 'sstd', 'ged', 'sged', 'jsu', 'ghyp', 'nig')) {
        # Fit the GARCH model and get the results as a dataframe
        cat("Realization:", n, "model_type:", model, ", p:", p, ", q:", q, ", dist:", dist, "\n")
        result_df <- fit_garch_model(crypto = "BTC",
                                     data_returns = merged_btc_data[, "X_returns"],
                                     data_vol = merged_btc_data[, "X_vol"],
                                     model_type = model,
                                     p = p,
                                     q = q,
                                     dist = dist)
        
        # Combine this result with the final_results dataframe
        final_results <- rbind(final_results, result_df)
        n = n + 1
      }
    }
  }
  filename <- paste0('final_garch_results_btc2_', model, '.csv')
  write.csv(final_results, filename, row.names = FALSE)
}




######### ETHEREUM ############

merged_eth_data <- data_to_garch(crypto = "ETH", 
                                 comparative = "USD", 
                                 timestamp_now = as.numeric(as.POSIXct("2023-08-31")))
# Initialize an empty dataframe to store results

n = 1
#"sGARCH", "realGARCH", "gjrGARCH", 'fiGARCH', 'csGARCH',
# Loop through different model types, p, q, and dist
for (model in c("sGARCH", "realGARCH", "gjrGARCH",  'fiGARCH', 'csGARCH',  "eGARCH")) {
  final_results <- data.frame()
  for (p in 1:4) {
    for (q in 1:4) {
      for (dist in c('norm', 'std', 'sstd', 'ged', 'sged', 'jsu', 'ghyp', 'nig')) {
        # Fit the GARCH model and get the results as a dataframe
        cat("Realization:", n, "model_type:", model, ", p:", p, ", q:", q, ", dist:", dist, "\n")
        result_df <- fit_garch_model(crypto = "ETH",
                                     data_returns = merged_eth_data[, "X_returns"],
                                     data_vol = merged_eth_data[, "X_vol"],
                                     model_type = model,
                                     p = p,
                                     q = q,
                                     dist = dist)
        
        # Combine this result with the final_results dataframe
        final_results <- rbind(final_results, result_df)
        n = n + 1
      }
    }
  }
  filename <- paste0('final_garch_results_eth2_', model, '.csv')
  write.csv(final_results, filename, row.names = FALSE)
}
######### XRP ############

merged_xrp_data <- data_to_garch(crypto = "XRP", 
                                 comparative = "USD", 
                                 timestamp_now = as.numeric(as.POSIXct("2023-08-31")))
# Initialize an empty dataframe to store results

n = 1
#"sGARCH", "realGARCH", "gjrGARCH", 'fiGARCH', 'csGARCH',
# Loop through different model types, p, q, and dist
for (model in c("sGARCH", "realGARCH", "gjrGARCH",  'fiGARCH', 'csGARCH',  "eGARCH")) {
  final_results <- data.frame()
  for (p in 1:4) {
    for (q in 1:4) {
      for (dist in c('norm', 'std', 'sstd', 'ged', 'sged', 'jsu', 'ghyp', 'nig')) {
        # Fit the GARCH model and get the results as a dataframe
        cat("Realization:", n, "model_type:", model, ", p:", p, ", q:", q, ", dist:", dist, "\n")
        result_df <- fit_garch_model(crypto = "XRP",
                                     data_returns = merged_xrp_data[, "X_returns"],
                                     data_vol = merged_xrp_data[, "X_vol"],
                                     model_type = model,
                                     p = p,
                                     q = q,
                                     dist = dist)
        
        # Combine this result with the final_results dataframe
        final_results <- rbind(final_results, result_df)
        n = n + 1
      }
    }
  }
  filename <- paste0('final_garch_results_xrp2_', model, '.csv')
  write.csv(final_results, filename, row.names = FALSE)
}

######### BNB ############

merged_bnb_data <- data_to_garch(crypto = "BNB", 
                                 comparative = "USD", 
                                 timestamp_now = as.numeric(as.POSIXct("2023-08-31")))
# Find the maximum finite value in the data frame
max_finite_value <- max(merged_bnb_data[merged_bnb_data$X_vol != Inf], na.rm = TRUE)

# Replace Inf with the maximum finite value
merged_bnb_data$X_vol[merged_bnb_data$X_vol == Inf] <- max_finite_value
# Initialize an empty dataframe to store results

n = 1
#"sGARCH", "realGARCH", "gjrGARCH", 'fiGARCH', 'csGARCH',
# Loop through different model types, p, q, and dist
for (model in c("sGARCH", "realGARCH", "gjrGARCH",  'fiGARCH', 'csGARCH',  "eGARCH")) {
  final_results <- data.frame()
  for (p in 1:4) {
    for (q in 1:4) {
      for (dist in c('norm', 'std', 'sstd', 'ged', 'sged', 'jsu', 'ghyp', 'nig')) {
        # Fit the GARCH model and get the results as a dataframe
        cat("Realization:", n, "model_type:", model, ", p:", p, ", q:", q, ", dist:", dist, "\n")
        result_df <- fit_garch_model(crypto = "BNB",
                                     data_returns = merged_bnb_data[, "X_returns"],
                                     data_vol = merged_bnb_data[, "X_vol"],
                                     model_type = model,
                                     p = p,
                                     q = q,
                                     dist = dist)
        
        # Combine this result with the final_results dataframe
        final_results <- rbind(final_results, result_df)
        n = n + 1
      }
    }
  }
  filename <- paste0('final_garch_results_bnb2_', model, '.csv')
  write.csv(final_results, filename, row.names = FALSE)
}

######### CARDANO ############

merged_ada_data <- data_to_garch(crypto = "ADA", 
                                 comparative = "USD", 
                                 timestamp_now = as.numeric(as.POSIXct("2023-08-31")))
# Initialize an empty dataframe to store results

n = 1
#"sGARCH", "realGARCH", "gjrGARCH", 'fiGARCH', 'csGARCH',
# Loop through different model types, p, q, and dist
for (model in c("sGARCH", "realGARCH", "gjrGARCH",  'fiGARCH', 'csGARCH',  "eGARCH")) {
  final_results <- data.frame()
  for (p in 1:4) {
    for (q in 1:4) {
      for (dist in c('norm', 'std', 'sstd', 'ged', 'sged', 'jsu', 'ghyp', 'nig')) {
        # Fit the GARCH model and get the results as a dataframe
        cat("Realization:", n, "model_type:", model, ", p:", p, ", q:", q, ", dist:", dist, "\n")
        result_df <- fit_garch_model(crypto = "ADA",
                                     data_returns = merged_ada_data[, "X_returns"],
                                     data_vol = merged_ada_data[, "X_vol"],
                                     model_type = model,
                                     p = p,
                                     q = q,
                                     dist = dist)
        
        # Combine this result with the final_results dataframe
        final_results <- rbind(final_results, result_df)
        n = n + 1
      }
    }
  }
  filename <- paste0('final_garch_results_ada2_', model, '.csv')
  write.csv(final_results, filename, row.names = FALSE)
}

######### DOGE ############

merged_doge_data <- data_to_garch(crypto = "DOGE", 
                                  comparative = "USD", 
                                  timestamp_now = as.numeric(as.POSIXct("2023-08-31")))
# Initialize an empty dataframe to store results

n = 1
#"sGARCH", "realGARCH", "gjrGARCH", 'fiGARCH', 'csGARCH',
# Loop through different model types, p, q, and dist
for (model in c("sGARCH", "realGARCH", "gjrGARCH",  'fiGARCH', 'csGARCH',  "eGARCH")) {
  final_results <- data.frame()
  for (p in 1:4) {
    for (q in 1:4) {
      for (dist in c('norm', 'std', 'sstd', 'ged', 'sged', 'jsu', 'ghyp', 'nig')) {
        # Fit the GARCH model and get the results as a dataframe
        cat("Realization:", n, "model_type:", model, ", p:", p, ", q:", q, ", dist:", dist, "\n")
        result_df <- fit_garch_model(crypto = "DOGE",
                                     data_returns = merged_doge_data[, "X_returns"],
                                     data_vol = merged_doge_data[, "X_vol"],
                                     model_type = model,
                                     p = p,
                                     q = q,
                                     dist = dist)
        
        # Combine this result with the final_results dataframe
        final_results <- rbind(final_results, result_df)
        n = n + 1
      }
    }
  }
  filename <- paste0('final_garch_results_doge2_', model, '.csv')
  write.csv(final_results, filename, row.names = FALSE)
}



# Read the CSV file into a dataframe
read_results <- read.csv("final_garch_results_btc.csv")

# View the first few rows of the dataframe
head(read_results)
