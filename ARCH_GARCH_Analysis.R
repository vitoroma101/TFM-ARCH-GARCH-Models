# Load necessary libraries
library(dplyr)
library(quantmod)
library(xts)
library(openxlsx)
library(moments)
library(tseries)
library(gt)
library(TSA)
library(forecast)
library(dynlm)
library(ggplot2)
library(gridExtra)
library(knitr)
library(rugarch)

# Load data from GitHub
# Since the file is large, it is not displayed directly on GitHub. You should first download the file from the provided URL and then replace the path with its location on your PC.

data <- read.xlsx(
  "https://github.com/vitoroma101/TFM-ARCH-GARCH-Models/blob/main/BVG-acciones.xlsx",
  startRow = 3,
  colNames = TRUE
)

# Convert the date column from numeric format to Date format
data$FECHA.NEGOCIACIÓN <- as.Date(data$FECHA.NEGOCIACIÓN, origin = "1899-12-30")

# Filter data for the company "CORPORACIÓN FAVORITA C.A."
# and select the last closing price per day
newdata <- data %>%
  filter(EMISOR == "CORPORACION FAVORITA C.A.") %>%
  group_by(FECHA.NEGOCIACIÓN) %>%
  slice_tail(n = 1) %>%
  ungroup()

# Convert filtered data to xts format for time series analysis
data_xts <- xts(newdata$PRECIO, order.by = newdata$FECHA.NEGOCIACIÓN)

# Plot Corporación Favorita stock price
plot(data_xts, main = "")

# Calculate daily log returns
log_returns_daily <- diff(log(data_xts))
log_returns_daily <- na.omit(log_returns_daily)

# Summary of daily log returns
summary(log_returns_daily)

# Convert daily log returns to a data frame for ggplot
log_returns_df <- data.frame(date = index(log_returns_daily),
                             log_returns = coredata(log_returns_daily))

# Histogram of daily log returns
p <- ggplot(log_returns_df, aes(x = log_returns)) +
  geom_histogram(color = "black", fill = "white") +
  labs(x = "Daily Returns", y = "Frequency") +
  theme_minimal()
p

# Create plots for daily, squared, and absolute log returns
par(mfrow=c(3,1))
plot(log_returns_daily,
     main = "Daily Log Returns for CORPORACIÓN FAVORITA",
     col = "blue")
plot(log_returns_daily^2,
     main = "Squared Log Returns for CORPORACIÓN FAVORITA",
     col = "red")
plot(abs(log_returns_daily),
     main = "Absolute Log Returns for CORPORACIÓN FAVORITA",
     col = "green")

# Create ACF and PACF plots
par(mfrow = c(3, 1))
acf(log_returns_daily,
    main = "ACF of Log Returns for CORPORACIÓN FAVORITA",
    col = "blue",
    lag.max = 20)
acf(log_returns_daily^2,
    main = "ACF of Squared Log Returns for CORPORACIÓN FAVORITA",
    col = "red",
    lag.max = 20)
acf(abs(log_returns_daily),
    main = "ACF of Absolute Log Returns for CORPORACIÓN FAVORITA",
    col = "green",
    lag.max = 20)

# Ljung-Box tests for different lags
ljung_box_results <- lapply(c(1, 5, 10), function(lag) {
  list(
    returns = Box.test(log_returns_daily, lag = lag, type = "Ljung-Box"),
    squared = Box.test(log_returns_daily^2, lag = lag, type = "Ljung-Box"),
    abs = Box.test(abs(log_returns_daily), lag = lag, type = "Ljung-Box")
  )
})

# Lagrange Multiplier (LM) tests for ARCH effects at different lags
arch_lm_results <- lapply(c(1, 5, 10), function(lag) {
  ArchTest(log_returns_daily, lags = lag)
})

# Fit ARIMA model
modelo_arima <- auto.arima(log_returns_daily)
modelo_arima

# Fit ARIMA(0,0,3) model with the ma2 coefficient fixed to 0
modelo_arima_sin_ma2 <- Arima(log_returns_daily, order = c(0, 0, 3),
                              fixed = c(NA, 0, NA, NA))
modelo_arima_sin_ma2

# Fit ARIMA(0,0,3) model without the mean
modelo_arima_sin_media <- Arima(log_returns_daily, order = c(0, 0, 3),
                                fixed = c(NA, 0, NA), include.mean = FALSE)
modelo_arima_sin_media

# Plot and check residuals
autoplot(modelo_arima_sin_media)
checkresiduals(modelo_arima_sin_media)

# Calculate errors and squared errors
Errors <- resid(modelo_arima_sin_media)
Squared_Errors <- resid(modelo_arima_sin_media)^2
chartSeries(Squared_Errors, name = "Squared Errors")

# Fit a regression model with squared errors
Regresion1 <- dynlm(Squared_Errors ~ L(Squared_Errors))
summary(Regresion1)

# Generate ACF and PACF for errors and squared errors
plot_acf_errors <- autoplot(acf(Errors, lag.max = 20, plot = FALSE)) +
  labs(title = "ACF of Errors", x = "Lags", y = "Autocorrelation")
plot_pacf_errors <- autoplot(pacf(Errors, lag.max = 20, plot = FALSE)) +
  labs(title = "PACF of Errors", x = "Lags", y = "Partial Autocorrelation")
plot_acf_squared_errors <- autoplot(acf(Squared_Errors, lag.max = 20, plot = FALSE)) +
  labs(title = "ACF of Squared Errors", x = "Lags", y = "Autocorrelation")
plot_pacf_squared_errors <- autoplot(pacf(Squared_Errors, lag.max = 20, plot = FALSE)) +
  labs(title = "PACF of Squared Errors", x = "Lags", y = "Partial Autocorrelation")

# Combine all plots into one graphic
grid.arrange(plot_acf_errors, plot_pacf_errors,
             plot_acf_squared_errors, plot_pacf_squared_errors, ncol = 2)

# Specification of ARCH and GARCH models
# Extract residuals from ARIMA model without mean
residuals_arima <- residuals(modelo_arima_sin_media)

# Fit ARCH(1) model
spec_arch1 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         variance.model = list(model = "sGARCH", garchOrder = c(1, 0)))
fit_arch1 <- ugarchfit(spec = spec_arch1, data = residuals_arima)
arch1_ic <- infocriteria(fit_arch1)
arch1_loglik <- likelihood(fit_arch1)
print(arch1_ic)
print(arch1_loglik)

# Fit ARCH(2) model
spec_arch2 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         variance.model = list(model = "sGARCH", garchOrder = c(2, 0)))
fit_arch2 <- ugarchfit(spec = spec_arch2, data = residuals_arima)
arch2_ic <- infocriteria(fit_arch2)
arch2_loglik <- likelihood(fit_arch2)
print(arch2_ic)
print(arch2_loglik)

# Extract standardized residuals from ARCH(2)
standardized_residuals_arch2 <- residuals(fit_arch2, standardize = TRUE)

# Ljung-Box test for ARCH(2) standardized residuals
ljung_box_arch2 <- Box.test(standardized_residuals_arch2, lag = 12, type = "Ljung-Box")

# ARCH LM test for ARCH(2)
arch_lm_arch2 <- ArchTest(standardized_residuals_arch2, lags = 12)

# Jarque-Bera test for ARCH(2)
jb_test_arch2 <- jarque.bera.test(standardized_residuals_arch2)

# Fit ARCH(3) model
spec_arch3 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         variance.model = list(model = "sGARCH", garchOrder = c(3, 0)))
fit_arch3 <- ugarchfit(spec = spec_arch3, data = residuals_arima)
arch3_ic <- infocriteria(fit_arch3)
arch3_loglik <- likelihood(fit_arch3)
print(arch3_ic)
print(arch3_loglik)

# Extract standardized residuals from ARCH(3)
standardized_residuals_arch3 <- residuals(fit_arch3, standardize = TRUE)

# Ljung-Box test for ARCH(3) standardized residuals
ljung_box_arch3 <- Box.test(standardized_residuals_arch3, lag = 12, type = "Ljung-Box")

# ARCH LM test for ARCH(3)
arch_lm_arch3 <- ArchTest(standardized_residuals_arch3, lags = 12)

# Jarque-Bera test for ARCH(3)
jb_test_arch3 <- jarque.bera.test(standardized_residuals_arch3)

# Fit ARCH(4) model
spec_arch4 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         variance.model = list(model = "sGARCH", garchOrder = c(4, 0)))
fit_arch4 <- ugarchfit(spec = spec_arch4, data = residuals_arima)
arch4_ic <- infocriteria(fit_arch4)
arch4_loglik <- likelihood(fit_arch4)
print(arch4_ic)
print(arch4_loglik)

# Fit ARCH(5) model
spec_arch5 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                         variance.model = list(model = "sGARCH", garchOrder = c(5, 0)))
fit_arch5 <- ugarchfit(spec = spec_arch5, data = residuals_arima)
arch5_ic <- infocriteria(fit_arch5)
arch5_loglik <- likelihood(fit_arch5)
print(arch5_ic)
print(arch5_loglik)

# Fit GARCH(1,1) model
spec_garch11 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                           variance.model = list(model = "sGARCH", garchOrder = c(1, 1)))
fit_garch11 <- ugarchfit(spec = spec_garch11, data = residuals_arima, solver = "hybrid")
garch11_ic <- infocriteria(fit_garch11)
garch11_loglik <- likelihood(fit_garch11)
print(garch11_ic)
print(garch11_loglik)

# Fit GARCH(2,1) model
spec_garch21 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                           variance.model = list(model = "sGARCH", garchOrder = c(2, 1)))
fit_garch21 <- ugarchfit(spec = spec_garch21, data = residuals_arima)
garch21_ic <- infocriteria(fit_garch21)
garch21_loglik <- likelihood(fit_garch21)
print(garch21_ic)
print(garch21_loglik)

# Extract standardized residuals from GARCH(2,1)
standardized_residuals_garch21 <- residuals(fit_garch21, standardize = TRUE)

# Ljung-Box test for GARCH(2,1) standardized residuals
ljung_box_garch21 <- Box.test(standardized_residuals_garch21, lag = 12, type = "Ljung-Box")

# ARCH LM test for GARCH(2,1)
arch_lm_garch21 <- ArchTest(standardized_residuals_garch21, lags = 12)

# Jarque-Bera test for GARCH(2,1)
jb_test_garch21 <- jarque.bera.test(standardized_residuals_garch21)

# Fit GARCH(1,2) model
spec_garch12 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                           variance.model = list(model = "sGARCH", garchOrder = c(1, 2)))
fit_garch12 <- ugarchfit(spec = spec_garch12, data = residuals_arima)
garch12_ic <- infocriteria(fit_garch12)
garch12_loglik <- likelihood(fit_garch12)
print(garch12_ic)
print(garch12_loglik)

# Fit GARCH(2,2) model
spec_garch22 <- ugarchspec(mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                           variance.model = list(model = "sGARCH", garchOrder = c(2, 2)))
fit_garch22 <- ugarchfit(spec = spec_garch22, data = residuals_arima)
garch22_ic <- infocriteria(fit_garch22)
garch22_loglik <- likelihood(fit_garch22)
print(garch22_ic)
print(garch22_loglik)

# Extract conditional volatilities from the models
vol_arch2 <- sigma(fit_arch2)
vol_arch3 <- sigma(fit_arch3)
vol_garch21 <- sigma(fit_garch21)

# Create dataframes for ggplot2
data_arch2 <- data.frame(Date = index(log_returns_daily),
                         Volatility = vol_arch2)
data_arch3 <- data.frame(Date = index(log_returns_daily),
                         Volatility = vol_arch3)
data_garch21 <- data.frame(Date = index(log_returns_daily),
                           Volatility = vol_garch21)
data_log_returns <- data.frame(Date = index(log_returns_daily),
                               LogReturns = as.numeric(log_returns_daily))

# Plot daily returns
plot_returns <- ggplot(data_log_returns, aes(x = Date, y = LogReturns)) +
  geom_line() +
  labs(title = "Daily Returns of La Favorita") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Plot conditional volatility ARCH(2)
plot_arch2 <- ggplot(data_arch2, aes(x = Date, y = Volatility)) +
  geom_line() +
  labs(title = "Conditional Volatility of La Favorita (ARCH(2))") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Plot conditional volatility ARCH(3)
plot_arch3 <- ggplot(data_arch3, aes(x = Date, y = Volatility)) +
  geom_line() +
  labs(title = "Conditional Volatility of La Favorita (ARCH(3))") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Plot conditional volatility GARCH(2,1)
plot_garch21 <- ggplot(data_garch21, aes(x = Date, y = Volatility)) +
  geom_line() +
  labs(title = "Conditional Volatility of La Favorita (GARCH(2,1))") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())

# Display all plots
grid.arrange(plot_returns, plot_arch2, plot_arch3, plot_garch21, ncol = 1)

# FINAL MODEL
# Extract standardized residuals and squared residuals from the ARCH(2) model
standardized_residuals <- residuals(fit_arch2, standardize = TRUE)
squared_residuals <- standardized_residuals^2

# Plot ACF of squared residuals
par(mfrow = c(2, 1))  # Set up graphical window for two plots
acf(squared_residuals, main = "Squared Residuals La Favorita ACF", lag.max = 20)

# Q-Q plot of standardized residuals
qqnorm(standardized_residuals, main = "", ylab = "Standardized Residuals of La Favorita", xlab = "Standard Normal Quantiles")
qqline(standardized_residuals, col = "red")

# Perform forecast
ug_forecast <- ugarchforecast(fit_arch2, n.ahead = 7)
ug_forecast
