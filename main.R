# Data Overview

## Check and install required packages if not already installed
list.of.packages <- c("FinTS", "ggplot2", "dplyr", "lubridate", "rjags", "reshape2",
                      "R2jags", "coda", "stargazer","gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load necessary libraries
library("FinTS")
library("ggplot2")
library("gridExtra")
library("dplyr")
library("lubridate")
library("rjags")
library("reshape2")
library("R2jags")
library("coda")
library("stargazer")

# Load the data
data("q.gnp4791", package = "FinTS")
# Convert the data to a dataframe
q.gnp.df <- data.frame("Date"=as.yearqtr(rownames(as.data.frame(q.gnp4791)), 
                                         format="%Y (%q)"), "GNP_US"=q.gnp4791)
q.gnp.df$Date <- as.Date(yearmon(q.gnp.df$Date))

# Plotting the Quarterly Growth Rate of U.S. Real GNP
png(file="out/fig_01.png", width=8, height=4, units="in", res=600)
q.gnp.df %>%
  ggplot() +
  geom_line(aes(x = Date, y = GNP_US), alpha = 0.8) +
  scale_color_manual(values = c("#004C99")) +
  geom_hline(yintercept = mean(q.gnp.df$GNP_US), linetype = "dashed", color = "red") +
  geom_hline(yintercept = quantile(q.gnp.df$GNP_US, 0.95), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = quantile(q.gnp.df$GNP_US, 0.05), linetype = "dashed", color = "blue") +
  ggtitle("U.S. Real GNP Over Time") +
  ylab("Price") + xlab("Date") +
  scale_y_continuous(n.breaks = 12) +
  theme_bw() +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) #transparent legend panel
dev.off()

# Classical ARIMA Approach for Time Series Modeling

## Selection of Optimal Model

# Set graphical parameters for the plot
png(file="out/fig_02.png", width=8, height=4, units="in", res=600)
par(mfrow=c(1,2), mar=c(3, 3, 0.1, 0.1), mgp=c(1.5, 0.5, 0), cex=0.7)
# Plot the autocorrelation function (ACF) of the time series 'q.gnp4791'
acf(q.gnp4791, main="\n")
pacf(q.gnp4791, main = "\n")
dev.off()

# Define a function to fit an ARIMA model with specified order 'p' 
# to the time series data 'q.gnp4791'
fit_AR_model <- function(p) {
  # Use the arima function with the specified order (p, 0, 0) 
  # and the maximum likelihood (ML) method
  model <- arima(q.gnp4791, order = c(p, 0, 0), method = "ML")
  # Return the fitted ARIMA model
  return(model)
}

# Generate a sequence of p values from 1 to 20
p_values <- 1:20  

# Apply the fit_AR_model function to each value of 'p' in the sequence, 
# creating a list of ARIMA models
models <- lapply(p_values, fit_AR_model)

# Extract the residuals from each fitted model and store them in a list
models_residuals <- lapply(models, residuals)

# Calculate AIC (Akaike Information Criterion) values for each fitted ARIMA model
AIC_values <- sapply(models, AIC)

# Calculate BIC (Bayesian Information Criterion) values for each fitted ARIMA model
BIC_values <- sapply(models, BIC)

# Set graphical parameters for the plot
png(file="out/fig_03.png", width=8, height=4, units="in", res=600)
par(mfrow=c(1,2), mar=c(3, 3, 2, 0.5), mgp=c(1.5, 0.5, 0), cex=0.7)
# Create a plot of AIC and BIC scores against the different values of 'p'
plot(p_values, AIC_values, xlab="p", ylab="AIC scores", type="b", pch=19)
plot(p_values, BIC_values, xlab="p", ylab="BIC scores", type="b", pch=19)
dev.off()


## Optimal Model

# Find the value of 'p' that minimizes the AIC values
optimal_p_AIC <- p_values[which.min(AIC_values)]
# Find the value of 'p' that minimizes the BIC values
optimal_p_BIC <- p_values[which.min(BIC_values)]

# Fit the optimal ARIMA model using the identified 'optimal_p_AIC' value
AIC_model <- arima(q.gnp4791, order = c(optimal_p_AIC, 0, 0), method = "ML")
# Fit the optimal ARIMA model using the identified 'optimal_p_BIC' value
BIC_model <- arima(q.gnp4791, order = c(optimal_p_BIC, 0, 0), method = "ML")

stargazer(AIC_model, BIC_model, column.labels = c("AR(p=3)","AR(p=1)"), header = FALSE, 
          dep.var.caption = "", dep.var.labels.include = FALSE, single.row = TRUE,
          type = "html", out = "out/table_01.html",
          title="Coefficient Estimates and Model Comparison for ARIMA Models.")

# Plot the autocorrelation function (ACF) of the residuals of the optimal models
png(file="out/fig_04.png", width=8, height=4, units="in", res=600)
par(mfrow=c(1,2), mar=c(3, 3, 2.8, 0.1), mgp=c(1.5, 0.5, 0), cex=0.7)
acf(AIC_model$residuals, main="AR(p=3)")
acf(BIC_model$residuals, main="AR(p=1)")
dev.off()

## Forecast 

# Generate fitted values 
fitted_values_AIC <- q.gnp.df$GNP_US - AIC_model$residuals

# Generate fitted values 
fitted_values_BIC <- q.gnp.df$GNP_US - BIC_model$residuals

# Generate forecasts with confidence intervals for the next 2 time points
forecast_values_AIC <- predict(AIC_model, n.ahead = 2, interval = "prediction", level = 0.95)

# Generate forecasts with confidence intervals for the next 2 time points
forecast_values_BIC <- predict(BIC_model, n.ahead = 2, interval = "prediction", level = 0.95)

# Set graphical parameters for the plot
png(file="out/fig_05.png", width=8, height=4, units="in", res=600)
par(mar=c(3, 3, 0.1, 0.1), mgp=c(1.5, 0.5, 0), cex=0.7)
# Plot the original time series
plot(q.gnp4791, type = "l", col = "blue", ylab = "Growth Rate")
# Plot fitted values
lines(fitted_values_AIC, col = "red")
# Plot forecasted values
lines(forecast_values_AIC$pred, col = "red", lty = 1)
# Plot confidence intervals
lines(forecast_values_AIC$pred + forecast_values_AIC$se, col = "red", lty = 3)  # Upper 
lines(forecast_values_AIC$pred - forecast_values_AIC$se, col = "red", lty = 3)  # Lower 
# Plot fitted values
lines(fitted_values_BIC, col = "darkgreen")
# Plot forecasted values
lines(forecast_values_BIC$pred, col = "darkgreen", lty = 1)
# Plot confidence intervals
lines(forecast_values_BIC$pred + forecast_values_BIC$se, col = "darkgreen", lty = 3)  # Upper 
lines(forecast_values_BIC$pred - forecast_values_BIC$se, col = "darkgreen", lty = 3)  # Lower 
# Add legend
legend("topright", ncol = 2,  col = c("red", "red", "darkgreen", "darkgreen", "blue"), 
       lty = c(1, 3, 1, 3, 1), bg="white",
       legend = c("Forecast (p=3)","95% Confidence Intervals (p=3)", "Forecast (p=1)", 
                  "95% Confidence Intervals (p=1)"))
dev.off()

# Bayesian Approach for Time Series Modeling

# Extracting the 'GNP_US' column from the 'q.gnp.df' data frame
gnp <- q.gnp.df$GNP_US

# Number of observations
N <- dim(q.gnp.df)[1]

# Defining the location to save the JAGS AR(1) model script
model.loc.ar1 <- "AR1_gnp.txt"
# Creating the JAGS script for the AR(1) model
jagsscript <- cat("
                  model {  
                     # priors on parameters
                     u ~ dnorm(0, 0.01); 
                     inv.q ~ dgamma(0.001,0.001); 
                     q <- 1/inv.q; 
                     b ~ dunif(-1,1);
                     X0 ~ dnorm(0, inv.q * (1 - b * b));
                     
                     # likelihood
                     X[1] ~ dnorm(b * X0 + u, inv.q);
                     
                     for(t in 2:N) {
                        EX[t] <- X[t-1] * b + u
                        X[t] ~ dnorm(b * X[t-1] + u, inv.q);
                     }
                  }  
                  ", file = model.loc.ar1)

# Defining the location to save the JAGS AR(1) model script
model.loc.ar3 <- ("AR3_gnp.txt")
# Creating the JAGS script for the AR(3) model
jagsscript_ar3 <- cat("
                  model {
                  
                     # priors on parameters
                     u ~ dnorm(0, 0.01); 
                     inv.q ~ dgamma(0.001,0.001); 
                     q <- 1/inv.q;
                     
                     # uniform priors for the model coefficients
                     b1 ~ dunif(-1,1);
                     b2 ~ dunif(-1,1);
                     b3 ~ dunif(-1,1);
                     
                     # normal prior for the coefficients
                     X0 ~ dnorm(0, inv.q * (1 - (b1 * b1 + b2 * b2 + b3 * b3)));
                     X1 ~ dnorm(0, inv.q * (1 - (b1 * b1 + b2 * b2 + b3 * b3)));
                     X2 ~ dnorm(0, inv.q * (1 - (b1 * b1 + b2 * b2 + b3 * b3)));
                     
                     # likelihood
                     X[1] ~ dnorm(b1 * X0 + b2 * X1 + b3 * X2 + u, inv.q);
                     X[2] ~ dnorm(b1 * X[1] + b2 * X0 + b3 * X1 + u, inv.q);
                     X[3] ~ dnorm(b1 * X[2] + b2 * X[1] + b3 * X0 + u, inv.q);
                     
                     for(t in 4:N) {
                        EX[t] <- X[t-1] * b1 + X[t-2] * b2 + X[t-3] * b3 + u
                        X[t] ~ dnorm(X[t-1] * b1 + X[t-2] * b2 + X[t-3] * b3 + u, inv.q);
                     }
                  }  
                  ", file = model.loc.ar3)

# Creating a list of data to pass to JAGS
jags.data <- list(X = gnp, N = N)

# Specifying parameters for AR(1) to be saved
jags.params.ar1 <- c("q", "u", "b", "EX", "X")

# Running the JAGS for AR(1) model
mod_ar1_intercept <- R2jags::jags(jags.data, parameters.to.save = jags.params.ar1, 
                                  model.file = model.loc.ar1, n.chains = 3, n.burnin = 5000, n.thin = 5, 
                                  n.iter = 10000, DIC = TRUE)

# Specifying parameters for AR(3) to be saved
jags.params.ar3 <- c("q", "u", "b1", "b2", "b3", "EX", "X")

# Running the JAGS for AR(3) model
mod_ar3_intercept <- R2jags::jags(jags.data, parameters.to.save = jags.params.ar3, 
                                  model.file = model.loc.ar3, n.chains = 3, n.burnin = 5000, n.thin = 5, 
                                  n.iter = 10000, DIC = TRUE)

## In-Sample Results

# Getting quantiles and means of the predictions of AR(1) model
EX_ar1 <- mod_ar1_intercept$BUGSoutput$sims.list$EX
X_ar1 <- mod_ar1_intercept$BUGSoutput$sims.list$X

# Creating a data frame for plotting AR(1) model
plot_df_ar1 <- as.data.frame(cbind(apply(EX_ar1, 2, quantile, 0.025), 
                                   apply(EX_ar1, 2, mean), 
                                   apply(EX_ar1, 2, quantile, 0.975)))
colnames(plot_df_ar1) <- c("lower", "mean_pred", "upper")
plot_df_ar1$fit <- apply(X_ar1, 2 , mean)[-1]

# Adding date and observed values to the data frame for AR(1) model
plot_df_ar1$Date <- q.gnp.df$Date[-1]
plot_df_ar1$obs <- q.gnp.df$GNP_US[-1]

# Creating the plot for AR(1) model
plot_fit_ar1 <- plot_df_ar1 %>% 
  ggplot() +
  geom_point(aes(x = Date, y = obs), alpha = 0.8, size = 0.5, col = "black") +
  geom_line(aes(x = Date, y = fit), alpha = 1, col = "black", linewidth = 0.3) +
  geom_line(aes(x = Date, y = mean_pred), alpha = 0.8, linewidth = 0.3,  col = "red") +
  ylab("Growth Rate") +
  scale_y_continuous(n.breaks = 12) +
  theme_bw() +
  theme(text = element_text(size = 8),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) +  #transparent legend panel
  ggtitle("AR(1) model")

# Getting quantiles and means of the predictions of AR(3) model
EX_ar3 <- mod_ar3_intercept$BUGSoutput$sims.list$EX
X_ar3 <- mod_ar3_intercept$BUGSoutput$sims.list$X

# Creating a data frame for plotting AR(3) model
plot_df_ar3 <- as.data.frame(cbind(apply(EX_ar3, 2, quantile, 0.025),
                                   apply(EX_ar3, 2, mean), 
                                   apply(EX_ar3, 2, quantile, 0.975)))
colnames(plot_df_ar3) <- c("lower", "mean_pred", "upper")
plot_df_ar3$fit <- apply(X_ar3, 2, mean)[4:N]

# Adding date and observed values to the data frame for AR(3) model
plot_df_ar3$Date <- q.gnp.df$Date[4:N]
plot_df_ar3$obs <- q.gnp.df$GNP_US[4:N]

plot_fit_ar3 <- plot_df_ar3 %>% 
  ggplot() +
  geom_point(aes(x = Date, y = obs), alpha = 0.8, size = 0.5,col = "black") +
  geom_line(aes(x = Date, y = fit), alpha = 1, col = "black", linewidth = 0.3) +
  geom_line(aes(x = Date, y = mean_pred), alpha = 0.8, linewidth = 0.3,  col = "red") +
  ylab("Growth Rate") +
  scale_y_continuous(n.breaks = 12) +
  theme_bw() +
  theme(text = element_text(size = 8),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent')) +  #transparent legend panel
  ggtitle("AR(3) model")

png(file="out/fig_06.png", width=8, height=4, units="in", res=600)
grid.arrange(plot_fit_ar1, plot_fit_ar3, ncol=2, widths=c(1,1))
dev.off()

# Extracting posterior parameters from the JAGS output
post.params.ar1 <- mod_ar1_intercept$BUGSoutput$sims.list

# Extracting posterior parameters from the JAGS output
post.params.ar3 <- mod_ar3_intercept$BUGSoutput$sims.list

# Creating a histogram for each parameter
png(file="out/fig_07.png", width=8, height=4, units="in", res=600)
par(mfcol=c(5,2), mar=c(3, 3, 1.5, 0.1), mgp=c(1.5, 0.5, 0), cex=0.5)
hist(post.params.ar3$u, 40, col = "grey", xlab = "u", main = "AR(3)")
hist(post.params.ar3$q, 40, col = "grey", xlab = "q", main = "AR(3)")
hist(post.params.ar3$b1, 40, col = "grey", xlab = "b1", main = "AR(3)")
hist(post.params.ar3$b2, 40, col = "grey", xlab = "b2", main = "AR(3)")
hist(post.params.ar3$b3, 40, col = "grey", xlab = "b3", main = "AR(3)")
hist(post.params.ar1$u, 40, col = "grey", xlab = "u", main = "AR(1)")
hist(post.params.ar1$q, 40, col = "grey", xlab = "q", main = "AR(1)")
hist(post.params.ar1$b, 40, col = "grey", xlab = "b", main = "AR(1)")
dev.off()

## Model Diagnostics

# Running another JAGS model for diagnostics of AR(1)
jags.params.ar1 <- c("q", "u", "b")
mod_diag_ar1 <- R2jags::jags(jags.data, parameters.to.save = jags.params.ar1, 
                             model.file = model.loc.ar1, n.chains = 3, n.burnin = 5000, n.thin = 5, 
                             n.iter = 10000, DIC = TRUE)

# Running another JAGS model for diagnostics of AR(3)
jags.params.ar3 <- c("q", "u", "b1", "b2", "b3")
mod_diag_ar3 <- R2jags::jags(jags.data, parameters.to.save = jags.params.ar3, 
                             model.file = model.loc.ar3, n.chains = 3, n.burnin = 5000, n.thin = 5, 
                             n.iter = 10000, DIC = TRUE)

# Function to create an MCMC list from JAGS model
createMcmcList <- function(jagsmodel) {
  McmcArray <- as.array(jagsmodel$BUGSoutput$sims.array)
  McmcList <- vector("list", length = dim(McmcArray)[2])
  for (i in 1:length(McmcList)) McmcList[[i]] <- as.mcmc(McmcArray[,i, ])
  McmcList <- mcmc.list(McmcList)
  return(McmcList)
}

# Creating an MCMC list from the diagnostics AR(1) model
McmcList_ar1 <- createMcmcList(mod_diag_ar1)

# Summarizing the MCMC list AR(1)
summary(McmcList_ar1)

# Creating an MCMC list from the diagnostics AR(3) model
McmcList_ar3 <- createMcmcList(mod_diag_ar3)

# Summarizing the MCMC list of AR(3)
summary(McmcList_ar3)

# Setting up the plot area for diagnostic plots
png(file="out/fig_08.png", width=8, height=4, units="in", res=600)
par(mar=c(3, 3, 2, 0.1), mgp=c(1.5, 0.5, 0), cex=0.7)
# Plotting diagnostic plots of AR(1) model
plot(McmcList_ar1)
dev.off()

# Setting up the plot area for diagnostic plots
png(file="out/fig_09.png", width=8, height=4, units="in", res=600)
par(mar=c(3, 3, 2, 0.1), mgp=c(1.5, 0.5, 0), cex=0.7)
# Plotting diagnostic plots of AR(3) model
plot(McmcList_ar3[,1:3])
dev.off()

# Setting up the plot area for diagnostic plots
png(file="out/fig_10.png", width=8, height=4, units="in", res=600)
par(mar=c(3, 3, 2, 0.1), mgp=c(1.5, 0.5, 0), cex=0.7)
# Plotting diagnostic plots of AR(3) model
plot(McmcList_ar3[,4:6])
dev.off()


## Forecast

# Updating JAGS data for two-steps ahead forecasts
jags.data <- list(X = c(gnp, NA, NA), N = N+2)

# Running the JAGS AR(1) model for two-steps ahead forecasts
jags.params.ar1 <- c("q", "b", "u", "EX", "X")
model.loc.ar3 <- ("AR1_gnp.txt")
mod_ss_forecast_ar1 <- jags(jags.data, parameters.to.save = jags.params.ar1, 
                            model.file = model.loc.ar1, n.chains = 3, n.burnin = 5000, n.thin = 5, 
                            n.iter = 10000, DIC = TRUE)

# Running the JAGS AR(3) model for two-steps ahead forecasts
jags.params.ar3 <- c("q", "u", "b1", "b2", "b3", "EX", "X")
model.loc.ar3 <- ("AR3_gnp.txt")
mod_ss_forecast_ar3 <- jags(jags.data, parameters.to.save = jags.params.ar3, 
                            model.file = model.loc.ar3, n.chains = 3, n.burnin = 5000, n.thin = 5, 
                            n.iter = 10000, DIC = TRUE)

# Getting quantiles and means of the predictions
EX_ar1 <- mod_ss_forecast_ar1$BUGSoutput$sims.list$EX 
X_ar1 <- mod_ss_forecast_ar1$BUGSoutput$sims.list$X

plot_df_ar1 <- as.data.frame(cbind(apply(EX_ar1, 2, quantile, 0.025), 
                                   apply(EX_ar1, 2, mean), 
                                   apply(EX_ar1, 2, quantile, 0.975)))
plot_df_ar1 <- rbind(matrix(rep(NA_real_, 3), ncol=3), plot_df_ar1)
colnames(plot_df_ar1) <- c("lower", "mean_pred", "upper")
plot_df_ar1$fit <- apply(X_ar1, 2 , mean)
plot_df_ar1$Date <- c(q.gnp.df$Date, q.gnp.df$Date[dim(q.gnp.df)[1]] %m+% months(3), 
                      q.gnp.df$Date[dim(q.gnp.df)[1]] %m+% months(6))
plot_df_ar1$obs <- c(q.gnp.df$GNP_US, NA, NA)

plot_forecast_ar1 <- plot_df_ar1 %>% 
  ggplot() +
  geom_point(aes(x = Date, y = obs), alpha = 0.8, size = 0.5, col = "black") +
  geom_line(aes(x = Date, y = fit), alpha = 1, linewidth = 0.3, col = "black") +
  geom_line(aes(x = Date, y = mean_pred), alpha = 0.8, linewidth = 0.3,  col = "red") +
  ylab("Growth Rate") +
  scale_y_continuous(n.breaks = 12) +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) +
  ggtitle("Forecast view of AR(1) model")

# Getting quantiles and means of the predictions
EX_ar3 <- mod_ss_forecast_ar3$BUGSoutput$sims.list$EX 
X_ar3 <- mod_ss_forecast_ar3$BUGSoutput$sims.list$X

plot_df_ar3 <- as.data.frame(cbind(apply(EX_ar3, 2, quantile, 0.025), 
                                   apply(EX_ar3, 2, mean), 
                                   apply(EX_ar3, 2, quantile, 0.975)))
plot_df_ar3 <- rbind(matrix(rep(NA_real_, 9), ncol=3), plot_df_ar3)
colnames(plot_df_ar3) <- c("lower", "mean_pred", "upper")
plot_df_ar3$fit <- apply(X_ar3, 2 , mean)
plot_df_ar3$Date <- c(q.gnp.df$Date, q.gnp.df$Date[dim(q.gnp.df)[1]] %m+% months(3), 
                      q.gnp.df$Date[dim(q.gnp.df)[1]] %m+% months(6))
plot_df_ar3$obs <- c(q.gnp.df$GNP_US, NA, NA)

plot_forecast_ar3 <- plot_df_ar3 %>% 
  ggplot() +
  geom_point(aes(x = Date, y = obs), alpha = 0.8, size = 0.5, col = "black") +
  geom_line(aes(x = Date, y = fit), alpha = 1, linewidth = 0.3, col = "black") +
  geom_line(aes(x = Date, y = mean_pred), alpha = 0.8, linewidth = 0.3,  col = "red") +
  ylab("Growth Rate") +
  scale_y_continuous(n.breaks = 12) +
  theme_bw() +
  theme(
    text = element_text(size = 8),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  ) +
  ggtitle("Forecast view of AR(3) model")

png(file="out/fig_11.png", width=8, height=4, units="in", res=600)
grid.arrange(plot_forecast_ar1, plot_forecast_ar3, ncol=2, widths=c(1,1))
dev.off()

