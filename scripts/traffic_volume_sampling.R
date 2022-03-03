library(HMCBRR)

seed <- 2022
nsamples <- 30000
warmup <- 20000
nu <- 6
sigma <- 1
prior_means <- rep(0, 8)
prior_variances <- rep(10, 8)
n_cores <- parallel::detectCores()

##### Loading in Data #####

original_data <- read.csv('scripts/traffic_volume.csv')
# holiday: 1 for yes (i.e. not "None") or 0 for no
traffic_volume <- data.frame(holiday = as.numeric(original_data$holiday!="None"))
# temp: continuous variable
temp_mean <- mean(original_data$temp)
temp_sd <- sd(original_data$temp)
# traffic_volume <- data.frame('temp' = (original_data$temp-temp_mean)/temp_sd)
traffic_volume$temp <- (original_data$temp-temp_mean)/temp_sd
# rain: continuous variable
rain_mean <- mean(original_data$rain_1h)
rain_sd <- sd(original_data$rain_1h)
traffic_volume$rain <- (original_data$rain_1h-rain_mean)/rain_sd
# snow: continuous variable
snow_mean <- mean(original_data$snow_1h)
snow_sd <- sd(original_data$snow_1h)
traffic_volume$snow <- (original_data$snow_1h-snow_mean)/snow_sd
# clouds: continuous variable
clouds_mean <- mean(original_data$clouds_all)
clouds_sd <- sd(original_data$clouds_all)
traffic_volume$clouds <- (original_data$clouds_all-clouds_mean)/clouds_sd
# weather_main: 1 for good (i.e. "Clear" or "Clouds") or 0 for bad
traffic_volume$weather <- as.numeric(original_data$weather_main %in% c("Clear", "Clouds"))
# time: 1 for rush hour (i.e. between 07:00-09:00 and 16:00-19:00) or 0 for no
original_data$hour <- format(strptime(original_data$date_time, format = "%Y-%m-%d %H:%M:%S"), format = "%H")
traffic_volume$rush_hour <- as.numeric(original_data$hour %in% c("07", "08", "09", "16", "17", "18", "19"))
# y (response variable): traffic volume
traffic_volume$y <- original_data$traffic_volume
rm(original_data)

##### Sampling from full posterior #####

traffic_volume_data <- list()
traffic_volume_data$y <- traffic_volume$y
traffic_volume_data$X <- as.matrix(cbind(rep(1, nrow(traffic_volume[,1:7])), traffic_volume[,1:7]))
colnames(traffic_volume_data$X)[1] <- 'intercept'

full_posterior <- hmc_sample_BRR(y = traffic_volume_data$y,
                                 X = traffic_volume_data$X,
                                 C = 1,
                                 nu = nu,
                                 sigma = sigma,
                                 prior_means = prior_means,
                                 prior_variances = prior_variances,
                                 iterations = nsamples + warmup,
                                 warmup = warmup,
                                 chains = 1,
                                 seed = seed,
                                 output = T)

##### Sampling from sub-posterior C=2 #####

data_split_2 <- split_data(traffic_volume,
                           y_col_index = 8,
                           X_col_index = 1:7,
                           C = 2,
                           as_dataframe = F)
sub_posteriors_2 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                         data_split = data_split_2,
                                         C = 2,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Sampling from sub-posterior C=4 #####

data_split_4 <- split_data(traffic_volume,
                           y_col_index = 8,
                           X_col_index = 1:7,
                           C = 4,
                           as_dataframe = F)
sub_posteriors_4 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                         data_split = data_split_4,
                                         C = 4,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Sampling from sub-posterior C=8 #####

data_split_8 <- split_data(traffic_volume,
                           y_col_index = 8,
                           X_col_index = 1:7,
                           C = 8,
                           as_dataframe = F)
sub_posteriors_8 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                         data_split = data_split_8,
                                         C = 8,
                                         nu = nu,
                                         sigma = sigma,
                                         prior_means = prior_means,
                                         prior_variances = prior_variances,
                                         warmup = warmup,
                                         seed = seed,
                                         output = T)

##### Sampling from sub-posterior C=16 #####

data_split_16 <- split_data(traffic_volume,
                            y_col_index = 8,
                            X_col_index = 1:7,
                            C = 16,
                            as_dataframe = F)
sub_posteriors_16 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                          data_split = data_split_16,
                                          C = 16,
                                          nu = nu,
                                          sigma = sigma,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = warmup,
                                          seed = seed,
                                          output = T)

##### Sampling from sub-posterior C=32 #####

data_split_32 <- split_data(traffic_volume,
                            y_col_index = 8,
                            X_col_index = 1:7,
                            C = 32,
                            as_dataframe = F)
sub_posteriors_32 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                          data_split = data_split_32,
                                          C = 32,
                                          nu = nu,
                                          sigma = sigma,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = warmup,
                                          seed = seed,
                                          output = T)

##### Sampling from sub-posterior C=64 #####

data_split_64 <- split_data(traffic_volume,
                            y_col_index = 8,
                            X_col_index = 1:7,
                            C = 64,
                            as_dataframe = F)
sub_posteriors_64 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                          data_split = data_split_64,
                                          C = 64,
                                          nu = nu,
                                          sigma = sigma,
                                          prior_means = prior_means,
                                          prior_variances = prior_variances,
                                          warmup = warmup,
                                          seed = seed,
                                          output = T)

##### Sampling from sub-posterior C=128 #####

data_split_128 <- split_data(traffic_volume,
                             y_col_index = 8,
                             X_col_index = 1:7,
                             C = 128,
                             as_dataframe = F)
sub_posteriors_128 <- hmc_base_sampler_BRR(nsamples = nsamples,
                                           data_split = data_split_128,
                                           C = 128,
                                           nu = nu,
                                           sigma = sigma,
                                           prior_means = prior_means,
                                           prior_variances = prior_variances,
                                           warmup = warmup,
                                           seed = seed,
                                           output = T)
