#' HMC sampler for Robust Regression model
#'
#' Sample from (sub-)posterior using Stan
#'
#' @param full_data_count a matrix or dataframe of the unique data with their
#'                        corresponding counts
#' @param C number of sub-posterior (default to 1)
#' @param nu degrees of freedom in t-distribution (must be greater than 2) to
#'           ensure existence of variance
#' @param sigma scale parameter in t-distribution
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param iterations number of iterations per chain
#' @param warmup number of burn in iterations
#' @param chains number of chains
#' @param seed seed number for random number generation
#' @param output boolean value: defaults to T, determines whether or not to
#'               print output to console
#'
#' @return samples from the (sub-)posterior target for the robust regression model
#'
#' @export
hmc_sample_BRR <- function(full_data_count,
                           C,
                           nu,
                           sigma,
                           prior_means,
                           prior_variances,
                           iterations,
                           warmup,
                           chains,
                           seed = sample.int(.Machine$integer.max, 1),
                           output = F) {
  if (!is.matrix(full_data_count) & !is.data.frame(full_data_count)) {
    stop("hmc_sample_BRR: full_data_count must be a matrix or data frame")
  } else if (!is.vector(prior_means)) {
    stop("hmc_sample_BRR: prior_means must be a vector")
  } else if (!is.vector(prior_variances)) {
    stop("hmc_sample_BRR: prior_variances must be a vector")
  } else if (nu <= 2) {
    stop("hmc_sample_BRR: nu must be greater than 2 to ensure existence of variance")
  } else if (sigma <= 0) {
    stop("hmc_sample_BRR: sigma must be greater than 0")
  }
  y <- full_data_count$y
  X <- as.matrix(subset(full_data_count, select = -c(y, count)))
  count <- full_data_count$count
  if (length(y) != nrow(X)) {
    stop("hmc_sample_BRR: y and X do not have the same number of samples")
  }
  dim <- ncol(X)
  if (length(prior_means)!=dim) {
    stop("hmc_sample_BRR: prior_means must be a vector of length ncol(X)")
  } else if (length(prior_variances)!=dim) {
    stop("hmc_sample_BRR: prior_variances must be a vector of length ncol(X)")
  }
  # check that the design matrix does have the intercept (i.e. first column not identical to a vector of 1s)
  # reset rownames so that we can compare the first column to a vector of 1s
  rownames(X) <- c()
  if (!identical(X[,1], rep(1, nrow(X)))) {
    X <- cbind(rep(1, nrow(X)), X)
    colnames(X)[1] <- 'intercept'
    warning("hmc_sample_BRR: the first column of the design matrix was not a column of 1s -
            the design matrix has been changed to include the intercept")
  }
  print("Sampling from robust regression model")
  training_data <- list(nsamples = length(y),
                        p = (ncol(X)-1),
                        y = y,
                        X = X,
                        count = count,
                        prior_means = prior_means,
                        prior_variances = prior_variances,
                        C = C,
                        nu = nu,
                        sigma = sigma)
  if (output) {
    model <- rstan::sampling(object = stanmodels$bayes_robust_reg,
                             data = training_data,
                             iter = iterations,
                             warmup = warmup,
                             chains = chains,
                             seed = seed,
                             control = list(adapt_delta = 0.99,
                                            max_treedepth = 20))
  } else {
    model <- rstan::sampling(object = stanmodels$bayes_robust_reg,
                             data = training_data,
                             iter = iterations,
                             warmup = warmup,
                             chains = chains,
                             verbose = FALSE,
                             refresh = 0,
                             seed = seed,
                             control = list(adapt_delta = 0.99,
                                            max_treedepth = 20))
  }
  print('Finished sampling from robust regression model')
  return(rstan::extract(model)$beta)
}

#' HMC sampler for base level for robust regression model
#'
#' Sample for base level for robust regression
#'
#' @param nsamples number of samples per node
#' @param warmup number of burn in iterations
#' @param data list of length C where each item is a list where for c=1,...,C,
#'             data_split[[c]]$full_data_count is a matrix or data frame of
#'             the unique data with their corresponding counts
#' @param C number of sub-posteriors (default to 1)]
#' @param nu degrees of freedom in t-distribution (must be greater than 2) to
#'           ensure existence of variance
#' @param sigma scale parameter in t-distribution
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param seed seed number for random number generation
#' @param output boolean value: defaults to T, determines whether or not
#'               to print output to console
#'
#' @return samples from the sub-posterior targets for the split data sets
#'         for the robust regression model
#'
#' @export
hmc_base_sampler_BRR <- function(nsamples,
                                 warmup,
                                 data_split,
                                 C,
                                 nu,
                                 sigma,
                                 prior_means,
                                 prior_variances,
                                 seed = sample.int(.Machine$integer.max, 1),
                                 output = F) {
  if (!is.list(data_split)) {
    stop("hmc_base_sampler_BRR: data_split must be a list")
  } else if (!all(sapply(1:length(data_split), function(i) is.data.frame(data_split[[i]]$full_data_count)))) {
    stop("hmc_base_sampler_BRR: for each i in 1:length(data_split), data_split[[i]]$full_data_count must be a data frame")
  } else if (!is.vector(prior_means)) {
    stop("hmc_base_sampler_BRR: prior_means must be a vector")
  } else if (!is.vector(prior_variances)) {
    stop("hmc_base_sampler_BRR: prior_variances must be a vector")
  }  else if (nu <= 2) {
    stop("hmc_base_sampler_BRR: nu must be greater than 2 to ensure existence of variance")
  } else if (sigma <= 0) {
    stop("hmc_base_sampler_BRR: sigma must be greater than 0")
  }
  cl <- parallel::makeCluster(length(data_split),
                              setup_strategy = "sequential",
                              outfile = 'output_hmc_sample_BRR.txt')
  parallel::clusterExport(cl, envir = environment(), varlist = c(ls(), "hmc_sample_BRR", "seed"))
  base_samples <- parallel::parLapply(cl, X = 1:length(data_split), fun = function(c) {
    hmc_sample_BRR(full_data_count = data_split[[c]]$full_data_count,
                   C = C,
                   nu = nu,
                   sigma = sigma,
                   prior_means = prior_means,
                   prior_variances = prior_variances,
                   iterations = nsamples+warmup,
                   warmup = warmup,
                   chains = 1,
                   seed = seed,
                   output = output)})
  parallel::stopCluster(cl)
  return(base_samples)
}
