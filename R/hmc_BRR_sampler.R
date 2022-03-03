#' HMC sampler for Robust Regression model
#'
#' Sample from (sub-)posterior using Stan
#'
#' @param noise_error the distribution for the residuals of the regression model.
#'                    Must be either "student_t" or "laplace" (student_t by default)
#' @param y vector of responses
#' @param X design matrix
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
hmc_sample_BRR <- function(noise_error = 'student_t',
                           y,
                           X,
                           C,
                           nu = NULL,
                           sigma,
                           prior_means,
                           prior_variances,
                           iterations,
                           warmup,
                           chains,
                           seed = sample.int(.Machine$integer.max, 1),
                           output = F) {
  if (!(noise_error %in% c('student_t', 'laplace'))) {
    stop("hmc_sample_BRR: noise_error must be \"student_t\" or \"laplace\"")
  } else if (!is.vector(prior_means)) {
    stop("hmc_sample_BRR: prior_means must be a vector")
  } else if (!is.vector(prior_variances)) {
    stop("hmc_sample_BRR: prior_variances must be a vector")
  } 
  if (noise_error == "student_t") {
    if (!is.numeric(nu)) {
      stop("hmc_sample_BRR: if noise_error is \"student_t\", nu must be a numeric value greater than 2")
    } else if (nu <= 2) {
      stop("hmc_sample_BRR: nu must be greater than 2 to ensure existence of variance")
    }
  }
  if (sigma <= 0) {
    stop("hmc_sample_BRR: sigma must be greater than 0")
  }
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
  if (noise_error == "student_t") {
    training_data <- list(nsamples = length(y),
                          p = (ncol(X)-1),
                          y = y,
                          X = X,
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
  } else if (noise_error == "laplace") {
    training_data <- list(nsamples = length(y),
                          p = (ncol(X)-1),
                          y = y,
                          X = X,
                          prior_means = prior_means,
                          prior_variances = prior_variances,
                          C = C,
                          sigma = sigma)
    if (output) {
      model <- rstan::sampling(object = stanmodels$bayes_quantile_reg,
                               data = training_data,
                               iter = iterations,
                               warmup = warmup,
                               chains = chains,
                               seed = seed,
                               control = list(adapt_delta = 0.99,
                                              max_treedepth = 20))
    } else {
      model <- rstan::sampling(object = stanmodels$bayes_quantile_reg,
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
  }
  print('Finished sampling from robust regression model')
  return(rstan::extract(model)$beta)
}

#' HMC sampler for base level for robust regression model
#'
#' Sample for base level for robust regression
#'
#' @param noise_error the distribution for the residuals of the regression model.
#'                    Must be either "student_t" or "laplace" (student_t by default)
#' @param nsamples number of samples per node
#' @param warmup number of burn in iterations
#' @param data_split list of length C where each item is a list where for c=1,...,C,
#'                   data_split[[c]]$y is a vector of responses and data_split[[c]]$X
#'                   is the design matrix
#' @param C number of sub-posteriors (default to 1)]
#' @param nu degrees of freedom in t-distribution (must be greater than 2) to
#'           ensure existence of variance
#' @param sigma scale parameter in t-distribution
#' @param prior_means prior for means of predictors
#' @param prior_variances prior for variances of predictors
#' @param seed seed number for random number generation
#' @param n_cores number of cores to use
#' @param output boolean value: defaults to T, determines whether or not
#'               to print output to console
#'
#' @return samples from the sub-posterior targets for the split data sets
#'         for the robust regression model
#'
#' @export
hmc_base_sampler_BRR <- function(noise_error = 'student_t',
                                 nsamples,
                                 warmup,
                                 data_split,
                                 C,
                                 nu = NULL,
                                 sigma,
                                 prior_means,
                                 prior_variances,
                                 seed = sample.int(.Machine$integer.max, 1),
                                 n_cores = parallel::detectCores(),
                                 output = F) {
  if (!(noise_error %in% c('student_t', 'laplace'))) {
    stop("hmc_base_sampler_BRR: noise_error must be \"student_t\" or \"laplace\"")
  } else if (!is.list(data_split)) {
    stop("hmc_base_sampler_BRR: data_split must be a list")
  } else if (!is.vector(prior_means)) {
    stop("hmc_base_sampler_BRR: prior_means must be a vector")
  } else if (!is.vector(prior_variances)) {
    stop("hmc_base_sampler_BRR: prior_variances must be a vector")
  } 
  if (noise_error == "student_t") {
    if (!is.numeric(nu)) {
      stop("hmc_base_sampler_BRR: if noise_error is \"student_t\", nu must be a numeric value greater than 2")
    } else if (nu <= 2) {
      stop("hmc_base_sampler_BRR: nu must be greater than 2 to ensure existence of variance")
    }
  }
  if (sigma <= 0) {
    stop("hmc_base_sampler_BRR: sigma must be greater than 0")
  }
  cl <- parallel::makeCluster(n_cores,
                              setup_strategy = "sequential",
                              outfile = 'output_hmc_sample_BRR.txt')
  parallel::clusterExport(cl, envir = environment(), varlist = c(ls(), "hmc_sample_BRR", "seed"))
  base_samples <- parallel::parLapply(cl, X = 1:length(data_split), fun = function(c) {
    hmc_sample_BRR(y = data_split[[c]]$y,
                   X = data_split[[c]]$X,
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
