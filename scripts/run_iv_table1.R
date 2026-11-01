# Run here: Rscript scripts/run_iv_table1.R 1 1000 results/results_iv_table1_07_n1000 0.07
# Arguments: simulation ID, sample size, output directory, event rate.
args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])
n <- if (length(args) >= 2) as.numeric(args[2]) else 1000
event_rate <- if (length(args) >= 4) as.numeric(args[4]) else .07
stopifnot(length(id) == 1, is.finite(id), id == floor(id), id %in% 1:1000,
          n %in% c(1000, 5000, 10000), event_rate %in% c(.07, .15, .30, .60))
id <- as.integer(id)
n <- as.integer(n)
output <- if (length(args) >= 3) args[3] else
  sprintf("results/results_iv_table1_%02d_n%d", round(100 * event_rate), n)

source("R/IV.R")
sources <- c("R/IV.R", "scripts/run_iv_table1.R")
settings <- list(n = n, event_rate = event_rate, generator = "IV.R",
                 estimator = "estimation_iv", n_grid = 10L,
                 variance_method = "sum(IF^2)/n^2",
                 source_md5 = tools::md5sum(sources))
file <- file.path(output, sprintf("simulation_%04d.rds", id))
if (file.exists(file) && !identical(readRDS(file)$settings, settings))
  stop("Existing result has different settings/code; use a new output directory.")

started <- proc.time()["elapsed"]
simulation <- simulation_iv(n = n, seed = id, event_rate = event_rate)
data <- simulation$observed
parameters <- simulation$parameters
rm(simulation)
targets <- c("B(0.5)", "B(1)", "B(1.5)", "beta (pooled)")
results <- data.frame(
  simulation = id, cause = rep(1:2, each = 4), parameter = rep(targets, 2),
  truth = rep(parameters$beta, each = 4) * rep(c(.5, 1, 1.5, 1), 2),
  estimate = NA_real_, std = NA_real_, variance = NA_real_,
  lower = NA_real_, upper = NA_real_, error = ""
)

# Fit nuisance once per dataset; both causes use the same fitted models.
nuisance <- tryCatch(fit_iv_nuisance(data), error = function(e) e)
columns <- c("estimate", "std", "variance", "lower", "upper")
for (cause in 1:2) {
  rows <- which(results$cause == cause)
  fit <- if (inherits(nuisance, "error")) nuisance else tryCatch(
    estimation_iv(data, cause = cause, n_grid = 10, nuisance = nuisance),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    results$error[rows] <- conditionMessage(fit)
    next
  }
  x <- fit$summary
  stopifnot(identical(x$parameter, targets), all(x$cause == cause))
  results[rows, columns] <- x[columns]
  invalid <- !apply(is.finite(as.matrix(x[columns])), 1, all) | x$std < 0 | x$variance < 0
  results[rows[invalid], columns] <- NA_real_
  results$error[rows[invalid]] <- "Non-finite result or negative variance/standard error"
}

dir.create(output, recursive = TRUE, showWarnings = FALSE)
run <- list(settings = settings, seed = id, parameters = parameters, results = results,
            seconds = unname(proc.time()["elapsed"] - started))
temporary <- tempfile(".simulation_", tmpdir = output)
saveRDS(run, temporary)
if (!file.rename(temporary, file)) stop("Could not save ", file)
print(results)
