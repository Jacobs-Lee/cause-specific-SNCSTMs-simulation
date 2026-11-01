# Rscript scripts/run_si_cumulative.R 1 1000 results/results_si_cumulative_07_n1000 0.07
# Arguments: simulation ID, sample size, output directory, event rate.
args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])
n <- if (length(args) >= 2) as.numeric(args[2]) else 1000
event_rate <- if (length(args) >= 4) as.numeric(args[4]) else .07
stopifnot(length(args) %in% 1:4, is.finite(id), id == floor(id), id %in% 1:1000,
          n %in% c(1000, 5000, 10000), event_rate %in% c(.07, .15, .30, .60))
id <- as.integer(id)
n <- as.integer(n)
output <- if (length(args) >= 3) args[3] else
  sprintf("results/results_si_cumulative_%02d_n%d", round(100 * event_rate), n)
stopifnot(nzchar(output))

source("R/SI.R")
sources <- c("R/SI.R", "scripts/run_si_cumulative.R")
settings <- list(n = n, event_rate = event_rate, generator = "SI.R",
                 estimator = "estimation_si", n_grid = 10L,
                 B_se = "sqrt(sum(psi^2)/n^2)", pooled_se = "sqrt(sum(psi^2)/n^2)",
                 variance_method = "sum(IF^2)/n^2",
                 source_md5 = tools::md5sum(sources))
file <- file.path(output, sprintf("simulation_%04d.rds", id))
if (file.exists(file) && !identical(readRDS(file)$settings, settings))
  stop("Existing result has different settings/code; use a new output directory.")

started <- proc.time()["elapsed"]
simulation <- simulation_si(n = n, seed = id, event_rate = event_rate)
data <- simulation$observed
parameters <- simulation$parameters
rm(simulation)
targets <- c("B_D(0.5)", "B_D(1)", "B_D(1.5)", "beta (pooled)")
results <- data.frame(
  simulation = id, cause = rep(1:2, each = 4), parameter = rep(targets, 2),
  truth = rep(parameters$beta, each = 4) * rep(c(.5, 1, 1.5, 1), 2),
  estimate = NA_real_, std = NA_real_, variance = NA_real_,
  lower = NA_real_, upper = NA_real_, error = ""
)

# Shared nuisance fit; each cause uses the same data and original row order.
nuisance <- tryCatch(fit_si_nuisance(data), error = function(e) e)
columns <- c("estimate", "std", "variance", "lower", "upper")
for (cause in 1:2) {
  rows <- which(results$cause == cause)
  fit <- if (inherits(nuisance, "error")) nuisance else tryCatch(
    estimation_si(data, cause = cause, n_grid = 10, nuisance = nuisance),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    results$error[rows] <- conditionMessage(fit)
    next
  }
  B <- fit$B_summary
  pooled <- fit$pooled
  stopifnot(identical(B$time, c(.5, 1, 1.5)), all(B$cause == cause),
            nrow(pooled) == 1, pooled$cause == cause)
  x <- data.frame(estimate = c(B$B, pooled$beta), std = c(B$std_B, pooled$std),
                  variance = c(B$var_B, pooled$var_beta))
  x$lower <- x$estimate - 1.96 * x$std
  x$upper <- x$estimate + 1.96 * x$std
  results[rows, columns] <- x
  invalid <- !apply(is.finite(as.matrix(x)), 1, all) | x$std < 0 | x$variance < 0
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
