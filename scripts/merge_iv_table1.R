# Run here: Rscript scripts/merge_iv_table1.R results/results_iv_table1_07_n1000 1000
# Arguments: output directory, number of simulations.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1) args[1] else "results/results_iv_table1_07_n1000"
M <- if (length(args) >= 2) as.numeric(args[2]) else 1000
stopifnot(is.finite(M), M == floor(M), M %in% 1:1000)
M <- as.integer(M)
targets <- c("B(0.5)", "B(1)", "B(1.5)", "beta (pooled)")
runs <- vector("list", M)
for (id in seq_len(M)) {
  file <- file.path(output, sprintf("simulation_%04d.rds", id))
  if (!file.exists(file)) stop("Missing simulation: ", file)
  run <- readRDS(file)
  if (id == 1) settings <- run$settings
  stopifnot(
    identical(settings$variance_method, "sum(IF^2)/n^2"),
    identical(run$settings, settings), settings$generator == "IV.R",
    settings$estimator == "estimation_iv", run$seed == id,
    nrow(run$results) == 8, all(run$results$simulation == id),
    identical(run$results$cause, rep(1:2, each = 4)),
    identical(run$results$parameter, rep(targets, 2)),
    isTRUE(all.equal(run$parameters$beta, c(-.25, -.125) * settings$event_rate)),
    isTRUE(all.equal(run$results$truth,
                    rep(run$parameters$beta, each = 4) * rep(c(.5, 1, 1.5, 1), 2)))
  )
  runs[[id]] <- run
}
names(runs) <- sprintf("simulation_%04d", seq_len(M))
results <- do.call(rbind, lapply(runs, function(x) x$results))
finite_stat <- function(x, fun) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  fun(x)
}
rows <- list()
for (cause in 1:2) for (target in targets) {
  x <- results[results$cause == cause & results$parameter == target, ]
  valid <- is.finite(x$lower) & is.finite(x$upper)
  rows[[length(rows) + 1]] <- data.frame(
    event_rate = settings$event_rate, n = settings$n, cause = cause,
    parameter = target, truth = x$truth[1], simulations = M,
    successful_estimates = sum(is.finite(x$estimate)),
    successful_intervals = sum(valid),
    mean_estimate = finite_stat(x$estimate, mean),
    Bias = finite_stat(x$estimate - x$truth, mean),
    SEE = finite_stat(x$estimate, sd), SD = finite_stat(x$std, mean),
    CP = if (any(valid)) 100 * mean(x$lower[valid] <= x$truth[valid] &
                                    x$truth[valid] <= x$upper[valid]) else NA_real_
  )
}
summary <- do.call(rbind, rows)
metrics <- c("Bias", "SEE", "SD", "CP")
table1 <- do.call(rbind, lapply(1:2, function(cause) {
  x <- summary[summary$cause == cause, ]
  values <- as.data.frame(t(as.matrix(x[metrics])))
  names(values) <- targets
  rownames(values) <- NULL
  data.frame(event_rate = settings$event_rate, n = settings$n, cause = cause,
             metric = c("Bias", "SEE", "SD", "95% CP"), values, check.names = FALSE)
}))
file <- file.path(output, sprintf("iv_table1_%02d_n%d.rds",
                                round(100 * settings$event_rate), settings$n))
temporary <- tempfile(".merged_", tmpdir = output)
saveRDS(list(settings = settings, runs = runs, summary = summary, table1 = table1), temporary)
if (!file.rename(temporary, file)) stop("Could not save ", file)
print(table1)
cat("Saved:", file, "\n")
