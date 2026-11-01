# Rscript scripts/merge_si_cumulative.R results/results_si_cumulative_07_n1000 1000
# Arguments: output directory, expected simulations (2 for a pilot).
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args) >= 1) args[1] else "results/results_si_cumulative_07_n1000"
M <- if (length(args) >= 2) as.numeric(args[2]) else 1000
stopifnot(length(args) <= 2, is.finite(M), M == floor(M), M %in% 1:1000)
M <- as.integer(M)
targets <- c("B_D(0.5)", "B_D(1)", "B_D(1.5)", "beta (pooled)")
runs <- vector("list", M)
for (id in seq_len(M)) {
  file <- file.path(output, sprintf("simulation_%04d.rds", id))
  if (!file.exists(file)) stop("Missing simulation: ", file)
  run <- readRDS(file)
  if (id == 1) settings <- run$settings
  stopifnot(
    identical(settings$variance_method, "sum(IF^2)/n^2"),
    identical(run$settings, settings), settings$generator == "SI.R",
    settings$estimator == "estimation_si", settings$n_grid == 10L,
    settings$B_se == "sqrt(sum(psi^2)/n^2)", settings$pooled_se == "sqrt(sum(psi^2)/n^2)",
    settings$n %in% c(1000, 5000, 10000), settings$event_rate %in% c(.07, .15, .30, .60),
    length(settings$source_md5) == 2, all(!is.na(settings$source_md5)), run$seed == id,
    nrow(run$results) == 8, all(run$results$simulation == id),
    identical(run$results$cause, rep(1:2, each = 4)),
    identical(run$results$parameter, rep(targets, 2)),
    isTRUE(all.equal(run$parameters$beta, c(-.25, -.125) * settings$event_rate)),
    isTRUE(all.equal(run$results$truth,
                    rep(run$parameters$beta, each = 4) * rep(c(.5, 1, 1.5, 1), 2)))
  )
  x <- run$results
  values <- as.matrix(x[c("estimate", "std", "variance", "lower", "upper")])
  success <- x$error == ""
  stopifnot(!anyNA(x$error), all(is.finite(values[success, , drop = FALSE])),
            all(is.na(values[!success, , drop = FALSE])),
            all(x$std[success] >= 0), all(x$variance[success] >= 0),
            isTRUE(all.equal(x$variance[success], x$std[success]^2)),
            isTRUE(all.equal(x$lower[success], x$estimate[success] - 1.96 * x$std[success])),
            isTRUE(all.equal(x$upper[success], x$estimate[success] + 1.96 * x$std[success])))
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
    successful_estimates = sum(is.finite(x$estimate)), successful_intervals = sum(valid),
    failed_estimates = sum(x$error != ""),
    mean_estimate = finite_stat(x$estimate, mean), Bias = finite_stat(x$estimate - x$truth, mean),
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
file <- file.path(output, sprintf("si_cumulative_%02d_n%d.rds",
                                round(100 * settings$event_rate), settings$n))
temporary <- tempfile(".merged_", tmpdir = output)
saveRDS(list(settings = settings, runs = runs, summary = summary, table1 = table1), temporary)
if (!file.rename(temporary, file)) stop("Could not save ", file)
print(table1)
cat("Saved:", file, "\n")
