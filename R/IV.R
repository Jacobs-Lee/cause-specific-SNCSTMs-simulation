# IV simulation: data generation, nuisance models, and the estimator used in the study.

## IV simulation with stagewise event realization

simulation_iv <- function(n = 10000, seed = 1, event_rate = 0.60) {
  set.seed(seed)

  times <- c(0, 0.5, 1, 1.5)
  tau <- 1.5
  beta <- c(-0.25, -0.125) * event_rate
  b <- c(-2.3307073, -1.4634162, -0.5399216, 0.8624681)[
    match(event_rate, c(.07, .15, .30, .60))
  ]
  if (is.na(b)) stop("event_rate must be 0.07, 0.15, 0.30, or 0.60")

  # CIF under treatment history D1,...,Dk
  cif <- function(t, L1, L2, U, D, k = 3) {
    rate <- c(0.025, 0.035)
    weight <- c(0.6, 0.4)
    m <- plogis(b + 0.05 * L1 + 0.05 * (L2 - 0.5) + 0.15 * U)

    F <- weight * m * (1 - exp(-rate * t)) / (1 - exp(-rate * tau))

    if (k == 0) return(F)
    if (t < times[k]) return(cif(t, L1, L2, U, D, k - 1))

    prev <- cif(t, L1, L2, U, D, k - 1)
    bound <- cif(times[k], L1, L2, U, D, k - 1)
    ell <- min(t, times[k + 1]) - times[k]
    e <- exp(-beta * D[k] * ell)

    prev + (1 - rev(bound) - prev) * (1 - e)
  }

  # Baseline variables and instrument
  L1 <- qnorm(runif(n, pnorm(-3), pnorm(3)))
  L2 <- rbinom(n, 1, 0.5)
  U <- sample(c(-1L, 1L), n, replace = TRUE)
  Z <- sample(0:2, n, replace = TRUE)

  # Treatment
  D <- matrix(0L, n, 3)
  treatment_lp <- 0.3 * L1 + 0.3 * (L2 - 0.5) + 0.8 * Z + 0.3 * U

  D[, 1] <- rbinom(n, 1, plogis(-1 + treatment_lp))
  for (k in 2:3) {
    pD <- plogis(
      -1.7 + 2 * D[, k - 1] +
        0.6 * rowMeans(D[, 1:(k - 1), drop = FALSE]) +
        treatment_lp
    )
    D[, k] <- rbinom(n, 1, pD)
  }

  # Stagewise realization of (T_tilde, V)
  event_time <- rep(Inf, n)
  event_type <- integer(n)
  U_E <- matrix(runif(3 * n), n, 3)
  U_V <- matrix(runif(3 * n), n, 3)

  for (i in 1:n) {
    for (k in 1:3) {
      t0 <- times[k]
      t1 <- times[k + 1]

      F0 <- cif(t0, L1[i], L2[i], U[i], D[i, ], k)
      F1 <- cif(t1, L1[i], L2[i], U[i], D[i, ], k)

      S0 <- 1 - sum(F0)
      p <- (F1 - F0) / S0

      if (U_E[i, k] > sum(p)) next
      j <- if (U_E[i, k] <= p[1]) 1L else 2L

      target <- F0[j] + U_V[i, k] * (F1[j] - F0[j])

      event_time[i] <- uniroot(
        function(t) cif(t, L1[i], L2[i], U[i], D[i, ], k)[j] - target,
        c(t0, t1), tol = 1e-9
      )$root

      event_type[i] <- j
      break
    }
  }

  # Cox censoring
  lp <- 0.3 * L1 + 0.3 * (L2 - 0.5)
  lambda0 <- 0.044591529489914
  C_star <- sqrt(-log(runif(n)) / (lambda0 * exp(lp)))
  C <- pmin(C_star, tau)

  T <- pmin(event_time, C)
  Delta <- as.integer(event_time <= C)
  V_obs <- ifelse(Delta == 1, event_type, 0L)

  full <- data.frame(
    id = 1:n, L1, L2, U, Z,
    D1 = D[, 1], D2 = D[, 2], D3 = D[, 3],
    event_time, event_type, C_star, C, T, Delta, V_obs
  )

  observed <- full[c(
    "id", "L1", "L2", "Z", "D1", "D2", "D3", "T", "Delta", "V_obs"
  )]

  for (k in 1:3)
    observed[observed$T < times[k], paste0("D", k)] <- NA_integer_

  list(
    observed = observed,
    full_data = full,
    parameters = list(
      times = times, tau = tau, beta = beta, b = b,
      censor_lambda0 = lambda0, censor_shape = 2
    )
  )
}

## Full IV influence functions: estimated instrument mean and Cox censoring.
## Source this file alone. Reuse fit_iv_nuisance(data) across causes/stages.
## All IFs have individual-contribution scale; variances are sum(IF^2) / n^2.

.iv_full_rows <- function(data) {
  columns <- c("Z", "L1", "L2", "T", "Delta", "V_obs")
  if (!is.data.frame(data) || nrow(data) < 2 || !all(columns %in% names(data)))
    stop("data must have at least two rows and columns: ", paste(columns, collapse = ", "))
  if (!all(vapply(data[columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data[columns]))))
    stop("Nuisance and outcome columns must be numeric and finite.")
  if (any(data$T < 0) || any(!data$Delta %in% 0:1) ||
      any(data$V_obs < 0 | data$V_obs != floor(data$V_obs)) ||
      any(data$Delta == 1 & data$V_obs == 0))
    stop("Invalid event times, event indicators, or observed causes.")
  data[c(columns, intersect("id", names(data)))]
}

.iv_full_inverse <- function(information, label) {
  if (any(!is.finite(information))) stop(label, " is non-finite.")
  tryCatch(solve(information), error = function(e) {
    stop(label, " is singular or numerically singular.", call. = FALSE)
  })
}

fit_iv_nuisance <- function(data) {
  rows <- .iv_full_rows(data)
  n <- nrow(data)
  mean_fit <- stats::lm(Z ~ L1 + L2, data = data)
  X <- stats::model.matrix(mean_fit)
  if (mean_fit$rank != ncol(X)) stop("The instrument mean design is singular.")
  alpha <- stats::coef(mean_fit)
  EZ <- drop(X %*% alpha)
  if_alpha <- (X * (data$Z - EZ)) %*%
    .iv_full_inverse(crossprod(X) / n, "Instrument mean information")

  x <- transform(data, censor_event = as.integer(Delta == 0 & T < 1.5))
  if (!any(x$censor_event)) stop("No censoring events are available to fit the Cox model.")
  censor_fit <- tryCatch(
    survival::coxph(survival::Surv(T, censor_event) ~ L1 + L2,
                   data = x, ties = "breslow", x = TRUE, y = TRUE),
    warning = function(w) stop("Censoring Cox fit: ", conditionMessage(w), call. = FALSE)
  )
  gamma <- stats::coef(censor_fit)
  if (any(!is.finite(gamma))) stop("The censoring Cox fit has undefined coefficients.")
  L <- censor_fit$x
  # Use the fitted times, including Cox's near-tie handling, for its risk sets.
  Tc <- censor_fit$y[, 1]
  censored <- censor_fit$y[, 2] == 1
  rr <- exp(drop(L %*% gamma))
  if (any(!is.finite(rr))) stop("Non-finite Cox relative risks.")
  censor_times <- sort(unique(Tc[censored]))
  order_c <- order(Tc)
  starts <- match(censor_times, Tc[order_c])
  risk_sum <- function(v) rev(cumsum(rev(v[order_c])))[starts]
  S0 <- risk_sum(rr) / n
  if (any(!is.finite(S0) | S0 <= 0)) stop("Unsupported Cox censoring risk set.")
  deaths <- tabulate(match(Tc[censored], censor_times), length(censor_times))
  dLambda <- deaths / (n * S0)
  Lbar <- cbind(risk_sum(rr * L[, 1]), risk_sum(rr * L[, 2])) / (n * S0)
  information <- matrix(0, 2, 2)
  for (j in 1:2) for (k in 1:2) {
    covariance <- risk_sum(rr * L[, j] * L[, k]) / (n * S0) - Lbar[, j] * Lbar[, k]
    information[j, k] <- sum(covariance * deaths / n)
  }

  # Zero-reference Breslow baseline; F = integral dLambda/S0, A = integral Lbar dLambda.
  Lambda <- c(0, cumsum(dLambda))
  F <- c(0, cumsum(dLambda / S0))
  A <- rbind(c(0, 0), cbind(cumsum(Lbar[, 1] * dLambda),
                           cumsum(Lbar[, 2] * dLambda)))
  index <- function(t) findInterval(t, censor_times) + 1L
  ic <- index(Tc)
  score <- -rr * (L * Lambda[ic] - A[ic, , drop = FALSE])
  score[censored, ] <- score[censored, , drop = FALSE] +
    L[censored, , drop = FALSE] - Lbar[match(Tc[censored], censor_times), , drop = FALSE]
  if_gamma <- score %*% .iv_full_inverse(information, "Censoring Cox information")
  if (any(!is.finite(if_alpha)) || any(!is.finite(if_gamma)))
    stop("Non-finite nuisance influence functions.")

  event_coefficient <- numeric(n)
  event_coefficient[censored] <- 1 / S0[match(Tc[censored], censor_times)]
  # Scalar t -> all subjects' baseline IFs, including the gamma correction.
  if_Lambda <- function(t) {
    event_coefficient * (Tc <= t) - rr * F[index(pmin(Tc, t))] -
      drop(if_gamma %*% A[index(t), ])
  }
  logG <- function(t, data) {
    covariates <- as.matrix(data[c("L1", "L2")])
    -Lambda[index(t)] * exp(drop(covariates %*% gamma))
  }

  # Cache original observation-time evaluations for weight derivatives.
  T <- data$T
  order_T <- order(T)
  lower <- findInterval(Tc, T[order_T], left.open = TRUE)  # Number of T_j < Tc_i.
  upper <- findInterval(Tc, T[order_T])                    # Number of T_j <= Tc_i.
  iT <- index(T)
  Lambda_T <- Lambda[iT]
  F_T <- F[iT]
  F_c <- F[ic]
  A_T <- A[iT, , drop = FALSE]

  # Return delta_i mean(b_j * a_j(t)) for every nuisance direction i.
  # b is held fixed; a already includes eligibility and competing-event weights.
  risk_if <- function(t, a, b = rep(1, n)) {
    v <- b * a * rr
    v[T >= t] <- 0  # Only prior competing events have nuisance-dependent weights.
    if (!any(v != 0)) return(numeric(n))
    total <- sum(v)
    prefix <- c(0, cumsum(v[order_T]))
    prefix_F <- c(0, cumsum((v * F_T)[order_T]))
    # sum_j v_j IF_Lambda_i(T_j), evaluated with prefix/suffix sums, not n x n.
    weighted_IF <- event_coefficient * (total - prefix[lower + 1L]) -
      rr * (prefix_F[upper + 1L] + F_c * (total - prefix[upper + 1L])) -
      drop(if_gamma %*% colSums(A_T * v))
    gamma_gradient <- colSums(L * (v * (Lambda[index(t)] - Lambda_T)))
    -(total * if_Lambda(t) - weighted_IF + drop(if_gamma %*% gamma_gradient)) / n
  }

  # Treat this object as immutable: the closures share the fitted IF information.
  structure(list(rows = rows, alpha = alpha, gamma = gamma, EZ = EZ, logG = logG,
                 if_alpha = if_alpha, if_gamma = if_gamma,
                 baseline = data.frame(time = censor_times, dLambda = dLambda,
                                       Lambda = Lambda[-1], S0 = S0),
                 if_Lambda = if_Lambda, risk_if = risk_if), class = "iv_nuisance")
}


## Projection weights are observed eligible * Y1, with no censoring adjustment.
## Treat these weights and the projection denominator as fixed when computing IFs.
## B_summary uses B(time-); a single stage starts from zero. Variance = sum(IF^2) / n^2.
.estimation_iv_core <- function(
    data, cause = 1, stage = NULL, n_grid = 10,
    nuisance = NULL, return_if_path = FALSE) {
  # Initialization: nuisance is still needed for B and its full IF.
  rows <- .iv_full_rows(data)
  n <- nrow(data)
  integer_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x)
  }
  if (!integer_scalar(cause) || cause < 1) stop("cause must be a positive integer.")
  if (!integer_scalar(n_grid) || n_grid < 1) stop("n_grid must be a positive integer.")
  if (!is.null(stage) && (!integer_scalar(stage) || !stage %in% 1:3))
    stop("stage must be NULL, 1, 2, or 3.")
  if (!is.logical(return_if_path) || length(return_if_path) != 1L || is.na(return_if_path))
    stop("return_if_path must be TRUE or FALSE.")
  if (!all(c("D1", "D2", "D3") %in% names(data))) stop("Missing treatment columns.")
  D <- as.matrix(data[c("D1", "D2", "D3")])
  if (!is.numeric(D)) stop("Treatment columns must be numeric.")
  if (is.null(nuisance)) nuisance <- fit_iv_nuisance(data)
  required <- c("rows", "alpha", "gamma", "EZ", "logG", "if_alpha", "if_gamma",
                "baseline", "if_Lambda", "risk_if")
  if (!inherits(nuisance, "iv_nuisance") || !all(required %in% names(nuisance)) ||
      !all(vapply(nuisance[c("logG", "if_Lambda", "risk_if")], is.function, logical(1))))
    stop("nuisance must be a complete object returned by fit_iv_nuisance().")
  if (!identical(rows, nuisance$rows)) stop("nuisance data or row order does not match data.")
  if (!identical(dim(nuisance$if_alpha), c(n, 3L)) ||
      !identical(dim(nuisance$if_gamma), c(n, 2L)) || length(nuisance$EZ) != n ||
      any(!is.finite(nuisance$if_alpha)) || any(!is.finite(nuisance$if_gamma)) ||
      any(!is.finite(nuisance$EZ)))
    stop("Invalid nuisance IF dimensions or values.")

  times <- c(0, 0.5, 1, 1.5)
  stages <- if (is.null(stage)) 1:3 else stage
  X <- stats::model.matrix(~ L1 + L2, data)
  Zc <- data$Z - nuisance$EZ
  target <- data$Delta == 1 & data$V_obs == cause
  competing_event <- data$Delta == 1 & data$V_obs != cause
  logG_T <- nuisance$logG(data$T, data)
  risk_weight <- function(t, eligible) {
    a <- as.numeric(data$T >= t)
    competing <- eligible & data$T < t & competing_event
    if (any(competing))
      a[competing] <- exp(nuisance$logG(t, data)[competing] - logG_T[competing])
    a[!eligible] <- 0
    a
  }
  projection_weight <- function(t, eligible) {
    Y1 <- data$T >= t | (data$T < t & competing_event)
    mean(eligible & Y1)
  }
  events <- lapply(stages, function(k) {
    sort(unique(data$T[target & data$T >= times[k] & data$T < times[k + 1]]))
  })
  M <- sum(lengths(events))
  path <- data.frame(time = numeric(M), stage = integer(M), dB = numeric(M),
                     B = numeric(M), std_B = numeric(M), std_B_fixed = numeric(M))
  summaries <- vector("list", length(stages))
  if (return_if_path) epsilon_B_path <- matrix(0, n, M)
  E <- U <- matrix(0, n, length(stages))
  effect <- epsilon_N <- U_N <- numeric(n)
  epsilon_B <- epsilon_B_fixed <- epsilon_B_nuisance <- numeric(n)
  B <- projection_num <- projection_den <- 0
  m <- 0L

  # Event recursion: unchanged B estimate and fixed/nuisance IF components.
  for (s in seq_along(stages)) {
    k <- stages[s]
    active <- seq_len(s)
    eligible <- data$T >= times[k]
    D_history <- D[eligible, stages[active], drop = FALSE]
    if (any(!is.finite(D_history))) stop("Missing treatment in stage ", k, ".")
    Dk <- D[eligible, k]

    for (t in events[[s]]) {
      a <- risk_weight(t, eligible)
      exp_effect <- exp(effect[eligible])
      q <- Zc[eligible] * exp_effect * a[eligible]
      denominator <- sum(q * Dk) / n
      scale <- sum(abs(q * Dk)) / n
      if (!is.finite(denominator) || abs(denominator) <= sqrt(.Machine$double.eps) * scale)
        stop("Zero, near-zero, or non-finite IV denominator at time ", t, ".")
      dN <- as.numeric(target[eligible] & data$T[eligible] == t)
      dB <- (sum(q * dN) / n) / denominator
      r <- dN - Dk * dB
      u <- numeric(n)
      u[eligible] <- (q / denominator) * r
      g <- drop(crossprod(D_history, u[eligible])) / n
      dE <- u + drop(E[, active, drop = FALSE] %*% g)

      alpha_gradient <- drop(crossprod(X[eligible, , drop = FALSE],
                                      exp_effect * a[eligible] * r)) / (n * denominator)
      b <- numeric(n)
      b[eligible] <- exp_effect * Zc[eligible] * r / denominator
      forcing <- -drop(nuisance$if_alpha %*% alpha_gradient) + nuisance$risk_if(t, a, b)
      dU <- forcing + drop(U[, active, drop = FALSE] %*% g)

      # Both IF recursions use states before this event's update.
      E[, s] <- E[, s] + dE
      U[, s] <- U[, s] + dU
      effect[eligible] <- effect[eligible] + Dk * dB
      B <- B + dB
      epsilon_B_fixed <- rowSums(E)
      epsilon_B_nuisance <- rowSums(U)
      epsilon_B <- epsilon_B_fixed + epsilon_B_nuisance
      if (!is.finite(B) || any(!is.finite(effect)) ||
          any(!is.finite(epsilon_B_fixed)) || any(!is.finite(epsilon_B)))
        stop("Non-finite event recursion at time ", t, ".")

      risk <- projection_weight(t, eligible)
      projection_num <- projection_num + risk * dB
      epsilon_N <- epsilon_N + risk * dE
      U_N <- U_N + risk * dU
      m <- m + 1L
      path[m, ] <- list(t, k, dB, B, sqrt(sum(epsilon_B^2) / n^2),
                       sqrt(sum(epsilon_B_fixed^2) / n^2))
      if (return_if_path) epsilon_B_path[, m] <- epsilon_B
    }

    var_B <- sum(epsilon_B^2) / n^2
    var_B_fixed <- sum(epsilon_B_fixed^2) / n^2
    summaries[[s]] <- data.frame(time = times[k + 1], stage = k, B = B,
                                 var_B = var_B, std_B = sqrt(var_B),
                                 var_B_fixed = var_B_fixed, std_B_fixed = sqrt(var_B_fixed))
    dt <- (times[k + 1] - times[k]) / n_grid
    grid <- times[k] + (seq_len(n_grid) - 0.5) * dt
    for (t in grid)
      projection_den <- projection_den + dt * projection_weight(t, eligible)
  }

  # Project increment IFs with fixed observed weights, as in ivsacim.
  if (!is.finite(projection_den) || projection_den <= 0)
    stop("The projection denominator must be finite and positive.")
  beta <- projection_num / projection_den
  epsilon_beta_fixed <- epsilon_N / projection_den
  epsilon_beta_nuisance <- U_N / projection_den
  epsilon_beta <- epsilon_beta_fixed + epsilon_beta_nuisance
  var_beta <- sum(epsilon_beta^2) / n^2
  var_beta_fixed <- sum(epsilon_beta_fixed^2) / n^2
  if (!is.finite(beta) || any(!is.finite(epsilon_beta)) ||
      !is.finite(var_beta) || !is.finite(var_beta_fixed))
    stop("Non-finite projection estimate or variance.")
  result <- list(
    beta = beta, var_beta = var_beta, std = sqrt(var_beta),
    var_beta_fixed = var_beta_fixed, std_fixed = sqrt(var_beta_fixed),
    epsilon_beta = epsilon_beta, epsilon_beta_fixed = epsilon_beta_fixed,
    epsilon_beta_nuisance = epsilon_beta_nuisance,
    epsilon_B = epsilon_B, epsilon_B_fixed = epsilon_B_fixed,
    epsilon_B_nuisance = epsilon_B_nuisance,
    B_summary = do.call(rbind, summaries), path = path, nuisance = nuisance
  )
  if (return_if_path) result$epsilon_B_path <- epsilon_B_path
  result
}

# Three cumulative B(t-) endpoints and pooled beta, with pointwise 95% intervals.
# Reuse the existing recursion and full nuisance-corrected standard errors.
estimation_iv <- function(
    data, cause = 1, n_grid = 10, nuisance = NULL, return_if_path = FALSE) {
  fit <- .estimation_iv_core(
    data, cause = cause, stage = NULL, n_grid = n_grid,
    nuisance = nuisance, return_if_path = return_if_path
  )
  B <- fit$B_summary
  fit$summary <- data.frame(
    cause = cause,
    parameter = c(paste0("B(", B$time, ")"), "beta (pooled)"),
    estimate = c(B$B, fit$beta),
    std = c(B$std_B, fit$std),
    variance = c(B$var_B, fit$var_beta)
  )
  fit$summary$lower <- fit$summary$estimate - 1.96 * fit$summary$std
  fit$summary$upper <- fit$summary$estimate + 1.96 * fit$summary$std
  fit
}
