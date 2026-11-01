# SI simulation: data generation, nuisance models, and the estimator used in the study.

## SI simulation with stagewise event realization

simulation_si <- function(n = 10000, seed = 1, event_rate = 0.60) {
  set.seed(seed)
  
  times <- c(0, 0.5, 1, 1.5)
  tau <- 1.5
  beta <- c(-0.25, -0.125) * event_rate
  b <- c(-2.4245972, -1.5648315, -0.6588387, 0.6722480)[
    match(event_rate, c(.07, .15, .30, .60))
  ]
  
  # CIF under treatment history D1,...,Dk
  cif <- function(t, L1, L2, D, k = 3) {
    rate <- c(0.025, 0.035)
    weight <- c(0.6, 0.4)
    m <- plogis(b + 0.05 * L1 + 0.05 * (L2 - 0.5))
    
    F <- weight * m * (1 - exp(-rate * t)) / (1 - exp(-rate * tau))
    f <- weight * m * rate * exp(-rate * t) / (1 - exp(-rate * tau))
    if (k == 0) return(list(F = F, f = f))
    if (t < times[k]) return(cif(t, L1, L2, D, k - 1))
    
    prev <- cif(t, L1, L2, D, k - 1)
    bound <- cif(times[k], L1, L2, D, k - 1)$F
    ell <- min(t, times[k + 1]) - times[k]
    e <- exp(-beta * D[k] * ell)
    de <- if (t < times[k + 1]) beta * D[k] * e else c(0, 0)
    
    list(
      F = prev$F + (1 - rev(bound) - prev$F) * (1 - e),
      f = prev$f * e + (1 - rev(bound) - prev$F) * de
    )
  }
  
  # L and treatment
  L1 <- qnorm(runif(n, pnorm(-3), pnorm(3)))
  L2 <- rbinom(n, 1, 0.5)
  D <- matrix(0L, n, 3)
  
  D[, 1] <- rbinom(n, 1, plogis(-1 + 0.3 * L1 + 0.3 * (L2 - 0.5)))
  for (k in 2:3) {
    pD <- plogis(
      -1.7 + 2 * D[, k - 1] +
        0.6 * rowMeans(D[, 1:(k - 1), drop = FALSE]) +
        0.3 * L1 + 0.3 * (L2 - 0.5)
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
      F0 <- cif(t0, L1[i], L2[i], D[i, ], k)$F
      F1 <- cif(t1, L1[i], L2[i], D[i, ], k)$F
      
      S0 <- 1 - sum(F0)
      p <- (F1 - F0) / S0   # P(event of each cause in stage k | T_tilde >= t_k)
      
      if (U_E[i, k] > sum(p)) next
      j <- if (U_E[i, k] <= p[1]) 1L else 2L
      
      target <- F0[j] + U_V[i, k] * (F1[j] - F0[j])
      event_time[i] <- uniroot(
        function(t) cif(t, L1[i], L2[i], D[i, ], k)$F[j] - target,
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
    id = 1:n, L1, L2, D1 = D[, 1], D2 = D[, 2], D3 = D[, 3],
    event_time, event_type, C_star, C, T, Delta, V_obs
  )
  observed <- full[c("id", "L1", "L2", "D1", "D2", "D3", "T", "Delta", "V_obs")]
  for (k in 1:3) observed[observed$T < times[k], paste0("D", k)] <- NA_integer_
  
  list(
    observed = observed,
    full_data = full,
    parameters = list(times = times, tau = tau, beta = beta, b = b,
                      censor_lambda0 = lambda0, censor_shape = 2)
  )
}

## Stagewise SI dA estimates with treatment and Cox censoring nuisance IFs.
## Source this file alone. stage = NULL returns three stage estimates and their pooled estimate.
## All IFs use the full sample scale: variance = sum(IF^2) / n^2.

.si_db_stages <- function(stage) {
  if (is.null(stage)) return(1:3)
  if (!is.numeric(stage) || length(stage) != 1L || is.na(stage) || !stage %in% 1:3)
    stop("stage must be NULL, 1, 2, or 3.")
  as.integer(stage)
}

.si_db_rows <- function(data) {
  outcomes <- c("L1", "L2", "T", "Delta", "V_obs")
  columns <- c(outcomes, "D1", "D2", "D3")
  if (!is.data.frame(data) || nrow(data) < 2 || !all(columns %in% names(data)))
    stop("data must have at least two rows and columns: ", paste(columns, collapse = ", "))
  if (!all(vapply(data[columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(data[outcomes]))))
    stop("Data columns must be numeric; covariates and outcomes must be finite.")
  if (any(data$T < 0) || any(!data$Delta %in% 0:1) ||
      any(data$V_obs < 0 | data$V_obs != floor(data$V_obs)) ||
      any(data$Delta == 1 & data$V_obs == 0))
    stop("Invalid event times, event indicators, or observed causes.")
  data[c(columns, intersect("id", names(data)))]
}

.si_db_inverse <- function(information, label) {
  if (any(!is.finite(information))) stop(label, " is non-finite.")
  tryCatch(solve(information), error = function(e) {
    stop(label, " is singular or numerically singular.", call. = FALSE)
  })
}

## Cox/Breslow censoring nuisance calculation.
.si_db_censoring <- function(data) {
  n <- nrow(data)
  x <- data
  x$censor_event <- as.integer(data$Delta == 0 & data$T < 1.5)
  if (!any(x$censor_event)) stop("No censoring events are available to fit the Cox model.")
  fit <- tryCatch(
    survival::coxph(survival::Surv(T, censor_event) ~ L1 + L2,
                   data = x, ties = "breslow", x = TRUE, y = TRUE),
    warning = function(w) stop("Censoring Cox fit: ", conditionMessage(w), call. = FALSE)
  )
  gamma <- stats::coef(fit)
  if (any(!is.finite(gamma))) stop("The censoring Cox fit has undefined coefficients.")
  L <- fit$x
  Tc <- fit$y[, 1]  # Preserve Cox's near-tie handling in the fitted risk sets.
  censored <- fit$y[, 2] == 1
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

  # Zero-reference Breslow baseline and its individual IFs.
  Lambda <- c(0, cumsum(dLambda))
  F <- c(0, cumsum(dLambda / S0))
  A <- rbind(c(0, 0), cbind(cumsum(Lbar[, 1] * dLambda),
                           cumsum(Lbar[, 2] * dLambda)))
  index <- function(t) findInterval(t, censor_times) + 1L
  ic <- index(Tc)
  score <- -rr * (L * Lambda[ic] - A[ic, , drop = FALSE])
  score[censored, ] <- score[censored, , drop = FALSE] +
    L[censored, , drop = FALSE] - Lbar[match(Tc[censored], censor_times), , drop = FALSE]
  if_gamma <- score %*% .si_db_inverse(information, "Censoring Cox information")
  if (any(!is.finite(if_gamma))) stop("Non-finite censoring influence functions.")
  event_coefficient <- numeric(n)
  event_coefficient[censored] <- 1 / S0[match(Tc[censored], censor_times)]
  if_Lambda <- function(t) {
    event_coefficient * (Tc <= t) - rr * F[index(pmin(Tc, t))] -
      drop(if_gamma %*% A[index(t), ])
  }
  logG <- function(t, data) {
    covariates <- as.matrix(data[c("L1", "L2")])
    -Lambda[index(t)] * exp(drop(covariates %*% gamma))
  }

  T <- data$T
  order_T <- order(T)
  lower <- findInterval(Tc, T[order_T], left.open = TRUE)
  upper <- findInterval(Tc, T[order_T])
  iT <- index(T)
  Lambda_T <- Lambda[iT]
  F_T <- F[iT]
  F_c <- F[ic]
  A_T <- A[iT, , drop = FALSE]

  # delta_i mean(b_r * a_r(t)); b is fixed, a includes stage eligibility.
  # Prefix sums avoid constructing an n-by-n matrix of baseline IFs.
  risk_if <- function(t, a, b = rep(1, n)) {
    v <- b * a * rr
    v[T >= t] <- 0
    if (!any(v != 0)) return(numeric(n))
    total <- sum(v)
    prefix <- c(0, cumsum(v[order_T]))
    prefix_F <- c(0, cumsum((v * F_T)[order_T]))
    weighted_IF <- event_coefficient * (total - prefix[lower + 1L]) -
      rr * (prefix_F[upper + 1L] + F_c * (total - prefix[upper + 1L])) -
      drop(if_gamma %*% colSums(A_T * v))
    gamma_gradient <- colSums(L * (v * (Lambda[index(t)] - Lambda_T)))
    -(total * if_Lambda(t) - weighted_IF + drop(if_gamma %*% gamma_gradient)) / n
  }
  list(gamma = gamma, if_gamma = if_gamma, logG = logG,
       baseline = data.frame(time = censor_times, dLambda = dLambda,
                             Lambda = Lambda[-1], S0 = S0),
       if_Lambda = if_Lambda, risk_if = risk_if)
}

## Fit once and reuse across causes or any of the fitted stages.
fit_si_nuisance <- function(data, stage = NULL) {
  rows <- .si_db_rows(data)
  stages <- .si_db_stages(stage)
  n <- nrow(data)
  times <- c(0, 0.5, 1, 1.5)
  formulas <- list(D1 ~ L1 + L2, D2 ~ D1 + L1 + L2,
                   D3 ~ D2 + I((D1 + D2) / 2) + L1 + L2)
  treatment <- setNames(vector("list", length(stages)), as.character(stages))
  for (k in stages) {
    keep <- data$T >= times[k]
    if (!any(keep)) stop("No eligible subjects in stage ", k, ".")
    history <- as.matrix(data[keep, paste0("D", seq_len(k)), drop = FALSE])
    if (any(!is.finite(history)) || any(!history %in% 0:1))
      stop("Stage ", k, " requires finite binary treatment and history for eligible subjects.")
    fit <- tryCatch(
      stats::glm(formulas[[k]], stats::binomial(), data = data[keep, , drop = FALSE]),
      warning = function(w) stop("Treatment model in stage ", k, ": ",
                                 conditionMessage(w), call. = FALSE)
    )
    X <- stats::model.matrix(fit)
    alpha <- stats::coef(fit)
    if (!fit$converged || fit$rank != ncol(X) || any(!is.finite(alpha)))
      stop("Treatment model did not converge or is singular in stage ", k, ".")
    eta <- drop(X %*% alpha)
    p <- stats::plogis(eta)
    information <- crossprod(X, X * (p * (1 - p))) / n
    if_alpha <- matrix(0, n, ncol(X), dimnames = list(NULL, colnames(X)))
    if_alpha[keep, ] <- (X * (data[[paste0("D", k)]][keep] - p)) %*%
      .si_db_inverse(information, paste("Treatment information in stage", k))
    if (any(!is.finite(eta)) || any(!is.finite(if_alpha)))
      stop("Non-finite treatment nuisance in stage ", k, ".")
    treatment[[as.character(k)]] <- list(alpha = alpha, X = X, eta = eta,
                                        if_alpha = if_alpha)
  }
  # Keep this object immutable: censoring closures share fitted quantities.
  structure(list(rows = rows, treatment = treatment, censor = .si_db_censoring(data)),
            class = "si_nuisance")
}

## Projection weights are observed eligible * Y1 (no censoring adjustment).
## Both projection weights and their integral are fixed in beta's IF.
## summary reports A(time-); path columns are stage-local, never cumulative across stages.
## pooled combines the selected stages; selecting one stage reproduces its beta and IF.
## beta_only skips A diagnostics and returns only summary and pooled, for batch runs.
.estimation_si_core <- function(
    data, cause = 1, stage = NULL, n_grid = 10,
    nuisance = NULL, return_if_path = FALSE, beta_only = FALSE) {
  rows <- .si_db_rows(data)
  stages <- .si_db_stages(stage)
  n <- nrow(data)
  integer_scalar <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x)
  }
  if (!integer_scalar(cause) || cause < 1) stop("cause must be a positive integer.")
  if (!integer_scalar(n_grid) || n_grid < 1) stop("n_grid must be a positive integer.")
  if (!is.logical(return_if_path) || length(return_if_path) != 1L || is.na(return_if_path))
    stop("return_if_path must be TRUE or FALSE.")
  if (!is.logical(beta_only) || length(beta_only) != 1L || is.na(beta_only))
    stop("beta_only must be TRUE or FALSE.")
  if (beta_only && return_if_path) stop("beta_only cannot be combined with return_if_path.")
  if (is.null(nuisance)) nuisance <- fit_si_nuisance(data, stage)
  if (!inherits(nuisance, "si_nuisance") ||
      !all(c("rows", "treatment", "censor") %in% names(nuisance)) ||
      !all(as.character(stages) %in% names(nuisance$treatment)))
    stop("nuisance must be returned by fit_si_nuisance() and contain all requested stages.")
  if (!identical(rows, nuisance$rows)) stop("nuisance data or row order does not match data.")
  censor <- nuisance$censor
  functions <- c("logG", "if_Lambda", "risk_if")
  if (!all(functions %in% names(censor)) ||
      !all(vapply(censor[functions], is.function, logical(1))))
    stop("Invalid censoring nuisance functions.")

  times <- c(0, 0.5, 1, 1.5)
  target <- data$Delta == 1 & data$V_obs == cause
  competing_event <- data$Delta == 1 & data$V_obs != cause
  logG_T <- censor$logG(data$T, data)
  risk_weight <- function(t, eligible) {
    a <- as.numeric(data$T >= t)
    competing <- eligible & data$T < t & competing_event
    if (any(competing))
      a[competing] <- exp(censor$logG(t, data)[competing] - logG_T[competing])
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
  if (!beta_only)
    path <- data.frame(time = numeric(M), stage = integer(M), dA = numeric(M),
                       A = numeric(M), std_A = numeric(M), std_A_fixed = numeric(M))
  summaries <- vector("list", length(stages))
  projection_den_by_stage <- setNames(numeric(length(stages)), paste0("stage", stages))
  epsilon_beta_fixed <- epsilon_beta_nuisance <-
    matrix(0, n, length(stages), dimnames = list(NULL, paste0("stage", stages)))
  if (!beta_only) epsilon_A_fixed <- epsilon_A_nuisance <- epsilon_beta_fixed
  if (return_if_path) epsilon_A_path <- matrix(0, n, M)
  j <- 0L

  for (s in seq_along(stages)) {
    k <- stages[s]
    eligible <- data$T >= times[k]
    Dk <- data[[paste0("D", k)]][eligible]
    tr <- nuisance$treatment[[as.character(k)]]
    if (!is.matrix(tr$X) || nrow(tr$X) != sum(eligible) ||
        length(tr$eta) != sum(eligible) ||
        !identical(dim(tr$if_alpha), c(n, ncol(tr$X))) ||
        any(!is.finite(c(tr$X, tr$eta, tr$if_alpha))))
      stop("Invalid treatment nuisance dimensions or values in stage ", k, ".")
    A <- projection_num <- 0
    E <- U <- numerator_fixed <- numerator_nuisance <- numeric(n)

    for (t in events[[s]]) {
      a <- risk_weight(t, eligible)
      m <- stats::plogis(tr$eta - A)
      V <- m * (1 - m)
      W <- Dk - m
      q <- W * a[eligible]
      denominator <- sum(q * Dk) / n
      scale <- sum(abs(q * Dk)) / n
      if (!is.finite(denominator) || abs(denominator) <= sqrt(.Machine$double.eps) * scale)
        stop("Zero, near-zero, or non-finite SI denominator in stage ", k, " at time ", t, ".")
      dN <- as.numeric(target[eligible] & data$T[eligible] == t)
      dA <- sum(q * dN) / (n * denominator)
      r <- dN - Dk * dA
      dR <- a[eligible] * r
      direct <- numeric(n)
      direct[eligible] <- W * dR / denominator
      feedback <- sum(V * dR) / (n * denominator)
      dE <- direct + feedback * E

      alpha_gradient <- drop(crossprod(tr$X, V * dR)) / (n * denominator)
      b <- numeric(n)
      b[eligible] <- W * r / denominator
      forcing <- -drop(tr$if_alpha %*% alpha_gradient) + censor$risk_if(t, a, b)
      dU <- forcing + feedback * U

      # Both IF increments use pre-update states, including at tied event times.
      E <- E + dE
      U <- U + dU
      A <- A + dA
      if (!is.finite(A) || any(!is.finite(E)) || any(!is.finite(U)))
        stop("Non-finite SI recursion in stage ", k, " at time ", t, ".")
      weight <- projection_weight(t, eligible)
      projection_num <- projection_num + weight * dA
      numerator_fixed <- numerator_fixed + weight * dE
      numerator_nuisance <- numerator_nuisance + weight * dU
      if (!beta_only) {
        j <- j + 1L
        path[j, ] <- list(t, k, dA, A, sqrt(sum((E + U)^2) / n^2),
                          sqrt(sum(E^2) / n^2))
        if (return_if_path) epsilon_A_path[, j] <- E + U
      }
    }

    dt <- (times[k + 1] - times[k]) / n_grid
    grid <- times[k] + (seq_len(n_grid) - 0.5) * dt
    projection_den <- dt * sum(vapply(grid, projection_weight, numeric(1), eligible = eligible))
    if (!is.finite(projection_den) || projection_den <= 0)
      stop("The projection denominator must be finite and positive in stage ", k, ".")
    projection_den_by_stage[s] <- projection_den
    beta <- projection_num / projection_den
    ef <- numerator_fixed / projection_den
    eu <- numerator_nuisance / projection_den
    var_beta <- sum((ef + eu)^2) / n^2
    var_beta_fixed <- sum(ef^2) / n^2
    if (any(!is.finite(c(beta, ef, eu, var_beta, var_beta_fixed))))
      stop("Non-finite SI projection estimate or variance in stage ", k, ".")
    epsilon_beta_fixed[, s] <- ef
    epsilon_beta_nuisance[, s] <- eu
    # No target events yield zero A, beta and empirical IFs; this is not a coverage guarantee.
    summaries[[s]] <- data.frame(
      cause = cause, stage = k, time = times[k + 1], beta = beta,
      var_beta = var_beta, std = sqrt(var_beta), var_beta_fixed = var_beta_fixed,
      std_fixed = sqrt(var_beta_fixed)
    )
    if (!beta_only) {
      epsilon_A_fixed[, s] <- E
      epsilon_A_nuisance[, s] <- U
      summaries[[s]]$A <- A
      summaries[[s]]$std_A <- sqrt(sum((E + U)^2) / n^2)
      summaries[[s]]$std_A_fixed <- sqrt(sum(E^2) / n^2)
    }
  }
  # Combine each subject's stage IFs before taking the variance, retaining cross-stage covariance.
  summary <- do.call(rbind, summaries)
  pooled_weights <- projection_den_by_stage / sum(projection_den_by_stage)
  pooled_fixed <- drop(epsilon_beta_fixed %*% pooled_weights)
  pooled_nuisance <- drop(epsilon_beta_nuisance %*% pooled_weights)
  pooled_if <- pooled_fixed + pooled_nuisance
  pooled_var <- sum(pooled_if^2) / n^2
  pooled_var_fixed <- sum(pooled_fixed^2) / n^2
  pooled <- data.frame(
    cause = cause, beta = sum(pooled_weights * summary$beta),
    var_beta = pooled_var, std = sqrt(pooled_var),
    var_beta_fixed = pooled_var_fixed, std_fixed = sqrt(pooled_var_fixed)
  )
  if (beta_only) return(list(summary = summary, pooled = pooled))
  result <- list(
    summary = summary, pooled = pooled, path = path,
    projection_den = projection_den_by_stage, pooled_weights = pooled_weights,
    epsilon_beta_pooled = pooled_if,
    epsilon_beta_pooled_fixed = pooled_fixed, epsilon_beta_pooled_nuisance = pooled_nuisance,
    epsilon_beta = epsilon_beta_fixed + epsilon_beta_nuisance,
    epsilon_beta_fixed = epsilon_beta_fixed, epsilon_beta_nuisance = epsilon_beta_nuisance,
    epsilon_A = epsilon_A_fixed + epsilon_A_nuisance,
    epsilon_A_fixed = epsilon_A_fixed, epsilon_A_nuisance = epsilon_A_nuisance,
    nuisance = nuisance
  )
  if (return_if_path) result$epsilon_A_path <- epsilon_A_path
  result
}

# Cumulative B(t-) at 0.5, 1, and 1.5, retaining all original stage/pooled outputs.
estimation_si <- function(data, cause = 1, n_grid = 10, nuisance = NULL) {
  fit <- .estimation_si_core(
    data, cause = cause, stage = NULL, n_grid = n_grid, nuisance = nuisance,
    return_if_path = FALSE, beta_only = FALSE
  )
  n <- nrow(data)
  psi <- fit$epsilon_A
  # Sum full IFs within each subject before squaring, preserving stage covariance.
  for (k in 2:3) psi[, k] <- psi[, k - 1] + fit$epsilon_A[, k]
  colnames(psi) <- paste0("B(", fit$summary$time, ")")
  asymptotic_variance <- unname(colMeans(psi^2))
  var_B <- asymptotic_variance / n
  B <- cumsum(fit$summary$A)
  std_B <- sqrt(var_B)

  fit$B_summary <- data.frame(
    cause = cause, time = fit$summary$time, B = B,
    asymptotic_variance = asymptotic_variance, var_B = var_B, std_B = std_B,
    lower = B - 1.96 * std_B, upper = B + 1.96 * std_B
  )
  fit$epsilon_B <- psi
  fit
}
