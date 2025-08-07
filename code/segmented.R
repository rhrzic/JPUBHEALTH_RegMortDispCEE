require(segmented)

segmented_wrapper <- function(y, x = NULL, se=NULL, npsi = 0:3) {
  
  # If x is not provided, default to sequence along y
  if (is.null(x)) {
    n <- length(y)
    x <- 1:n
  }
  
  # Define a restricted range for candidate joinpoints: exclude the first 2 and last 2 years
  lower_bound <- min(x) + 2
  upper_bound <- max(x) - 2
  min_gap <- 2  # minimum gap between joinpoints (in the same units as x)
  
  # If standard errors are not provided, use an unweighted model
  if (is.null(se)) {
    initial_model <- glm(y ~ x, family = gaussian)
  } else {
    initial_model <- glm(y ~ x, weights = 1/(se^2), family = gaussian)
  }
  
  ris.ok <- list()
  successful_npsi <- numeric(0)
  
  cat("Running segmented models ...\n")
    
  # Loop over each candidate number of breakpoints
    for (i in npsi) {
      modelCandidate <- tryCatch(
        withCallingHandlers(
          segmented(initial_model, npsi = i,
                    control = seg.control(n.boot = 100, seed = 9999, it.max = 100)),
          warning = function(w) {
            # If any warning occurs, trigger a stop that exits the loop.
            stop("segmented warning: ", w$message, call. = FALSE)
          }
        ),
        error = function(e) { NULL }
      )
      
    # Only process successful segmented models
    if (!is.null(modelCandidate) && inherits(modelCandidate, "segmented")) {
      
      # Retrieve estimated breakpoints as a vector
      psi_est <- NULL
      if (!is.null(modelCandidate$psi)) {
        if (is.matrix(modelCandidate$psi)) {
          psi_est <- sort(as.vector(modelCandidate$psi[, "Est."]))
        } else {
          psi_est <- sort(as.vector(modelCandidate$psi))
        }
      }
      
      # Check if any joinpoint is too close to the boundaries
      if (!is.null(psi_est)) {
        if(any(psi_est < lower_bound) || any(psi_est > upper_bound)) {
          next  # Skip this candidate
        }
        # If more than one joinpoint, enforce a minimum gap between them
        if (length(psi_est) > 1 && any(diff(psi_est) < min_gap)) {
          next  # Skip this candidate
        }
      }
      
      # Try extracting the confidence intervals using confint.segmented
      candidate_ci <- try(confint.segmented(modelCandidate), silent = TRUE)
      
      # If confint fails or any CI is extremely large, skip this candidate
      if (inherits(candidate_ci, "try-error")) {
        next
      }
      
      # Ensure candidate_ci is a matrix (it might be a vector if only one breakpoint)
      if (!is.matrix(candidate_ci)) {
        candidate_ci <- matrix(candidate_ci, nrow = 1, dimnames = list(NULL, names(candidate_ci)))
      }
      
      covmat <- try(vcov(modelCandidate), silent = TRUE)
      if (inherits(covmat, "try-error") || any(is.na(diag(covmat)))) {
        next
      }
      
      # Check if any of the CI bounds are larger than 1e09 in absolute value
      if (any(candidate_ci[, "CI(95%).up"] > 2030, na.rm = TRUE) ||
          any(candidate_ci[, "CI(95%).low"] < 1990, na.rm = TRUE)) {
        next
      }
      
      # Candidate passes all checks; store it
      ris.ok[[length(ris.ok) + 1]] <- modelCandidate
      successful_npsi <- c(successful_npsi, i)
    }
  }

  
  if (length(ris.ok) == 0) {
    cat("No successful segmented models found; returning linear model.\n")
    return(initial_model)
  }  
  
  # Use BIC to select the best model from the successful fits
  r <- sapply(ris.ok, BIC)
  id.ok <- which.min(r)
  npsi_selected <- successful_npsi[id.ok]
  
  o <- ris.ok[[id.ok]]
  cat(paste(npsi_selected, "breakpoints selected by BIC\n"))
  
  names(r) <- successful_npsi
  o$bic <- r
  class(o) <- c("segmented", "glm", "lm")
  return(o)
}

# Define a function to fit the segmented model with bootstrapping for one subgroup of data
fit_segmented_boot <- function(data, response, se) {
 
  data <- data[order(data$Year), ]
  y    <- data[[response]]
  w    <- 1 / (data[[se]]^2)   # weights = 1 / Var
  x    <- data$Year
  
  # Fit the segmented model (trying 0, 1, 2, or 3 breakpoints)
  seg_model <- segmented_wrapper(y = y, x = x, se = data[[se]], npsi = 0:3)
  
  # Extract fitted values and record the chosen model type
  data$fitted <- fitted(seg_model)
  
  if (is.null(seg_model$psi)) {
    num_breakpoints <- 0
  } else if (is.matrix(seg_model$psi)) {
    num_breakpoints <- nrow(seg_model$psi)
  } else {
    num_breakpoints <- length(seg_model$psi)
  }
  data$model_type <- paste(num_breakpoints, "breakpoints selected by BIC")
  
  # Extract breakpoint confidence intervals if breakpoints are present
  
  if(num_breakpoints == 0) {
    data$bp1_est <- NA
    data$bp1_lower <- NA
    data$bp1_upper <- NA
    data$bp2_est <- NA
    data$bp2_lower <- NA
    data$bp2_upper <- NA
    data$bp3_est <- NA
    data$bp3_lower <- NA
    data$bp3_upper <- NA
  } else if(num_breakpoints == 1) {
    ci <- confint.segmented(seg_model)
    if(is.null(dim(ci))) {
      data$bp1_est <- ci["Est."]
      data$bp1_lower <- ci["CI(95%).low"]
      data$bp1_upper <- ci["CI(95%).up"]
    } else {
      data$bp1_est <- round(ci[1, "Est."], 1)
      data$bp1_lower <- round(ci[1, "CI(95%).low"], 1)
      data$bp1_upper <- round(ci[1, "CI(95%).up"], 1)
    }
    data$bp2_est <- NA
    data$bp2_lower <- NA
    data$bp2_upper <- NA
    data$bp2_est <- NA
    data$bp2_lower <- NA
    data$bp2_upper <- NA
    data$bp3_est <- NA
    data$bp3_lower <- NA
    data$bp3_upper <- NA
  } else if(num_breakpoints == 2) {
    ci <- confint.segmented(seg_model)
    data$bp1_est <- round(ci[1, "Est."], 1)
    data$bp1_lower <- round(ci[1, "CI(95%).low"], 1)
    data$bp1_upper <- round(ci[1, "CI(95%).up"], 1)
    data$bp2_est <- round(ci[2, "Est."], 1)
    data$bp2_lower <- round(ci[2, "CI(95%).low"], 1)
    data$bp2_upper <- round(ci[2, "CI(95%).up"], 1)
    data$bp3_est <- NA
    data$bp3_lower <- NA
    data$bp3_upper <- NA
  } else if(num_breakpoints == 3) {
    ci <- confint.segmented(seg_model)
    data$bp1_est <- round(ci[1, "Est."], 1)
    data$bp1_lower <- round(ci[1, "CI(95%).low"], 1)
    data$bp1_upper <- round(ci[1, "CI(95%).up"], 1)
    data$bp2_est <- round(ci[2, "Est."], 1)
    data$bp2_lower <- round(ci[2, "CI(95%).low"], 1)
    data$bp2_upper <- round(ci[2, "CI(95%).up"], 1)
    data$bp3_est <- round(ci[3, "Est."], 1)
    data$bp3_lower <- round(ci[3, "CI(95%).low"], 1)
    data$bp3_upper <- round(ci[3, "CI(95%).up"], 1)
    } 
  
  # ---- Extract APC results ----
  # Attempt to extract slopes for each segment using slope()
  slopes_list <- try(slope(seg_model, parm = "x"), silent = TRUE)
  
  
  # If slopes_list is a list, convert it to a matrix
  if (!inherits(slopes_list, "try-error") && is.list(slopes_list)) {
    slopes_mat <- do.call(rbind, slopes_list)
  } else {
    slopes_mat <- slopes_list
  }

  
  # Create APC columns for up to 6 segments (0 breakpoints = 1 segment, up to 5 breakpoints = 6 segments)
  for(i in 1:4) {
    data[[paste0("AAPC", i)]] <- NA
    data[[paste0("AAPC", i, "_lower")]] <- NA
    data[[paste0("AAPC", i, "_upper")]] <- NA
  }
  
  tol <- 1e-6
  
  
  if (!inherits(slopes_mat, "try-error")) {
    # Ensure slopes_mat is a matrix
    if (!is.matrix(slopes_mat)) {
      slopes_mat <- matrix(slopes_mat, nrow = 1, 
                           dimnames = list(NULL, names(slopes_mat)))
    }
    num_segments <- nrow(slopes_mat)  # This equals num_breakpoints + 1
    for (i in 1:num_segments) {
      est <- slopes_mat[i, "Est."]
      low <- slopes_mat[i, "CI(95%).l"]
      up <- slopes_mat[i, "CI(95%).u"]
      
      # Compute APC: (exp(slope) - 1)*100
      
      # Compute APC while avoiding numerical issues when slope is near zero
      apc_val <- if (is.na(est)) {
        NA
      } else if (abs(est) < tol) {
        0
      } else {
        (exp(est) - 1) * 100
      }
      apc_low <- if (is.na(low)) {
        NA
      } else if (abs(low) < tol) {
        0
      } else {
        (exp(low) - 1) * 100
      }
      apc_up  <- if (is.na(up)) {
        NA
      } else if (abs(up) < tol) {
        0
      } else {
        (exp(up) - 1) * 100
      }
      
      # Replace any infinite or NaN values with NA
      if(is.infinite(apc_val) || is.nan(apc_val)) { apc_val <- NA }
      if(is.infinite(apc_low) || is.nan(apc_low)) { apc_low <- NA }
      if(is.infinite(apc_up)  || is.nan(apc_up))  { apc_up <- NA }
      
      data[[paste0("AAPC", i)]] <- round(apc_val, 2)
      data[[paste0("AAPC", i, "_lower")]] <- round(apc_low, 2)
      data[[paste0("AAPC", i, "_upper")]] <- round(apc_up, 2)
    }
  }
  
  #    a) time span
  yr_min  <- min(data$Year)
  yr_max  <- max(data$Year)
  #    b) predict on the log‐scale at endpoints
  newdat   <- data.frame(x = c(yr_min, yr_max))
  log_preds <- predict(seg_model, newdata = newdat, type = "response")
  #    c) total change & CI
  delta_log   <- log_preds[2] - log_preds[1]
  data$total_pct_change <- (exp(delta_log) - 1) * 100
  
  #    d) approximate CI for delta_log using delta method on predict()
  #       get var–cov of predictions:
  Vpred <- vcov(seg_model)
  
  #       build contrast vector C = [−1 at row1, +1 at row2] over the intercept+slopes
  #       but simpler: use predict(..., se.fit=TRUE) to get se(delta_log)
  se_preds <- predict(seg_model, newdata = newdat, se.fit = TRUE)$se.fit
  se_delta <- sqrt(se_preds[1]^2 + se_preds[2]^2)  # assuming independence of fits
  ci_low_log  <- delta_log - 1.96 * se_delta
  ci_high_log <- delta_log + 1.96 * se_delta
  
  data$total_pct_change_lower <- (exp(ci_low_log)  - 1) * 100
  data$total_pct_change_upper <- (exp(ci_high_log) - 1) * 100
  
  duration = yr_max-yr_min
  
  annual_beta <- delta_log / duration
  ci_low_beta <- ci_low_log  / duration
  ci_high_beta<- ci_high_log / duration
  
  data$AAPC_total        <- (exp(annual_beta)  - 1) * 100
  data$AAPC_total_lower  <- (exp(ci_low_beta)  - 1) * 100
  data$AAPC_total_upper  <- (exp(ci_high_beta) - 1) * 100
  
  return(data)
}
