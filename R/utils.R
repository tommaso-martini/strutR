#' Pretty-print a Struthers fronts list (debug utility)
#'
#' @description
#' Prints the current set of wetting fronts as a compact table with indices,
#' water contents \eqn{\theta} and depths \eqn{x} (m). Intended for interactive
#' tracing and debugging.
#'
#' @param fronts list of fronts in the form \code{list(list(theta=..., x=...), ...)}.
#' @param label optional character string printed as a header.
#'
#' @return
#' Invisibly returns \code{NULL}. Side effect: prints to the console.
#'
#' @examples
#' pretty_print_fronts(list(list(theta = 0.20, x = 1.2), list(theta = 0.30, x = 0.6)),
#'   label = "Current fronts"
#' )
#' @export
pretty_print_fronts <- function(fronts, label = NULL) {
  if (!is.null(label)) cat("=== ", label, " ===\n", sep = "")
  if (!length(fronts)) {
    cat("[no fronts]\n")
    return(invisible(NULL))
  }
  df <- data.frame(
    i = seq_along(fronts),
    theta = vapply(fronts, function(f) f$theta, numeric(1)),
    x = vapply(fronts, function(f) f$x, numeric(1))
  )
  print(df, row.names = FALSE)
}

#' Brooks–Corey unsaturated hydraulic conductivity K(θ)
#'
#' @description
#' Computes the Brooks–Corey conductivity \eqn{K(\theta)} with clamping of
#' \eqn{\theta} to \eqn{[\theta_r,\theta_s]}:
#' \deqn{K(\theta) = K_s \left(\frac{\theta - \theta_r}{\theta_s - \theta_r}\right)^{1/\beta}.}
#'
#' @param theta numeric vector; volumetric water content (dimensionless).
#' @param theta_s,theta_r numeric scalars; saturated and residual water contents (–).
#' @param beta numeric scalar; Brooks–Corey shape parameter (\eqn{1/\beta} is the exponent).
#' @param Ks numeric scalar; saturated conductivity [m/day].
#'
#' @return numeric vector of conductivities [m/day], same length as \code{theta}.
#'
#' @examples
#' Kfun(
#'   theta = c(0.10, 0.20, 0.30),
#'   theta_s = 0.45, theta_r = 0.05, beta = 1 / 8, Ks = 1.0
#' )
#' @export
Kfun <- function(theta, theta_s, theta_r, beta, Ks) {
  theta_c <- pmin(pmax(theta, theta_r), theta_s)
  Ks * ((theta_c - theta_r) / (theta_s - theta_r))^(1 / beta)
}

#' Recursively remove NULL/empty elements from a nested list
#'
#' @description
#' Cleans a nested list by removing \code{NULL} elements and empty sub-lists.
#' Useful for housekeeping of front lists after merges.
#'
#' @param x an arbitrary R object (typically a list).
#'
#' @return A list with empty elements pruned. Non-list inputs are returned unchanged.
#'
#' @examples
#' clean_list(list(1, NULL, list(), list(a = NULL, b = 2)))
#' @keywords internal
#' @noRd
clean_list <- function(x) {
  if (!is.list(x)) {
    return(x)
  }
  x <- lapply(x, clean_list)
  x <- x[!sapply(x, function(y) is.null(y) || (is.list(y) && length(y) == 0))]
  x
}

#' Semi-implicit Euler update for interior fronts (Struthers core)
#'
#' @description
#' Computes a robust semi-implicit increment \eqn{\mathrm{d}\theta/\mathrm{d}t}
#' for an interior segment \eqn{i \ge 2} using a mass-balance–motivated denominator.
#' This helper is designed to be called sequentially from deep to shallow layers.
#'
#' @details
#' The update uses the storage at the beginning of the step and the net supply
#' \eqn{q_i} over \eqn{\Delta t}, with a corrective term driven by conductivity
#' differences \eqn{K_i - K_{i-1}} and the newly updated \eqn{\theta_{i-1}^{new}}.
#' The denominator is safeguarded against near-zero values to prevent blow-ups.
#'
#' @param dt numeric; sub-step duration [day].
#' @param t numeric; kept for backward compatibility (not used in balance).
#' @param old_xs numeric; thickness of the segment i (m) relative to its lower bound.
#' @param old_thetas_i numeric; previous \eqn{\theta_i}.
#' @param old_thetas_i_meno_1 numeric; previous \eqn{\theta_{i-1}}.
#' @param q_vec numeric; net flux entering segment i [m/day].
#' @param Kvals_i,Kvals_i_meno_1 numeric; conductivities at i and i-1 [m/day].
#' @param dtheta_dt_i_meno_1 numeric; already-computed rate for i-1.
#' @param new_thetas_i_meno_1 numeric; updated \eqn{\theta_{i-1}^{new}}.
#' @param debug logical; if TRUE, prints detailed diagnostics.
#'
#' @return numeric; \eqn{\mathrm{d}\theta_i/\mathrm{d}t}.
#'
#' @examples
#' f_euler(
#'   dt = 1 / 24, old_xs = 0.2, old_thetas_i = 0.25, old_thetas_i_meno_1 = 0.20,
#'   q_vec = 0.01, Kvals_i = 0.5, Kvals_i_meno_1 = 0.3,
#'   dtheta_dt_i_meno_1 = 0, new_thetas_i_meno_1 = 0.205
#' )
#' @keywords internal
#' @noRd
f_euler <- function(dt,
                    t = 1,
                    old_xs,
                    old_thetas_i,
                    old_thetas_i_meno_1,
                    q_vec,
                    Kvals_i,
                    Kvals_i_meno_1,
                    dtheta_dt_i_meno_1,
                    new_thetas_i_meno_1,
                    debug = FALSE) {
  # Initial storage of segment i with respect to its lower bound (non-negative)
  stor_initial <- old_xs * max(old_thetas_i - old_thetas_i_meno_1, 0)

  # Denominator of the mass-balance term over dt
  denom <- stor_initial + q_vec * dt

  if (debug) {
    cat(
      "DEBUG f_euler:\n",
      "  dt =", dt, "\n",
      "  old_xs =", old_xs, "\n",
      "  stor_initial =", stor_initial, "\n",
      "  q_vec =", q_vec, "\n",
      "  denom =", denom, "\n",
      "  K_i - K_{i-1} =", Kvals_i - Kvals_i_meno_1, "\n"
    )
  }

  if (denom <= 1e-14) {
    dtheta_dt <- 0
  } else {
    # Net term: (q_i - (K_i - K_{i-1})) * (theta_i - theta_{i-1}^{new})
    num_part <- (q_vec - (Kvals_i - Kvals_i_meno_1)) * (old_thetas_i - new_thetas_i_meno_1)
    dtheta_dt <- num_part / denom + dtheta_dt_i_meno_1
    if (debug) {
      cat(
        "  num_part =", num_part, "\n",
        "  dtheta_dt_i_meno_1 =", dtheta_dt_i_meno_1, "\n",
        "  dtheta_dt =", dtheta_dt, "\n"
      )
    }
  }
  dtheta_dt
}

#' Merge two consecutive fronts by conservative storage recombination
#'
#' @description
#' Forms a single front from two adjacent ones (lower and upper) by conserving
#' storage relative to the preceding moisture \eqn{\theta_{\mathrm{prev}}}.
#' The resulting \eqn{\theta} is the maximum of the two (guaranteeing a positive jump).
#'
#' @param front_lower,front_upper lists with fields \code{$theta}, \code{$x}.
#' @param theta_r numeric; residual water content (used if no previous front).
#' @param old_lower,old_upper numeric; storages of the two segments (m).
#' @param theta_prima numeric; \eqn{\theta} below the lower front (or \eqn{\theta_r}).
#' @param delta numeric; absolute top reference [m].
#' @param debug logical; if TRUE, prints diagnostics.
#'
#' @return a list \code{list(theta=..., x=...)} for the merged front.
#'
#' @examples
#' merge_fronts(list(theta = 0.25, x = 1.2), list(theta = 0.27, x = 0.8),
#'   theta_r = 0.05, old_lower = 0.05, old_upper = 0.01,
#'   theta_prima = 0.10, delta = 0
#' )
#' @keywords internal
#' @noRd
merge_fronts <- function(front_lower,
                         front_upper,
                         theta_r,
                         old_lower,
                         old_upper,
                         theta_prima,
                         delta,
                         debug = FALSE) {
  theta_new <- max(front_lower$theta, front_upper$theta)
  storage_total <- max(old_lower, 0) + max(old_upper, 0)

  denom <- max(theta_new - theta_prima, 0)
  if (denom < 1e-14) {
    x_new <- max(front_lower$x, front_upper$x, delta)
  } else {
    x_new <- delta + storage_total / denom
  }

  if (debug) {
    cat(
      "DEBUG (merge_fronts): theta_new =", theta_new,
      " theta_prima =", theta_prima,
      " storage_total =", storage_total,
      " x_new =", x_new, "\n"
    )
  }
  list(theta = theta_new, x = x_new)
}

#' Global merging routine over a sequence of fronts
#'
#' @description
#' Applies conservative merging to the full front list until no geometric
#' overtake or near-equality in \eqn{\theta} remains. Enforces monotone depths
#' and removes inversions by recomputing a single effective front.
#'
#' @details
#' Merge triggers:
#' (i) geometric overtake \eqn{x_i \ge x_{i-1}};
#' (ii) near-equality \eqn{|\theta_i-\theta_{i-1}| < \texttt{merge\_tol}};
#' (iii) inversion \eqn{\theta_i < \theta_{i-1}}.
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param theta_r numeric; residual water content.
#' @param merge_tol numeric; absolute tolerance for merging in \eqn{\theta}.
#' @param old_storages reserved (pass \code{NULL}); kept for API compatibility.
#' @param debug logical; print diagnostics if TRUE.
#' @param delta numeric; absolute top reference [m].
#'
#' @return updated list of fronts after repeated merges.
#'
#' @examples
#' fs <- list(list(theta = 0.20, x = 1.2), list(theta = 0.21, x = 0.9), list(theta = 0.205, x = 0.85))
#' merging(fs, theta_r = 0.05, merge_tol = 1e-3, delta = 0)
#' @keywords internal
#' @noRd
merging <- function(fronts, theta_r, merge_tol, old_storages = NULL, debug = FALSE, delta) {
  debug <- isTRUE(debug)

  pp <- function(fs, label) {
    if (debug) {
      if (exists("pretty_print_fronts", mode = "function")) {
        pretty_print_fronts(fs, label = label)
      } else {
        cat(label, "\n")
        df <- data.frame(
          i     = seq_along(fs),
          theta = vapply(fs, `[[`, numeric(1), "theta"),
          x     = vapply(fs, `[[`, numeric(1), "x")
        )
        print(df, row.names = FALSE)
      }
    }
  }

  if (debug) {
    cat("DEBUG (merging): start with", length(fronts), "fronts\n")
    pp(fronts, "=== INITIAL FRONTS ===")
  }

  repeat {
    merged_any <- FALSE
    i <- 2L
    while (i <= length(fronts)) {
      x_lo <- fronts[[i - 1]]$x
      x_up <- fronts[[i]]$x
      th_lo <- fronts[[i - 1]]$theta
      th_up <- fronts[[i]]$theta

      sorpasso <- (x_up >= x_lo)
      quasi_ug <- (abs(th_up - th_lo) < merge_tol)
      inversione <- (th_up < th_lo)

      if (debug) {
        cat(sprintf(
          "DEBUG (merging): compare i-1=%d and i=%d | x[i-1]=%g, x[i]=%g, th[i-1]=%g, th[i]=%g, merge=%s\n",
          i - 1, i, x_lo, x_up, th_lo, th_up,
          ifelse(sorpasso || quasi_ug || inversione, "YES", "no")
        ))
      }

      if (sorpasso || quasi_ug || inversione) {
        theta_prima <- if (i == 2L) theta_r else fronts[[i - 2]]$theta

        S_lower <- (x_lo - delta) * max(th_lo - theta_prima, 0)

        # Upper storage is signed if there is an inversion, else non-negative
        if (inversione) {
          S_upper <- (x_up - delta) * (th_up - th_lo) # can be < 0
        } else {
          S_upper <- (x_up - delta) * max(th_up - th_lo, 0)
        }

        S_tot <- S_lower + S_upper
        theta_new <- max(th_lo, th_up)
        denom <- max(theta_new - theta_prima, 1e-14)

        x_new <- delta + S_tot / denom
        x_new <- max(x_new, delta) # do not go above the surface
        x_new <- min(x_new, max(x_lo, x_up)) # mild numerical stabilization

        if (debug) {
          cat("DEBUG (merging): Merging fronts", i - 1, "and", i, "\n")
          cat(sprintf(
            "DEBUG (merging): theta_prima=%g, S_lower=%g, S_upper=%g, S_tot=%g -> theta_new=%g, x_new=%g\n",
            theta_prima, S_lower, S_upper, S_tot, theta_new, x_new
          ))
        }

        fronts[[i - 1]] <- list(theta = theta_new, x = x_new)
        fronts <- fronts[-i]
        merged_any <- TRUE
      } else {
        i <- i + 1L
      }
    }
    if (!merged_any) break
  }

  if (debug) {
    cat("DEBUG (merging): end with", length(fronts), "fronts\n")
    pp(fronts, "=== MERGED FRONTS ===")
  }

  # Optional: post-merge monotonicity warning
  if (length(fronts) >= 2L) {
    th_all <- vapply(fronts, `[[`, numeric(1), "theta")
    if (any(diff(th_all) < -1e-12) && debug) {
      warning("Post-merge: non-monotone theta; check forcings/parameters.")
    }
  }

  fronts
}

#' Split a front list at a given absolute depth
#'
#' @description
#' Divides a monotone list of fronts into those at or above depth \code{Z}
#' (returned in \code{$upper}) and those below \code{Z} (\code{$lower}). If the
#' split point falls within a segment, a synthetic front at \code{x=Z} is
#' inserted to preserve piecewise-constant structure.
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param Z numeric; absolute depth [m] at which to split.
#'
#' @return a list with elements \code{upper} and \code{lower}, each a fronts list.
#'
#' @examples
#' divide_fronts(list(list(theta = 0.2, x = 0.6), list(theta = 0.3, x = 1.0)), Z = 0.8)
#' @keywords internal
#' @noRd
divide_fronts <- function(fronts, Z) {
  if (length(fronts) == 1) {
    if (fronts[[1]]$x <= Z) {
      return(list(upper = fronts, lower = list()))
    } else {
      upper_front <- list(theta = fronts[[1]]$theta, x = Z)
      lower_front <- list(theta = fronts[[1]]$theta, x = fronts[[1]]$x)
      return(list(upper = list(upper_front), lower = list(lower_front)))
    }
  }

  lower_fronts <- list()
  upper_fronts <- list()
  for (front in fronts) {
    if (front$x > Z) {
      lower_fronts <- c(lower_fronts, list(front))
    } else {
      upper_fronts <- c(upper_fronts, list(front))
    }
  }

  if (length(lower_fronts) > 0 && length(upper_fronts) > 0) {
    last_lower <- lower_fronts[[length(lower_fronts)]]
    if (!isTRUE(all.equal(last_lower$x, Z))) {
      dividing_front <- list(theta = last_lower$theta, x = Z)
      upper_fronts <- c(list(dividing_front), upper_fronts)
    }
  } else if (length(lower_fronts) > 0 && length(upper_fronts) == 0) {
    last_lower <- lower_fronts[[length(lower_fronts)]]
    dividing_front <- list(theta = last_lower$theta, x = Z)
    upper_fronts <- list(dividing_front)
  }

  list(upper = upper_fronts, lower = lower_fronts)
}

#' Plot a stepwise θ(z) profile with fixed axes
#' @importFrom graphics abline lines
#'
#' @description
#' Draws a step-plot of a piecewise-constant water content profile \eqn{\theta(z)}
#' over depth, using fixed \eqn{[\theta_r,\theta_s]} on the x-axis and an inverted
#' y-axis (depth increases downward). Optional horizontal markers can be added
#' for \code{Z_1} and \code{Z_2}.
#'
#' @param theta_r,theta_s numeric scalars; residual and saturated water contents.
#' @param profile list of points \code{list(list(theta=..., x=...), ...)} ordered
#'   by absolute depth \code{x} [m].
#' @param soil_depth numeric; bottom of the plotted domain [m].
#' @param timestep optional character; label shown in the plot title.
#' @param Z_1,Z_2 optional numeric; reference depths [m] to be drawn as dashed lines.
#' @param colore character; line color (default "black").
#'
#' @return Invisibly returns \code{NULL}. Side effect: produces a base R plot.
#'
#' @examples
#' prof <- list(list(theta = 0.18, x = 0.4), list(theta = 0.25, x = 0.8))
#' plot_single_profile(0.05, 0.45, prof, soil_depth = 1.0, timestep = "Example")
#' @export
plot_single_profile <- function(theta_r, theta_s, profile, soil_depth,
                                timestep = NULL, Z_1 = NA, Z_2 = NA, colore = NULL) {
  if (is.null(colore)) colore <- "black"

  # Extract and order points by increasing depth (surface -> deep)
  theta_i <- numeric()
  x_i <- numeric()
  for (fronts in profile) {
    x_i <- c(x_i, fronts$x)
    theta_i <- c(theta_i, fronts$theta)
  }

  ord <- order(x_i)
  theta_i <- theta_i[ord]
  x_i <- x_i[ord]

  # Build step endpoints and add surface/bottom anchors
  theta_i <- c(theta_i[1], theta_i, theta_r, theta_r)
  x_i <- c(0, x_i, x_i[length(x_i)], max(x_i))

  main <- if (is.null(timestep)) " " else paste0(" Profile at t: ", timestep)

  plot(theta_i, x_i,
    type = "s", lwd = 2, col = colore,
    xlab = expression(theta), ylab = "z [m]",
    main = main,
    xlim = c(theta_r, theta_s), ylim = c(soil_depth, 0)
  )
  if (!is.na(Z_1)) abline(h = Z_1, col = "red", lty = 2, lwd = 2)
  if (!is.na(Z_2)) abline(h = Z_2, col = "red", lty = 2, lwd = 2)
  lines(theta_i, x_i, type = "s", lwd = 2, col = colore)
}

#' Depth-averaged water content \eqn{\bar{\theta}(0\!-\!z^*)} for a stepwise profile
#'
#' @description
#' Computes the arithmetic, thickness-weighted average of a piecewise-constant
#' \eqn{\theta(z)} profile from the absolute top \eqn{\delta} to a depth limit
#' \eqn{z^*} (clamped to \code{L} if provided). Segments are interpreted with
#' respect to absolute depths; uncovered portions default to \eqn{\theta_r}.
#'
#' @param profile list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param depth_limit numeric; target depth \eqn{z^*} [m] for averaging.
#' @param theta_r numeric; residual water content (used for uncovered spans).
#' @param delta numeric; absolute depth origin [m] (default 0).
#' @param L optional numeric; lower bound of the domain [m]; used to cap \eqn{z^*}.
#' @param tol numeric; small tolerance for geometric comparisons.
#'
#' @return numeric scalar \eqn{\bar{\theta}(0\!-\!z^*)}; \code{NA} if \eqn{z^* \le \delta}.
#'
#' @examples
#' prof <- list(list(theta = 0.18, x = 0.2), list(theta = 0.25, x = 0.6))
#' compute_weighted_avg_theta_by_depth(prof,
#'   depth_limit = 0.5,
#'   theta_r = 0.05, delta = 0, L = 1.0
#' )
#' @export
compute_weighted_avg_theta_by_depth <- function(profile,
                                                depth_limit,
                                                theta_r,
                                                delta = 0,
                                                L = NULL,
                                                tol = 1e-12) {
  # Normalise depth_limit to domain bottom if provided
  if (!is.null(L)) depth_limit <- min(depth_limit, L)
  if (depth_limit <= delta + tol) {
    return(NA_real_)
  }

  # Empty profile implies theta_r everywhere
  if (!is.list(profile) || length(profile) == 0) {
    return(theta_r)
  }

  # Extract and order by depth (surface -> deep)
  x <- vapply(profile, function(fr) fr$x, numeric(1))
  th <- vapply(profile, function(fr) fr$theta, numeric(1))
  ord <- order(x)
  x <- x[ord]
  th <- th[ord]

  # Clamp fronts to [delta, z_star]
  z_star <- depth_limit
  x <- pmin(pmax(x, delta), z_star)

  # Layer bounds (z_top[i], z_bot[i]] with z_top[1]=delta, z_bot[i]=x[i]
  z_top <- c(delta, x[-length(x)])
  z_bot <- x

  # Fully covered layers
  fully_in <- which(z_bot <= z_star + tol)

  weighted_sum <- 0
  covered_thickness <- 0

  if (length(fully_in)) {
    layer_thick <- pmax(z_bot[fully_in] - z_top[fully_in], 0)
    weighted_sum <- weighted_sum + sum(th[fully_in] * layer_thick)
    covered_thickness <- covered_thickness + sum(layer_thick)
  }

  # Partial layer containing z_star (if any)
  part_idx <- which((z_top < z_star - tol) & (z_bot > z_star + tol))
  if (length(part_idx)) {
    i <- part_idx[1]
    rem <- z_star - max(z_top[i], delta)
    if (rem > tol) {
      weighted_sum <- weighted_sum + th[i] * rem
      covered_thickness <- covered_thickness + rem
    }
  }

  # Remaining uncovered thickness (delta -> z_star) defaults to theta_r
  total_span <- z_star - delta
  if (covered_thickness < total_span - tol) {
    weighted_sum <- weighted_sum + theta_r * (total_span - covered_thickness)
    covered_thickness <- total_span
  }

  if (covered_thickness <= tol) {
    return(NA_real_)
  }
  weighted_sum / total_span
}
