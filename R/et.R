#' Build stepwise layers from a fronts profile (deep -> shallow)
#'
#' @description
#' Converts a Struthers-style list of fronts (deep to shallow) into a
#' stepwise-layer representation that is convenient for top-down operations
#' (e.g., root uptake and bare-soil evaporation limited to a shallow zone).
#'
#' @details
#' Given fronts \eqn{(x_i,\theta_i)} with \eqn{x} increasing with depth and
#' the convention "front 1 = deepest", the routine returns ascending bounds
#' \eqn{[\delta=b_0 < b_1 < \dots < b_{n_F}]} and the corresponding layer
#' water contents \eqn{\theta_j} in \emph{top-down} order.
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param theta_r numeric; residual water content \eqn{\theta_r}.
#' @param delta numeric; absolute top reference [m].
#'
#' @return A list with:
#' \item{bounds}{numeric vector of length \eqn{n_F+1} in ascending order.}
#' \item{theta_layer}{numeric vector of length \eqn{n_F} (top layer first).}
#'
#' @examples
#' .build_layers_from_fronts(list(list(theta = 0.2, x = 0.6), list(theta = 0.3, x = 1.0)), theta_r = 0.05)
#' @keywords internal
#' @noRd
.build_layers_from_fronts <- function(fronts, theta_r, delta = 0) {
  nF <- length(fronts)
  x <- vapply(fronts, `[[`, numeric(1), "x")
  th <- vapply(fronts, `[[`, numeric(1), "theta")

  # front 1 = deepest; return bounds from surface to deep
  bounds <- c(delta, rev(x)) # length nF+1
  theta_layer <- rev(th) # layer 1 is shallow/top

  list(bounds = bounds, theta_layer = theta_layer)
}

#' Effective thickness of layer j within an active top zone
#'
#' @description
#' Computes the intersection thickness of the \code{j}-th layer with the
#' active top zone \eqn{[\delta, \delta + Z_{\mathrm{lim}}]}.
#'
#' @param j integer; layer index in top-down order (1..nL).
#' @param bounds numeric; layer bounds as returned by \code{.build_layers_from_fronts()}.
#' @param Zlim numeric; active zone thickness [m].
#' @param delta numeric; absolute top reference [m].
#'
#' @return numeric thickness [m] of layer \code{j} within the active zone.
#'
#' @examples
#' # internal helper
#' @keywords internal
#' @noRd
.layer_thickness_in_zone <- function(j, bounds, Zlim, delta = 0) {
  if (Zlim <= 0) {
    return(0)
  }
  ztop <- delta
  zbot <- delta + Zlim
  top_j <- max(bounds[j - 1 + 1], ztop) # bounds is 1..nF+1
  bot_j <- min(bounds[j + 1], zbot)
  max(bot_j - top_j, 0)
}

#' Storage above a base water content within a shallow top zone
#'
#' @description
#' Computes the total volumetric storage above a specified base threshold
#' within the active top zone. Typical bases are \eqn{\theta_{\mathrm{wp}}}
#' for transpiration and \eqn{\theta_r} for bare-soil evaporation.
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param Zlim numeric; active zone thickness [m] from \code{delta}.
#' @param base numeric; base water content threshold.
#' @param theta_r numeric; residual water content for continuity constraints.
#' @param delta numeric; absolute top reference [m].
#'
#' @return numeric storage [m] available above \code{base} in the active zone.
#'
#' @examples
#' # internal helper
#' @keywords internal
#' @noRd
.storage_above_base_in_zone <- function(fronts, Zlim, base, theta_r, delta = 0) {
  if (length(fronts) == 0 || Zlim <= 0) {
    return(0)
  }
  layers <- .build_layers_from_fronts(fronts, theta_r, delta)
  bounds <- layers$bounds
  tl <- layers$theta_layer
  nL <- length(tl)

  # Below-layer thetas to enforce piecewise-constant continuity
  tl_below <- c(theta_r, tl[-length(tl)])

  vol <- 0
  for (j in seq_len(nL)) {
    t_eff <- .layer_thickness_in_zone(j, bounds, Zlim, delta)
    if (t_eff <= 0) next
    base_j <- max(base, tl_below[j])
    vol <- vol + max(tl[j] - base_j, 0) * t_eff
  }
  vol
}

#' Remove a target water volume from the top zone with lower bounds on theta
#'
#' @description
#' Applies a top-down removal of a prescribed volume within
#' \eqn{[\delta, \delta + Z_{\mathrm{lim}}]}, respecting a minimum admissible
#' water content per layer (e.g., \eqn{\theta_{\mathrm{wp}}} for transpiration
#' or \eqn{\theta_r} for evaporation).
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param vol numeric; target volume to remove [m].
#' @param Zlim numeric; active zone thickness [m].
#' @param min_base_theta numeric; per-layer lower bound for water content.
#' @param theta_r numeric; residual water content for continuity.
#' @param delta numeric; absolute top reference [m].
#' @param debug logical; print stepwise diagnostics if TRUE.
#'
#' @return A list with:
#' \item{fronts}{updated fronts list after removal.}
#' \item{removed}{actual removed volume [m].}
#'
#' @examples
#' # internal helper
#' @keywords internal
#' @noRd
.remove_water_from_top_zone <- function(fronts, vol, Zlim, min_base_theta, theta_r,
                                        delta = 0, debug = FALSE) {
  if (vol <= 0 || Zlim <= 0 || length(fronts) == 0) {
    return(list(fronts = fronts, removed = 0))
  }
  layers <- .build_layers_from_fronts(fronts, theta_r, delta)
  bounds <- layers$bounds
  tl <- layers$theta_layer
  nL <- length(tl)

  tl_below <- c(theta_r, tl[-length(tl)]) # continuity

  removed <- 0
  for (j in seq_len(nL)) { # top-down: layer 1 is shallowest
    if (removed >= vol - 1e-14) break
    t_eff <- .layer_thickness_in_zone(j, bounds, Zlim, delta)
    if (t_eff <= 0) next

    min_theta_j <- max(min_base_theta, tl_below[j])
    cap_j <- max(tl[j] - min_theta_j, 0) * t_eff
    if (cap_j <= 0) next

    take <- min(vol - removed, cap_j)
    dtheta <- take / max(t_eff, 1e-14)
    tl[j] <- max(tl[j] - dtheta, min_theta_j)
    removed <- removed + take

    if (debug) {
      cat(sprintf(
        "ET remove: layer=%d (top), t_eff=%.6g, cap=%.6g, take=%.6g, theta->%.6g\n",
        j, t_eff, cap_j, take, tl[j]
      ))
    }
  }

  # Map back to fronts ordering (deep -> shallow)
  th_new <- rev(tl)
  for (i in seq_along(fronts)) fronts[[i]]$theta <- th_new[i]

  list(fronts = fronts, removed = removed)
}

#' Mean water content in the root zone \eqn{\bar{\theta}(0\!-\!Z_{\mathrm{root}})}
#'
#' @description
#' Computes the thickness-weighted mean \eqn{\bar{\theta}} within the top root-zone
#' depth \eqn{Z_{\mathrm{root}}}, defaulting uncovered spans to \eqn{\theta_r}.
#' This is consistent with the stepwise representation used in the ET operators.
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param Z_root numeric; root-zone depth [m] from \code{delta}.
#' @param theta_r numeric; residual water content.
#' @param delta numeric; absolute top reference [m].
#'
#' @return numeric scalar \eqn{\bar{\theta}} over \eqn{[0,Z_{\mathrm{root}}]}.
#' Returns \code{NA} if \eqn{Z_root \le 0}.
#'
#' @examples
#' mean_theta_in_root_zone(list(list(theta = 0.20, x = 0.3), list(theta = 0.28, x = 0.8)),
#'   Z_root = 0.5, theta_r = 0.05
#' )
#' @export
mean_theta_in_root_zone <- function(fronts, Z_root, theta_r, delta = 0) {
  if (Z_root <= 0) {
    return(NA_real_)
  }
  if (length(fronts) == 0) {
    return(theta_r)
  }

  # Deep->shallow to top->down layers
  x <- vapply(fronts, `[[`, numeric(1), "x")
  th <- vapply(fronts, `[[`, numeric(1), "theta")
  bounds <- c(delta, rev(x)) # b0=delta, b1=x[nF], ..., b_nF=x[1]
  tl <- rev(th) # top layer first

  ztop <- delta
  zbot <- delta + Z_root

  numer <- 0
  covered <- 0
  nL <- length(tl)

  for (j in seq_len(nL)) {
    top_j <- max(bounds[j], ztop)
    bot_j <- min(bounds[j + 1], zbot)
    t_eff <- max(bot_j - top_j, 0)
    if (t_eff > 0) {
      numer <- numer + tl[j] * t_eff
      covered <- covered + t_eff
    }
    if (bounds[j + 1] >= zbot - 1e-14) break
  }

  # Uncovered remainder defaults to theta_r
  if (covered < Z_root - 1e-14) {
    numer <- numer + theta_r * (Z_root - covered)
  }

  numer / Z_root
}

#' Paper-style ET operator on a fronts profile (single time step)
#'
#' @description
#' Applies a transpiration-evaporation partition consistent with common
#' ecohydrological practice: transpiration acts on a prescribed root-zone
#' depth with a \eqn{\theta^\star/\theta_{\mathrm{wp}}} scheme (potential phase
#' down to \eqn{\theta^\star}, then an exponential stress regime toward
#' \eqn{\theta_{\mathrm{wp}}}); bare-soil evaporation acts on a shallow
#' \eqn{Z_{\mathrm{evap}}} and is limited below by \eqn{\theta_r}.
#'
#' @details
#' Let the step durations be \eqn{\Delta t}, the potential demand be
#' \eqn{e_p} [m/day], and the transpiring fraction \eqn{C_v \in [0,1]}.
#' The potential step demands are:
#' \deqn{D_t = e_p C_v \Delta t,\quad D_{bs} = e_p (1-C_v) \Delta t.}
#' For transpiration, the available potential removal to reach
#' \eqn{\bar{\theta}\to\theta^\star} is \eqn{S_\star = \max(\bar{\theta}-\theta^\star,0) Z_{\mathrm{root}}}.
#' The stressed contribution is \eqn{S_{\mathrm{cap}} (1 - e^{-D_{\mathrm{rem}}/S_{\mathrm{cap}}})}
#' with \eqn{S_{\mathrm{cap}}=(\theta^\star-\theta_{\mathrm{wp}})Z_{\mathrm{root}}}.
#' Both TR and EV are capped by the actually available storage above their
#' respective lower bounds, and water is removed top-down.
#'
#' @param fronts list of fronts \code{list(list(theta=..., x=...), ...)}.
#' @param ep numeric; potential evaporation/ET [m/day] for the step.
#' @param Cv numeric in \eqn{[0,1]}; transpiring fraction.
#' @param dt numeric; step duration [day].
#' @param theta_star,theta_wp numeric; transpiration thresholds.
#' @param theta_r,theta_s numeric; residual and saturated water contents.
#' @param Z_root numeric; active root-zone depth [m] from \code{delta}.
#' @param Z_evap numeric; active bare-soil evaporation depth [m] from \code{delta}.
#' @param delta numeric; absolute top reference [m].
#' @param debug logical; print detailed diagnostics if TRUE.
#' @param tol_merge numeric; merge tolerance passed to the fronts merging routine.
#'
#' @return A list with:
#' \item{fronts}{updated fronts list after ET.}
#' \item{E_t}{actual transpiration volume removed [m].}
#' \item{E_bs}{actual bare-soil evaporation volume removed [m].}
#'
#' @examples
#' prof <- list(list(theta = 0.25, x = 0.2), list(theta = 0.30, x = 0.6))
#' apply_ET_paper_style(prof,
#'   ep = 4 / 1000, Cv = 0.7, dt = 1,
#'   theta_star = 0.20, theta_wp = 0.08,
#'   theta_r = 0.05, theta_s = 0.45,
#'   Z_root = 0.4, Z_evap = 0.03
#' )
#' @export
apply_ET_paper_style <- function(fronts,
                                 ep, Cv, dt,
                                 theta_star, theta_wp,
                                 theta_r, theta_s,
                                 Z_root, Z_evap,
                                 delta = 0,
                                 debug = FALSE,
                                 tol_merge = 1e-6) {
  if (length(fronts) == 0) {
    return(list(fronts = fronts, E_t = 0, E_bs = 0))
  }

  # --- Potential demands for the step ---
  Dt <- max(ep, 0) * max(min(Cv, 1), 0) * max(dt, 0) # m
  Dbs <- max(ep, 0) * max(1 - Cv, 0) * max(dt, 0) # m

  # =========================
  # (A) TRANSPIRATION (theta_star, theta_wp)
  # =========================
  E_t_target <- 0
  E_t_eff <- 0

  if (Dt > 0 && Z_root > 0) {
    th_bar <- mean_theta_in_root_zone(fronts, Z_root, theta_r, delta)

    # Potential phase to theta_star
    S_star <- max((th_bar - theta_star), 0) * Z_root
    if (debug) cat(sprintf("TR: th_bar=%.5f, S_star=%.6g m\n", th_bar, S_star))

    if (S_star >= Dt) {
      E_t_target <- Dt
    } else {
      D_rem <- Dt - S_star
      S_stress_cap <- max(theta_star - theta_wp, 0) * Z_root
      E_stress <- if (S_stress_cap > 0) S_stress_cap * (1 - exp(-D_rem / S_stress_cap)) else 0
      E_t_target <- S_star + E_stress
    }

    # Cap by actually available storage above theta_wp in the root zone
    S_avail_wp <- .storage_above_base_in_zone(fronts, Z_root, theta_wp, theta_r, delta)
    if (debug) cat(sprintf("TR: E_t_target=%.6g, S_avail_wp=%.6g m\n", E_t_target, S_avail_wp))
    E_t_target <- min(E_t_target, S_avail_wp)

    # Apply top-down with theta >= theta_wp
    rem_tr <- .remove_water_from_top_zone(fronts,
      vol = E_t_target,
      Zlim = Z_root,
      min_base_theta = theta_wp,
      theta_r = theta_r,
      delta = delta,
      debug = debug
    )
    fronts <- rem_tr$fronts
    E_t_eff <- rem_tr$removed

    if (debug) cat(sprintf("TR: removed=%.6g m\n", E_t_eff))
  }

  # =========================
  # (B) BARE-SOIL EVAPORATION (thetar)
  # =========================
  E_bs_target <- 0
  E_bs_eff <- 0

  if (Dbs > 0 && Z_evap > 0) {
    S_evap <- .storage_above_base_in_zone(fronts, Z_evap, theta_r, theta_r, delta)
    if (S_evap > 0) {
      # Exponential removal: S * (1 - exp(-D / S))
      E_bs_target <- S_evap * (1 - exp(-Dbs / S_evap))
      rem_ev <- .remove_water_from_top_zone(fronts,
        vol = E_bs_target,
        Zlim = Z_evap,
        min_base_theta = theta_r,
        theta_r = theta_r,
        delta = delta,
        debug = debug
      )
      fronts <- rem_ev$fronts
      E_bs_eff <- rem_ev$removed
    }
    if (debug) {
      cat(sprintf(
        "EV: S_evap=%.6g, target=%.6g, removed=%.6g m\n",
        S_evap, E_bs_target, E_bs_eff
      ))
    }
  }

  # =========================
  # (C) Light merge & cleanup
  # =========================
  if (length(fronts) >= 2L) {
    fronts <- merging(
      fronts = fronts, theta_r = theta_r, delta = delta,
      merge_tol = max(tol_merge, 1e-8), old_storages = NULL, debug = debug
    )
  }

  # Optional: re-check monotonicity
  if (length(fronts) >= 2L) {
    th <- vapply(fronts, `[[`, numeric(1), "theta")
    if (any(diff(th) < -1e-12) && debug) {
      cat("WARNING: theta inversion detected post-ET; re-merging.\n")
      fronts <- merging(
        fronts = fronts, theta_r = theta_r, delta = delta,
        merge_tol = max(tol_merge, 1e-8), old_storages = NULL, debug = debug
      )
    }
  }

  list(
    fronts = fronts,
    E_t = E_t_eff,
    E_bs = E_bs_eff
  )
}
