#' @title Semi-implicit Picard solver for the simultaneous update of moisture fronts
#' @description
#' Solves, within a single sub-step, the coupled update of \eqn{\theta} at all
#' Struthers-like fronts using a fixed-point (Picard) iteration. At each
#' iteration it recomputes the hydraulic conductivity \eqn{K(\theta)}, partitions
#' the inter-layer flux vector \eqn{\mathbf{q}}, integrates the ODEs for
#' \eqn{\theta}, and clamps the result within \eqn{[\theta_r,\theta_s]}.
#'
#' Convergence is declared when \eqn{\max_i |\theta_i^{(k+1)}-\theta_i^{(k)}| <
#' \mathrm{tol\_theta}} or when the maximum number of iterations is reached.
#'
#' @details
#' **Algorithmic outline**
#'
#' 1. Given the current iterate \eqn{\theta^{(k)}}, compute \eqn{K(\theta^{(k)})}.
#' 2. Compute the flux partition \eqn{\mathbf{q}} among layers from the
#'    boundary forcing \eqn{f_t} and local conductivities (top-down allocation).
#' 3. Advance \eqn{\theta} over \code{dt_sub}:
#'    - For the top layer, explicit Euler on storage balance.
#'    - For inner layers, use \code{f_euler} which couples \eqn{i} and \eqn{i-1}
#'      consistently with the chosen discretisation.
#' 4. Clamp \eqn{\theta} to \eqn{[\theta_r,\theta_s]} and test convergence.
#'
#' If the caller detects non-convergence at the controller level, time
#' sub-stepping and retries can be applied externally (see
#' \code{struthers_redistr_under}).
#'
#' @param fronts list of fronts, each as \code{list(theta=..., x=...)} with
#'   \code{x} the depth-like coordinate and \code{theta} the volumetric content.
#' @param theta_r,theta_s residual and saturation water contents (-).
#' @param beta Brooks-Corey shape parameter used in \eqn{K(\theta)}.
#' @param Ks saturated hydraulic conductivity [m/day].
#' @param f_t top boundary supply (e.g., infiltration rate) [m/day].
#' @param dt_sub sub-step duration [day].
#' @param delta offset (thickness above the first Struthers layer) [m].
#' @param max_iter maximum Picard iterations.
#' @param tol_theta absolute convergence tolerance on \eqn{\theta}.
#' @param debug logical; if \code{TRUE} prints per-iteration diagnostics.
#'
#' @return A list with
#' \itemize{
#'   \item \code{new_theta}: numeric vector of updated \eqn{\theta} at fronts;
#'   \item \code{q_vec}: numeric vector of layer fluxes used during the solve;
#'   \item \code{converged}: logical flag;
#'   \item \code{n_iter}: number of iterations performed.
#' }
#'
#' @seealso \code{\link{struthers_redistr_under}}
#' @examples
#' ## Internal solver; see struthers_redistr_under() for end-to-end examples.
#' @export
.solve_q_step_picard <- function(fronts, theta_r, theta_s, beta, Ks,
                                 f_t, dt_sub, delta,
                                 max_iter = 20, tol_theta = 1e-8,
                                 debug = FALSE) {
  nF <- length(fronts)
  K_of <- function(th) Kfun(theta = th, theta_r = theta_r, theta_s = theta_s, beta = beta, Ks = Ks)

  old_theta <- vapply(fronts, `[[`, numeric(1), "theta")
  old_x <- vapply(fronts, `[[`, numeric(1), "x")

  theta_cur <- old_theta
  q_vec <- numeric(nF)
  dtheta_dt <- numeric(nF)

  for (it in seq_len(max_iter)) {
    # 1) K(theta) and flux partition
    Kvals <- vapply(seq_len(nF), function(i) K_of(theta_cur[i]), numeric(1))
    Kbelow <- c(0, Kvals[-nF])
    for (i in seq_len(nF)) q_vec[i] <- min(f_t, Kvals[i]) - min(f_t, Kbelow[i])

    # 2) Coupled ODEs (same structure as the production code)
    new_theta <- theta_cur # use theta_{i-1}^{new} when advancing i
    # i = 1
    thick1 <- max(old_x[1] - delta, 0)
    dtheta_dt[1] <- if (thick1 > 1e-12) (q_vec[1] - Kvals[1]) / thick1 else 0
    new_theta[1] <- old_theta[1] + dtheta_dt[1] * dt_sub
    # i >= 2
    if (nF >= 2) {
      for (i in 2:nF) {
        dtheta_dt[i] <- f_euler(
          dt = dt_sub, t = dt_sub,
          old_xs = max(old_x[i] - delta, 0),
          old_thetas_i = old_theta[i],
          old_thetas_i_meno_1 = old_theta[i - 1],
          q_vec = q_vec[i],
          Kvals_i = Kvals[i],
          Kvals_i_meno_1 = Kvals[i - 1],
          dtheta_dt_i_meno_1 = dtheta_dt[i - 1],
          new_thetas_i_meno_1 = new_theta[i - 1],
          debug = debug
        )
        new_theta[i] <- old_theta[i] + dtheta_dt[i] * dt_sub
      }
    }

    # 3) clamp + convergence test
    new_theta <- pmin(pmax(new_theta, theta_r), theta_s)
    maxdiff <- max(abs(new_theta - theta_cur))
    if (debug) cat(sprintf("DEBUG (Picard it=%d): max|Delta_theta|=%.3e\n", it, maxdiff))
    theta_cur <- new_theta
    if (maxdiff < tol_theta) {
      return(list(
        new_theta = theta_cur,
        q_vec = q_vec,
        converged = TRUE,
        n_iter = it
      ))
    }
  }

  # not converged
  list(new_theta = theta_cur, q_vec = q_vec, converged = FALSE, n_iter = max_iter)
}


#' @title Infiltration-redistribution-drainage step (Struthers-style fronts)
#' @description
#' Advances the Struthers front-based soil column by one sub-step under a given
#' top flux \eqn{f_t}, performing (i) possible creation of a new top front,
#' (ii) pre-merge of adjacent fronts, (iii) simultaneous update of \eqn{\theta}
#' at all fronts via a Picard solve, (iv) update of the front depths \eqn{x} by
#' exact storage balance, (v) post-merge, and (vi) deep drainage on the
#' shallowest (topmost) front with clamping \eqn{x\le L}.
#'
#' @details
#' **Step order (faithful to the i-r-d sequence in the original formulation)**
#'
#' 1. *Optional new top front:* if \eqn{f_t > K(\theta_\mathrm{top})}, create a
#'    thin head-front at \eqn{x\approx \delta} with
#'    \eqn{\theta_\mathrm{new} = \theta_r + (\theta_s-\theta_r)(f_t/K_s)^\beta}.
#' 2. *Pre-merge:* merge negligible layers (tolerance \code{tol_merge}).
#' 3. *Simultaneous \eqn{\theta} solve:* use Picard iteration
#'    (\code{.solve_q_step_picard}) unless \code{infill_all=TRUE}, in which case
#'    all the boundary flux is imposed at the topmost layer for one fixed update.
#' 4. *Update front depths:* from exact mass balance (Eq. 1 in code comments),
#'    using previous and updated \eqn{\theta}.
#' 5. *Apply state and post-merge:* write back \eqn{\theta,x} and merge again.
#' 6. *Deep drainage:* only on the shallowest stored layer; excess thickness
#'    \eqn{x_1-L} is converted into drainage volume proportional to the local
#'    contrast in \eqn{\theta} and \eqn{x_1} is clamped to \eqn{L}.
#'
#' Field capacity corrections are intentionally **not** applied here to keep the
#' sequencing faithful to i-r-d; if desired, they should be applied *after*
#' redistribution/merge and *before* ET at the caller level.
#'
#' @param fronts list of fronts, each \code{list(theta=..., x=...)} (at least one).
#' @param theta_r,theta_s residual and saturation water contents (-).
#' @param beta Brooks-Corey shape parameter used in \eqn{K(\theta)}.
#' @param Ks saturated hydraulic conductivity [m/day].
#' @param L total soil depth [m].
#' @param delta geometric offset above the first stored layer [m].
#' @param f_t imposed top flux [m/day] for this sub-step.
#' @param dt_sub sub-step duration [day].
#' @param infill_all logical; if \code{TRUE} all \eqn{f_t} is injected at the
#'   shallowest layer in one fixed update (debug/compatibility mode).
#' @param debug logical; if \code{TRUE} prints detailed diagnostics.
#' @param theta_field_capacity,Z_1 kept for backward compatibility (unused here).
#' @param tol_merge merging tolerance for fronts.
#' @param solver_max_iter,solver_tol_theta,solver_max_halvings Picard controls
#'   and time-halving retry budget used when \code{infill_all=FALSE}.
#'
#' @return A list with
#' \itemize{
#'   \item \code{fronts}: updated fronts after post-merge;
#'   \item \code{fronts_infil}: snapshot of fronts immediately after the
#'         "geometric" infiltration and pre-merge (useful for plotting);
#'   \item \code{drainage}: deep drainage volume [m] released in this sub-step.
#' }
#'
#' @references
#' Struthers-style front models for infiltration-redistribution-drainage;
#' Brooks, R.H. & Corey, A.T. (1966) *Hydraulic properties of porous media*.
#'
#' @seealso \code{\link{.solve_q_step_picard}}
#' @examples
#' ## Minimal call (using predefined Kfun, merging, f_euler utilities):
#' # out <- struthers_redistr_under(fronts, theta_r=0.05, theta_s=0.45, beta=1/8,
#' #                                Ks=1.0, L=2.0, delta=0.2, f_t=0.01, dt_sub=1/24)
#' @export
struthers_redistr_under <- function(
    fronts,
    theta_r, theta_s, beta, Ks,
    L = 2, delta = 0,
    f_t,
    dt_sub = 1,
    infill_all = FALSE, # kept for compatibility; if TRUE, inject all f_t at the top
    debug = FALSE,
    theta_field_capacity = NA_real_,
    Z_1 = NA_real_,
    tol_merge = 1e-6,
    # [solver] convergence controls
    solver_max_iter = 20,
    solver_tol_theta = 1e-8,
    solver_max_halvings = 2) {
  stopifnot(length(fronts) >= 1L)
  nF <- length(fronts)
  K_of <- function(th) Kfun(theta = th, theta_r = theta_r, theta_s = theta_s, beta = beta, Ks = Ks)

  # ---------- (1) Optional new top front ----------
  top_idx <- length(fronts)
  K_top <- K_of(fronts[[top_idx]]$theta)
  if (debug) cat("DEBUG K_top =", K_top, "  f_t =", f_t, "\n")

  if (f_t > K_top && Ks > 0 && f_t > 0) {
    theta_new <- theta_r + (theta_s - theta_r) * (f_t / Ks)^beta
    theta_new <- min(theta_new, theta_s)
    tol_theta <- 1e-8
    if (theta_new > fronts[[top_idx]]$theta + tol_theta) {
      fronts[[top_idx + 1]] <- list(theta = theta_new, x = delta + max(tol_merge, 1e-8))
      if (debug) cat("DEBUG: created NEW TOP FRONT theta =", theta_new, " at x approx delta \n")
    } else if (debug) {
      cat("DEBUG: skip NEW TOP FRONT (theta_new approx theta_top)\n")
    }
  }

  # ---------- (2) Pre-merge ----------
  fronts <- merging(
    fronts = fronts, theta_r = theta_r, delta = delta,
    merge_tol = max(tol_merge, 1e-8), old_storages = NULL, debug = debug
  )
  nF <- length(fronts)
  fronts_infil <- fronts # snapshot after "geometric" infiltration and pre-merge

  # ---------- (3) Simultaneous theta solve (Picard) ----------
  if (infill_all) {
    # force all q at the shallowest layer in one fixed update
    old_theta <- vapply(fronts, `[[`, numeric(1), "theta")
    old_x <- vapply(fronts, `[[`, numeric(1), "x")
    Kvals <- vapply(fronts, function(fr) K_of(fr$theta), numeric(1))
    q_vec <- rep(0, nF)
    q_vec[nF] <- f_t
    dtheta_dt <- numeric(nF)
    # i = 1
    thick1 <- max(old_x[1] - delta, 0)
    dtheta_dt[1] <- if (thick1 > 1e-12) (q_vec[1] - Kvals[1]) / thick1 else 0
    new_theta <- old_theta
    new_theta[1] <- old_theta[1] + dtheta_dt[1] * dt_sub
    converged <- TRUE
    n_iter <- 1
    if (nF >= 2) {
      for (i in 2:nF) {
        dtheta_dt[i] <- f_euler(
          dt = dt_sub, t = dt_sub,
          old_xs = max(old_x[i] - delta, 0),
          old_thetas_i = old_theta[i],
          old_thetas_i_meno_1 = old_theta[i - 1],
          q_vec = q_vec[i],
          Kvals_i = Kvals[i],
          Kvals_i_meno_1 = Kvals[i - 1],
          dtheta_dt_i_meno_1 = dtheta_dt[i - 1],
          new_thetas_i_meno_1 = new_theta[i - 1],
          debug = debug
        )
        new_theta[i] <- old_theta[i] + dtheta_dt[i] * dt_sub
      }
    }
    new_theta <- pmin(pmax(new_theta, theta_r), theta_s)
  } else {
    # fallback with sub-stepping if Picard fails to converge
    attempts <- 0
    dt_eff <- dt_sub
    fronts_eff <- fronts
    repeat {
      sol <- .solve_q_step_picard(fronts_eff, theta_r, theta_s, beta, Ks,
        f_t, dt_eff, delta,
        max_iter = solver_max_iter, tol_theta = solver_tol_theta,
        debug = debug
      )
      if (sol$converged || attempts >= solver_max_halvings) {
        new_theta <- sol$new_theta
        q_vec <- sol$q_vec
        converged <- sol$converged
        n_iter <- sol$n_iter
        break
      } else {
        if (debug) cat("DEBUG: Picard did not converge, halving dt and performing two substeps.\n")
        # perform two consecutive half-steps to consume dt_eff
        dt_half <- dt_eff / 2
        for (rep2 in 1:2) {
          sol2 <- .solve_q_step_picard(fronts_eff, theta_r, theta_s, beta, Ks,
            f_t, dt_half, delta,
            max_iter = solver_max_iter, tol_theta = solver_tol_theta,
            debug = debug
          )
          th2 <- sol2$new_theta
          for (i in seq_along(fronts_eff)) fronts_eff[[i]]$theta <- th2[i]
        }
        attempts <- attempts + 1
        new_theta <- vapply(fronts_eff, `[[`, numeric(1), "theta")
        q_vec <- sol2$q_vec
        converged <- sol2$converged
        n_iter <- sol2$n_iter
        break
      }
    }
  }

  # ---------- (4) Update x by exact storage balance (Eq. 1) ----------
  old_theta <- vapply(fronts, `[[`, numeric(1), "theta")
  old_x <- vapply(fronts, `[[`, numeric(1), "x")
  new_x <- old_x

  for (i in seq_len(nF)) {
    if (i == 1) {
      S_old <- max(old_x[1] - delta, 0) * max(old_theta[1] - theta_r, 0)
      S_new <- S_old + q_vec[1] * dt_sub
      denom <- max(new_theta[1] - theta_r, 0)
      if (S_new <= 1e-12 && denom <= 1e-12) {
        new_x[1] <- L
        new_theta[1] <- theta_r
      } else {
        new_x[1] <- delta + S_new / max(denom, 1e-14)
      }
    } else {
      S_old <- max(old_x[i] - delta, 0) * max(old_theta[i] - old_theta[i - 1], 0)
      S_new <- S_old + q_vec[i] * dt_sub
      denom <- max(new_theta[i] - new_theta[i - 1], 0)
      if (S_new <= 1e-12 && denom <= 1e-12) {
        # empty layer: keep numerically stable thickness
        new_x[i] <- new_x[i - 1]
      } else {
        new_x[i] <- delta + S_new / max(denom, 1e-14)
      }
    }
  }

  # ---------- (5) Apply state and post-merge ----------
  for (i in seq_len(nF)) {
    fronts[[i]]$theta <- new_theta[i]
    fronts[[i]]$x <- new_x[i]
  }

  fronts <- merging(
    fronts = fronts, theta_r = theta_r, delta = delta,
    merge_tol = max(tol_merge, 1e-8), old_storages = NULL, debug = debug
  )
  nF <- length(fronts)

  # ---------- (6) Drainage only on the shallowest stored layer ----------
  drainage <- 0
  if (nF >= 1) {
    x_last <- fronts[[1]]$x
    if (x_last > L) {
      th_last <- fronts[[1]]$theta
      th_below <- if (nF == 1) theta_r else fronts[[nF - 1]]$theta
      thdiff <- max(th_last - th_below, 0)
      excess_h <- x_last - L
      drainage <- drainage + max(excess_h, 0) * thdiff
      fronts[[1]]$x <- L
    }
  }

  # Remove tiny internal layers at theta_r
  nF <- length(fronts)
  if (nF > 1) {
    keep <- rep(TRUE, nF)
    if (nF > 2) {
      for (i in 2:(nF - 1)) if (abs(fronts[[i]]$theta - theta_r) < 1e-12) keep[i] <- FALSE
    }
    fronts <- fronts[keep]
  }
  fronts <- clean_list(fronts)

  # Field Capacity note: intentionally omitted here to preserve the i-r-d-ET cycle.

  list(
    fronts = fronts,
    fronts_infil = fronts_infil,
    drainage = drainage
  )
}
