## ======================================================================
## Experiment — Sensitivity to hydraulic shape (β) and merge tolerance
## ----------------------------------------------------------------------
## Purpose
##   Compare two contrasting settings of the Brooks–Corey shape parameter
##   (via β, where the conductivity scales as ((θ-θr)/(θs-θr))^(1/β))
##   together with different merge tolerances, under a single 30 mm pulse
##   on day 1 followed by drydown.
##
## What this shows
##   • How a “sharper” K(θ) curve (smaller β ⇒ larger 1/β) accelerates early
##     percolation and modifies subsequent redistribution.
##   • How an aggressive merge tolerance can smooth (and sometimes overdamp)
##     the front structure, whereas a weak tolerance preserves fine layering.
##   • Consequences for depth-averaged θ̄(0–z) drydown curvature vs. linearity.
##
## Requirements (assumed exported & available in the workspace)
##   struthers_redistr_under, Kfun, f_euler, merging, merge_fronts,
##   clean_list, plot_single_profile, compute_weighted_avg_theta_by_depth
## ======================================================================

set.seed(1)

## -----------------------------
## Geometry, units, and forcing
## -----------------------------
L       <- 2.0           # column depth [m]
delta   <- 0.0           # reference elevation (surface) [m]
theta_r <- 0.05          # residual water content [-]
theta_s <- 0.45          # saturated water content [-]
Ks      <- 1.0           # saturated conductivity [m/day] (kept same)
dt_sub  <- 1/24          # hourly sub-step [day]
sub_steps <- as.integer(1/dt_sub)

ndays   <- 30
rain_mm <- rep(0, ndays); rain_mm[1] <- 30
f_day   <- rain_mm / 1000  # convert rainfall to [m/day]

## -----------------------------
## Initial condition
## -----------------------------
# Single deep front at the base L; moderately dry profile to favour motion
fronts0 <- list(list(theta = 0.18, x = L))

## ----------------------------------------------------------------------
## Scenario runner
##  beta: Brooks–Corey exponent as used in K(θ) = Ks * (...)^(1/beta)
##  tol_merge: θ-contrast tolerance for merging adjacent fronts
##  days_to_plot: selection for profile snapshots
##  label: descriptor propagated to figure titles and legends
## ----------------------------------------------------------------------
run_scenario <- function(beta,
                         tol_merge,
                         days_to_plot = c(1, 2, 3, 20, 30),
                         label = "Scenario",
                         verbose = FALSE) {
  
  fronts <- fronts0
  profiles <- vector("list", ndays)
  theta_bar_0_40  <- numeric(ndays)
  theta_bar_0_100 <- numeric(ndays)
  
  for (d in seq_len(ndays)) {
    f_t <- f_day[d]
    day_drain <- 0
    
    # Hourly IRD integration over one day
    for (k in seq_len(sub_steps)) {
      step <- struthers_redistr_under(
        fronts = fronts,
        theta_r = theta_r, theta_s = theta_s, beta = beta, Ks = Ks,
        L = L, delta = delta,
        f_t = f_t, dt_sub = dt_sub,
        infill_all = FALSE, debug = FALSE,
        theta_field_capacity = NA_real_,
        Z_1 = NA_real_,
        tol_merge = tol_merge
      )
      fronts    <- step$fronts
      day_drain <- day_drain + step$drainage
    }
    
    if (isTRUE(verbose)) {
      cat("\n")
      if (exists("pretty_print_fronts")) pretty_print_fronts(fronts, label = paste(label, "— day", d))
      cat("\n")
    }
    
    profiles[[d]] <- fronts
    theta_bar_0_40[d]  <- compute_weighted_avg_theta_by_depth(
      profile = profiles[[d]], depth_limit = 0.40,
      theta_r = theta_r, delta = delta, L = L
    )
    theta_bar_0_100[d] <- compute_weighted_avg_theta_by_depth(
      profile = profiles[[d]], depth_limit = 1.00,
      theta_r = theta_r, delta = delta, L = L
    )
  }
  
  list(beta = beta,
       tol_merge = tol_merge,
       profiles = profiles,
       theta_bar_0_40 = theta_bar_0_40,
       theta_bar_0_100 = theta_bar_0_100,
       label = label,
       days_to_plot = days_to_plot)
}

## -----------------------------
## Two contrasting scenarios
## -----------------------------
# Note: with beta = 1/5, the conductivity exponent is 1/beta = 5 (sharper K(θ)).
#       with beta = 1/10, exponent = 10 (even sharper).
scA <- run_scenario(beta = 1/5,  tol_merge = 1e-4, label = "A: 1/β = 5  | tol_merge = 1e-4")
scB <- run_scenario(beta = 1/10, tol_merge = 1e-8, label = "B: 1/β = 10 | tol_merge = 1e-8")

## ----------------------------------------------------------------------
## Profiles on selected days
##  - Upper row: Scenario A (steelblue)
##  - Lower row: Scenario B (firebrick)
##  - Each panel: top-right legend with white background; rainfall in mm
##  - Long labels are wrapped with strwrap() so they fit the box
## ----------------------------------------------------------------------

wrap_lines <- function(txt, width = 20) {
  # returns a character vector of wrapped lines
  strwrap(txt, width = width)
}

op <- par(no.readonly = TRUE)
on.exit(par(op), add = TRUE)

par(mfrow = c(2, length(scA$days_to_plot)), mar = c(4, 4, 3, 1))

## --- Scenario A (top row) ---
for (dd in scA$days_to_plot) {
  plot_single_profile(
    theta_r, theta_s, scA$profiles[[dd]],
    soil_depth = L, timestep = paste0("A-", dd),
    Z_1 = NA, Z_2 = NA, colore = "steelblue"
  )
  
  # Wrap scenario label and build per-line colours
  lab_lines_A <- wrap_lines(scA$label, width = 22)
  legend_text_A <- c(lab_lines_A, sprintf("Rain: %g mm", rain_mm[dd]))
  legend_cols_A <- c(rep("steelblue", length(lab_lines_A)), "gray20")
  
  legend("bottomright",
         legend   = legend_text_A,
         text.col = legend_cols_A,
         bty = "o", bg = "white", inset = 0.02, cex = 0.8, xpd = NA)
}

## --- Scenario B (bottom row) ---
for (dd in scB$days_to_plot) {
  plot_single_profile(
    theta_r, theta_s, scB$profiles[[dd]],
    soil_depth = L, timestep = paste0("B-", dd),
    Z_1 = NA, Z_2 = NA, colore = "firebrick"
  )
  
  # Wrap scenario label and build per-line colours
  lab_lines_B <- wrap_lines(scB$label, width = 22)
  legend_text_B <- c(lab_lines_B, sprintf("Rain: %g mm", rain_mm[dd]))
  legend_cols_B <- c(rep("firebrick", length(lab_lines_B)), "gray20")
  
  legend("bottomright",
         legend   = legend_text_B,
         text.col = legend_cols_B,
         bty = "o", bg = "white", inset = 0.02, cex = 0.8, xpd = NA)
}

## ----------------------------------------------------------------------
## Depth-averaged series — contrasting drydown shapes
##  Left: θ̄(0–0.40 m) ; Right: θ̄(0–1.00 m)
##  Expectation:
##    - Larger 1/β ⇒ steeper K(θ) ⇒ faster early redistribution and
##      potentially more concave early-time decay.
##    - Weak merging retains finer structure ⇒ slightly noisier θ̄.
## ----------------------------------------------------------------------
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

ylim_all <- range(
  c(scA$theta_bar_0_40,  scB$theta_bar_0_40,
    scA$theta_bar_0_100, scB$theta_bar_0_100),
  na.rm = TRUE
)

plot(1:ndays, scA$theta_bar_0_40, type = "l", lwd = 2, col = "steelblue",
     ylim = ylim_all, xlab = "Day",
     ylab = expression(bar(theta)),
     main = expression(paste("Mean ", bar(theta), " (0–0.40 m)")))
lines(1:ndays, scB$theta_bar_0_40, lwd = 2, col = "firebrick")
legend("topright",
       legend = c(scA$label, scB$label),
       col = c("steelblue", "firebrick"), lwd = 2,
       bg = "white", cex = 0.9, inset = 0.02)

plot(1:ndays, scA$theta_bar_0_100, type = "l", lwd = 2, col = "steelblue",
     ylim = ylim_all, xlab = "Day",
     ylab = expression(bar(theta)),
     main = expression(paste("Mean ", bar(theta), " (0–1.00 m)")))
lines(1:ndays, scB$theta_bar_0_100, lwd = 2, col = "firebrick")
legend("topright",
       legend = c(scA$label, scB$label),
       col = c("steelblue", "firebrick"), lwd = 2,
       bg = "white", cex = 0.9, inset = 0.02)

# return to a single panel for any subsequent plots
par(mfrow = c(1, 1))