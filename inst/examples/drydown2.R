# =============================================================================
# Nonlinear drydown demo (Struthers-only) with variable ET
# Units: metres (depths, storages) and days (time); fluxes in m/day
# -----------------------------------------------------------------------------
# Purpose
#   Demonstrate how nonlinearity in depth-averaged drydown ( \bar{θ}(0–z) )
#   emerges from (i) gravitational infiltration–redistribution–drainage (IRD)
#   with moving wetting fronts and (ii) a two-threshold ET scheme (θ★/θ_wp)
#   applied after each daily IRD cycle, under variable potential ET.
# =============================================================================

# --- Assumed available in the workspace --------------------------------------
#   struthers_redistr_under(), apply_ET_paper_style(), plot_single_profile(),
#   compute_weighted_avg_theta_by_depth(), pretty_print_fronts() [optional]

# --- Soil & model parameters (accentuate nonlinearity) -----------------------
theta_r   <- 0        # residual water content [–]
theta_s   <- 0.45        # saturation [–]
beta      <- 1/8         # Brooks–Corey shape exponent (1/beta = 5)
Ks        <- 0.1         # saturated hydraulic conductivity [m/day]
L         <- 2.0         # column depth [m]
delta     <- 0.0         # surface reference [m]
tol_merge <- 1e-6        # front merge tolerance [–]

# --- ET thresholds and active zones ------------------------------------------
theta_wp <- 0.02         # wilting threshold for transpiration [–]
theta_st <- 0.20         # stress onset θ★ [–]
Z_root   <- 0.60         # active rooting depth from the surface [m]
Z_evap   <- 0.05         # surface layer active for bare-soil evaporation [m]
Cv       <- 1.0          # all potential ET allocated to transpiration (max drydown signal)

# --- Forcings (rainfall & potential ET) --------------------------------------
set.seed(42)
ndays <- 40

# Rainfall: three consecutive bursts (30 + 20 + 30 mm), then drydown
rain_mm        <- rep(0, ndays)
rain_mm[1:3]   <- c(30, 20, 30)
f_day          <- rain_mm / 1000  # convert to m/day

# Potential ET (variable; mm/day → m/day). Amplitude and short periods
# create visible variability in realised ET and soil water state.
t      <- 1:ndays
ETp_mm <- pmax(0, 6 + 2.5*sin(2*pi*t/9) + 0.6*cos(2*pi*t/5) + rnorm(ndays, 0, 0.3))
ETp    <- ETp_mm / 1000

# --- Time stepping ------------------------------------------------------------
dt_sub   <- 1/24              # hourly sub-steps [day]
subSteps <- as.integer(1/dt_sub)
dt_day   <- 1.0               # daily ET application [day]

# --- Initial profile (Day 0): uniform moist column (one deepest front) -------
theta0  <- 0.3
fronts0 <- list(list(theta = theta0, x = L))  # front 1 = deepest, at x = L

# --- Containers ---------------------------------------------------------------
profiles_daily_preET  <- vector("list", ndays + 1)  # include Day 0 at index 1
profiles_daily_postET <- vector("list", ndays + 1)
drainage_daily <- numeric(ndays)                   # daily drainage (m)
E_t_daily      <- numeric(ndays)                   # daily transpiration removed (m)
E_bs_daily     <- numeric(ndays)                   # daily bare-soil evaporation removed (m)

# Save Day 0
fronts <- fronts0
profiles_daily_preET[[1]]  <- fronts
profiles_daily_postET[[1]] <- fronts

# --- Helper: overlay a step-profile (post-ET) on an existing panel -----------
.add_profile_overlay <- function(theta_r, profile, col = "firebrick", lwd = 2) {
  th <- vapply(profile, function(fr) fr$theta, numeric(1))
  xx <- vapply(profile, function(fr) fr$x,     numeric(1))
  ord <- order(xx); th <- th[ord]; xx <- xx[ord]
  # close the step shape explicitly
  th <- c(th[1], th, theta_r, theta_r)
  xx <- c(0, xx, xx[length(xx)], max(xx))
  lines(th, xx, type = "s", lwd = lwd, col = col)
}

# =============================================================================
# Main loop: daily IRD (hourly) + daily ET (θ★/θ_wp)
# =============================================================================
for (d in 1:ndays) {
  
  # --- IRD over the day (hourly sub-steps) -----------------------------------
  f_t <- f_day[d]
  day_drain <- 0
  for (k in 1:subSteps) {
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
    day_drain <- day_drain + step$drainage      # accumulate volume per unit area [m]
  }
  drainage_daily[d] <- day_drain
  profiles_daily_preET[[d + 1]] <- fronts       # state before ET
  
  # --- Apply ET (once per day) with θ★/θ_wp scheme ---------------------------
  et <- apply_ET_paper_style(
    fronts      = fronts,
    ep          = ETp[d],     # m/day potential ET
    Cv          = Cv,         # fraction to transpiration
    dt          = dt_day,     # one-day ET step
    theta_star  = theta_st,
    theta_wp    = theta_wp,
    theta_r     = theta_r,
    theta_s     = theta_s,
    Z_root      = Z_root,
    Z_evap      = Z_evap,
    delta       = delta,
    debug       = FALSE,
    tol_merge   = tol_merge
  )
  fronts        <- et$fronts
  E_t_daily[d]  <- et$E_t
  E_bs_daily[d] <- et$E_bs
  
  profiles_daily_postET[[d + 1]] <- fronts      # state after ET
}

# =============================================================================
# Plots — profiles & depth-averaged θ
# =============================================================================

# --- (1) Profiles over time: pre- vs post-ET ---------------------------------
days_to_plot <- c(0, 1, 2, 3, 5, 10, 15, 25, 35)  # include Day 0

op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
par(mfrow = c(3, 3), mar = c(4, 4, 2.8, 1))

for (dd in days_to_plot) {
  idx <- dd + 1
  prof_pre  <- profiles_daily_preET[[idx]]
  prof_post <- profiles_daily_postET[[idx]]
  
  # Base panel: pre-ET (steelblue)
  if (dd == 0) {
    # --- Initial condition (start) ---
    plot_single_profile(
      theta_r, theta_s, prof_pre,
      soil_depth = L,
      timestep   = "Start",
      Z_1 = NA, Z_2 = NA, colore = "steelblue"
    )
    legend("bottomright",
           legend = "Initial state (Day 0)",
           col = "steelblue", lwd = 2,
           bg = "white", bty = "o", cex = 0.9)
  } else {
    # --- Normal days ---
    plot_single_profile(
      theta_r, theta_s, prof_pre,
      soil_depth = L,
      timestep   = sprintf("Day %d (pre/post ET)", dd),
      Z_1 = NA, Z_2 = NA, colore = "steelblue"
    )
    
    # Overlay: post-ET (firebrick)
    .add_profile_overlay(theta_r, prof_post, col = "firebrick", lwd = 2)
    
    # Annotate top with rainfall and ET forcing
    mtext(sprintf("Rain = %g mm   |   ETp = %.1f mm",
                  rain_mm[dd], ETp_mm[dd]),
          side = 3, line = 0.2, cex = 0.85)
    
    legend("bottomright",
           legend = c("Pre-ET", "Post-ET"),
           col = c("steelblue", "firebrick"), lwd = 2,
           bg = "white", bty = "o", cex = 0.85)
  }
}

# --- (2) Depth-averaged θ time series (post-ET) ------------------------------
depths <- c(0.05, 0.10, 0.20, 0.40, 1.0)
depths <- depths[depths <= L]
nd_all <- ndays + 1  # Day 0..ndays

Theta_bar <- matrix(NA_real_, nrow = nd_all, ncol = length(depths))
colnames(Theta_bar) <- paste0("0–", depths, " m")

for (i in 1:nd_all) {
  prof_i <- profiles_daily_postET[[i]]
  for (j in seq_along(depths)) {
    Theta_bar[i, j] <- compute_weighted_avg_theta_by_depth(
      profile     = prof_i,
      depth_limit = depths[j],
      theta_r     = theta_r,
      delta       = delta,
      L           = L
    )
  }
}

par(mfrow = c(1, 1), mar = c(4, 4, 2.2, 1))
matplot(0:ndays, Theta_bar, type = "l", lwd = 2, lty = 1,
        xlab = "Day",
        ylab = expression(bar(theta)~"(0–z)"),
        main = expression(paste("Nonlinear drydown of ", bar(theta), " with variable ET")))
legend("topright", legend = colnames(Theta_bar),
       col = seq_len(ncol(Theta_bar)), lwd = 2, cex = 0.9, bg = "white")

# --- (3) Forcing and realised ET series --------------------------------------
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

plot(1:ndays, rain_mm, type = "h", lwd = 2,
     xlab = "Day", ylab = "Rain [mm]", main = "Rain")
abline(h = 0, lwd = 1)

plot(1:ndays, ETp_mm, type = "l", lwd = 2,
     xlab = "Day", ylab = "ETp [mm d^-1]", main = "Variable potential ET")

plot(1:ndays, (E_t_daily + E_bs_daily) * 1000, type = "h", lwd = 2,
     xlab = "Day", ylab = "ET removed [mm]", main = "Realised ET (daily)")
abline(h = 0, lwd = 1)
par(mfrow = c(1, 1))