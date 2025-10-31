## ======================================================================
## Example — Daily simulation including evapotranspiration (θ★, θwp)
## ----------------------------------------------------------------------
## Units: metres and days
## Fluxes and Ks in m/day, depths in m, rainfall and ET in mm/day (converted internally).
## ----------------------------------------------------------------------
## This example demonstrates a 40-day synthetic simulation of the
## Struthers wetting-front model with daily rainfall, hourly
## infiltration–redistribution–drainage (IRD) integration, and daily
## evapotranspiration (ET) following the θ★–θwp scheme.
##
## It illustrates:
##   - dynamic infiltration and redistribution of wetting fronts,
##   - daily drainage generation,
##   - transpiration and bare-soil evaporation removal,
##   - evolution of depth-averaged soil moisture.
## ======================================================================

set.seed(1)

# --- 1) Soil hydraulic parameters (Brooks–Corey) ----------------------
theta_r <- 0.00         # residual water content [-]
theta_s <- 0.45         # saturated water content [-]
beta    <- 1/5          # Brooks–Corey shape exponent (b⁻¹)
Ks      <- 10.0         # saturated hydraulic conductivity [m/day]
L       <- 2.0          # total soil depth [m]
delta   <- 0.0          # surface reference level [m]
tol_merge <- 1e-6       # tolerance for merging fronts

# (Optional) Field capacity used only for post-processing illustration
theta_fc <- NA_real_

# --- 2) Initial state --------------------------------------------------
# Start from a moderately wet profile represented by a single front at L.
theta0  <- 0.35
fronts0 <- list(list(theta = theta0, x = L))   # front 1 = deepest, at base L

# --- 3) Daily forcing: rainfall, potential ET, and vegetation cover ----
ndays <- 40

# Rainfall [mm/day], converted to [m/day]
rain_mm <- rep(0, ndays)
rain_mm[3:4]   <- 30
rain_mm[10]    <- 15
rain_mm[15]    <- 40
rain_mm[28:30] <- 8
rain_mm[35]    <- 12
f_day <- rain_mm / 1000  # [m/day]

# Potential ET (ep_day) [m/day] and vegetation cover (Cv_day, 0–1)
ep_day <- rep(0.004, ndays)  # ~4 mm/day constant potential demand
Cv_day <- pmin(1, pmax(0, seq(0.3, 0.8, length.out = ndays)))  # linearly increasing canopy cover

# Transpiration stress thresholds
theta_star <- 0.16  # stress onset θ★ [-]
theta_wp   <- 0.01  # wilting θwp [-]

# Active depths for transpiration and evaporation
Z_root <- 1.0    # root-zone depth [m]
Z_evap <- 0.10   # evaporative surface layer [m]

# --- 4) Sub-daily integration and ET coupling -------------------------
dt_sub <- 1/24                  # sub-step = 1 hour = 1/24 day
sub_steps_per_day <- as.integer(1 / dt_sub)

profiles_daily_preET  <- vector("list", ndays)   # profile before ET
profiles_daily_postET <- vector("list", ndays)   # profile after ET
drainage_daily <- numeric(ndays)                 # daily drainage [m]
E_t_daily      <- numeric(ndays)                 # daily transpiration [mm]
E_bs_daily     <- numeric(ndays)                 # daily bare-soil evaporation [mm]

fronts <- fronts0

# --- 5) Main loop: daily IRD + ET ------------------------------------
for (d in seq_len(ndays)) {
  f_t <- f_day[d]  # infiltration flux for the current day
  
  # --- (a) Hourly infiltration–redistribution–drainage integration -----
  for (k in seq_len(sub_steps_per_day)) {
    step <- struthers_redistr_under(
      fronts = fronts,
      theta_r = theta_r, theta_s = theta_s, beta = beta, Ks = Ks,
      L = L, delta = delta,
      f_t = f_t,                      # [m/day]
      dt_sub = dt_sub,                # [day]
      infill_all = FALSE,
      debug = FALSE,
      theta_field_capacity = theta_fc,
      Z_1 = NA_real_,
      tol_merge = tol_merge
    )
    fronts <- step$fronts
    drainage_daily[d] <- drainage_daily[d] + step$drainage
  }
  
  # Snapshot before ET application
  profiles_daily_preET[[d]] <- fronts
  
  # --- (b) Apply daily ET using θ★–θwp scheme --------------------------
  et <- apply_ET(
    fronts,
    ep = ep_day[d], Cv = Cv_day[d], dt = 1.0,   # daily step
    theta_star = theta_star, theta_wp = theta_wp,
    theta_r = theta_r, theta_s = theta_s,
    Z_root = Z_root, Z_evap = Z_evap,
    delta = delta, debug = FALSE, tol_merge = tol_merge
  )
  
  # Update fronts and daily fluxes
  fronts <- et$fronts
  E_t_daily[d]  <- et$E_t  * 1000  # mm/day
  E_bs_daily[d] <- et$E_bs * 1000  # mm/day
  
  profiles_daily_postET[[d]] <- fronts
}

# --- 6) Daily profiles before and after ET -----------------------------
# Selected days for visual comparison
days_to_plot <- c(2, 3, 4, 5, 10, 15, 16, 30, 33)

op <- par(no.readonly = TRUE)
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

for (dd in days_to_plot) {
  # (a) Pre-ET profile in black
  plot_single_profile(
    theta_r = theta_r, theta_s = theta_s,
    profile = profiles_daily_preET[[dd]],
    soil_depth = L, timestep = paste0("Day ", dd),
    Z_1 = NA, Z_2 = NA, colore = "black"
  )
  
  # (b) Post-ET profile in red (overlaid)
  theta_i <- vapply(profiles_daily_postET[[dd]], function(fr) fr$theta, numeric(1))
  x_i     <- vapply(profiles_daily_postET[[dd]], function(fr) fr$x, numeric(1))
  ord <- order(x_i)
  theta_i <- theta_i[ord]; x_i <- x_i[ord]
  theta_i <- c(theta_i[1], theta_i, theta_r, theta_r)
  x_i     <- c(0, x_i, x_i[length(x_i)], max(x_i))
  lines(theta_i, x_i, type = "s", lwd = 2, col = "red")
  
  # (c) Top-right legend with white background and multi-line text
  legend(
    "topright",
    legend = c(
      "pre-ET",
      "post-ET",
      sprintf("Rain: %g mm",  rain_mm[dd]),
      sprintf("E_t:  %.1f mm", E_t_daily[dd]),
      sprintf("E_bs: %.1f mm", E_bs_daily[dd])
    ),
    col   = c("black", "red", NA, NA, NA),
    lwd   = c(2, 2, NA, NA, NA),
    bty   = "o", bg = "white", cex = 0.85, inset = 0.02,
    text.col = c("black", "red", rep("gray20", 3))
  )
}
par(op)

# --- 7) Time series of drainage and ET components ---------------------
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

plot(seq_len(ndays), drainage_daily * 1000, type = "h", lwd = 2,
     xlab = "Day", ylab = "Drainage [mm/day]",
     main = "Daily drainage series")
abline(h = 0, lwd = 1)

plot(seq_len(ndays), E_t_daily, type = "h", lwd = 2, col = "darkgreen",
     xlab = "Day", ylab = expression(E[t]~"[mm/day]"),
     main = "Daily transpiration (θ★, θ[wp])")
abline(h = 0, lwd = 1)

plot(seq_len(ndays), E_bs_daily, type = "h", lwd = 2, col = "orange",
     xlab = "Day", ylab = expression(E[bs]~"[mm/day]"),
     main = "Daily bare-soil evaporation")
abline(h = 0, lwd = 1)
par(op)

# --- 8) Depth-averaged soil moisture θ̄(0–z) --------------------------
# Compute the post-ET mean water content over several depths
stopifnot(exists("profiles_daily_postET"), length(profiles_daily_postET) > 0)

depths <- c(0.05, 0.10, 0.20, 0.40, 1.0)
depths <- depths[depths <= L]

Theta_bar <- matrix(NA_real_, nrow = ndays, ncol = length(depths))
colnames(Theta_bar) <- paste0("0–", depths, " m")

for (d in seq_len(ndays)) {
  prof <- profiles_daily_postET[[d]]
  for (j in seq_along(depths)) {
    Theta_bar[d, j] <- compute_weighted_avg_theta_by_depth(
      profile     = prof,
      depth_limit = depths[j],
      theta_r     = theta_r,
      delta       = delta,
      L           = L
    )
  }
}

# Plot the multi-depth time series of θ̄(0–z)
cols <- c("black", "blue", "darkgreen", "orange", "red")


par(mfrow = c(1, 1))

plot(1:ndays, Theta_bar[, 1], type = "l", lwd = 2, col = cols[1],
     xlab = "Day", ylab = expression(bar(theta)~"(0–z)"),
     ylim = range(Theta_bar, na.rm = TRUE),
     main = expression(paste("Post-ET mean ", bar(theta),
                             " over selected depths")))
if (ncol(Theta_bar) > 1) {
  for (j in 2:ncol(Theta_bar)) {
    lines(1:ndays, Theta_bar[, j], lwd = 2, col = cols[j])
  }
}

legend("topright", inset = 0.02, bg = "white",
       legend = colnames(Theta_bar),
       col = cols[seq_len(ncol(Theta_bar))], lwd = 2, cex = 0.9)
par(op)

# --- 9) Diagnostic plots: rainfall, potential ET, and canopy cover ----
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))

plot(seq_len(ndays), rain_mm, type = "h", lwd = 2,
     xlab = "Day", ylab = "Rainfall [mm/day]",
     main = "Rainfall forcing")
abline(h = 0, lwd = 1)

plot(seq_len(ndays), ep_day * 1000, type = "l", lwd = 2, col = "gray30",
     xlab = "Day", ylab = "Potential ET [mm/day]",
     main = expression(paste("Potential evapotranspiration ", e[p])))

plot(seq_len(ndays), Cv_day, type = "l", lwd = 2, col = "gray30",
     xlab = "Day", ylab = expression(C[v]~"[-]"),
     main = "Vegetation cover fraction")
par(op)