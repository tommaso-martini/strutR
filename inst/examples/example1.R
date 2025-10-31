## ======================================================================
## Example — Daily simulation of infiltration–redistribution–drainage
## ----------------------------------------------------------------------
## Units: metres and days (K_s and fluxes in m/day; depths in m)
## ----------------------------------------------------------------------
## This example demonstrates a 40-day synthetic simulation of the
## Struthers wetting-front model with hourly sub-stepping.
## It illustrates infiltration pulses, front progression, drainage
## generation, and depth-averaged moisture dynamics.
## ======================================================================

# --- 0) Check that all required routines are loaded -------------------
req_funs <- c("struthers_redistr_under", "Kfun", "f_euler", "merging",
              "merge_fronts", "clean_list", "plot_single_profile")
missing <- req_funs[!sapply(req_funs, exists)]
if (length(missing))
  stop("Missing functions in workspace: ", paste(missing, collapse = ", "))

set.seed(1)  # reproducibility

# --- 1) Soil and column parameters (Brooks–Corey) ---------------------
theta_r <- 0.00              # residual water content (-)
theta_s <- 0.45              # saturation (-)
beta    <- 1/5               # Brooks–Corey exponent (b in Struthers et al.)
Ks      <- 10.0              # saturated conductivity [m/day]
L       <- 4.0               # soil column depth [m]
delta   <- 0.0               # reference surface elevation [m]
tol_merge <- 1e-6            # tolerance for merging wetting fronts

# Optional: field capacity, for illustrative post-processing
theta_fc <- NA_real_          # e.g. 0.25 if you wish to fix FC on the lowest front

# --- 2) Initial condition ---------------------------------------------
# Uniformly dry column represented by one front at depth L
theta0  <- 0.12
fronts0 <- list(list(theta = theta0, x = L))  # front 1 = deepest, at base L

# --- 3) Daily rainfall forcing ----------------------------------------
ndays <- 40
rain_mm <- rep(0, ndays)
rain_mm[3:4]   <- 20   # 20 mm/day on days 3–4
rain_mm[10]    <- 30
rain_mm[15]    <- 40
rain_mm[30:32] <- 15
rain_mm[35]    <- 25
f_day <- rain_mm / 1000  # convert to m/day

# --- 4) Time stepping configuration -----------------------------------
dt_sub <- 1/24                     # hourly sub-steps (1 hour = 1/24 day)
sub_steps_per_day <- as.integer(1 / dt_sub)

profiles_daily <- vector("list", ndays)  # daily θ(z) profiles
drainage_daily <- numeric(ndays)         # cumulative drainage [m/day]
fronts <- fronts0

# --- 5) Core loop: integrate day-by-day -------------------------------
for (d in seq_len(ndays)) {
  f_t <- f_day[d]  # daily infiltration rate (constant within the day)
  for (k in seq_len(sub_steps_per_day)) {
    step <- struthers_redistr_under(
      fronts = fronts,
      theta_r = theta_r, theta_s = theta_s, beta = beta, Ks = Ks,
      L = L, delta = delta,
      f_t = f_t,                    # [m/day]
      dt_sub = dt_sub,              # [day]
      infill_all = FALSE,
      debug = FALSE,
      theta_field_capacity = theta_fc,
      Z_1 = NA_real_,
      tol_merge = tol_merge
    )
    fronts <- step$fronts
    drainage_daily[d] <- drainage_daily[d] + step$drainage  # m
  }
  profiles_daily[[d]] <- fronts  # store end-of-day profile
}

# --- 6) Daily profile visualization -----------------------------------
# Select representative days to display the infiltration–drydown sequence
days_to_plot <- c(1, 3, 4, 5, 10, 15, 16, 30, 33)

op <- par(no.readonly = TRUE)

par(mfrow = c(3,3), mar = c(4,4,3.2,1))  # un po' più spazio sopra

for (dd in days_to_plot) {
  prof <- profiles_daily[[dd]]
  plot_single_profile(theta_r, theta_s, prof, soil_depth = L,
                      timestep = paste("Day", dd),
                      Z_1 = NA, Z_2 = NA, colore = "black")
  legend("topright",
         legend = sprintf("Rain: %g mm", rain_mm[dd]),
         bty = "o", bg = "white", cex = 0.9, inset = 0.02, text.col = "gray20")
}
par(op)

# --- 7) Overlay of multiple profiles on one plot ----------------------
# Compare wetting front progression over a short event sequence
overlay_days <- c(6, 7, 8, 9, 10)
cols <- c("black", "blue", "darkgreen", "orange", "red")

# Base profile
d0 <- overlay_days[1]
plot_single_profile(theta_r, 0.16, # for a closer look
                    profiles_daily[[d0]], soil_depth = L,
                    timestep = paste("Day", d0),
                    colore = cols[1])

# Helper to overlay additional step plots
add_profile_lines <- function(theta_r, profile, colore = "red") {
  th <- vapply(profile, function(fr) fr$theta, numeric(1))
  x  <- vapply(profile, function(fr) fr$x, numeric(1))
  ord <- order(x)
  th <- th[ord]; x <- x[ord]
  th <- c(th[1], th, theta_r, theta_r)
  x  <- c(0, x, x[length(x)], max(x))
  lines(th, x, type = "s", lwd = 2, col = colore)
}

for (j in 2:length(overlay_days)) {
  dj <- overlay_days[j]
  add_profile_lines(theta_r, profiles_daily[[dj]], colore = cols[j])
}

legend("bottomright",
       legend = paste0("Day ", overlay_days, " (", rain_mm[overlay_days], " mm)"),
       col = cols, lwd = 2, cex = 0.8, bg = "white")

# Restore graphical parameters
par(op)

# --- 8) Daily drainage series -----------------------------------------
plot(seq_len(ndays), drainage_daily * 1000, type = "h", lwd = 2,
     xlab = "Day", ylab = "Drainage [mm/day]",
     main = "Daily drainage series")
abline(h = 0, lwd = 1)

# --- 9) Depth-averaged soil moisture θ̄(0–z) --------------------------
# Compute daily depth-weighted mean water contents for selected depths

depths <- c(0.05, 0.30, 0.60, 1, 2)
depths <- depths[depths <= L]

Theta_bar <- matrix(NA_real_, nrow = ndays, ncol = length(depths))
colnames(Theta_bar) <- paste0("0–", depths, " m")

for (d in seq_len(ndays)) {
  prof <- profiles_daily[[d]]
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

# --- Plot θ̄(0–z) time series ----------------------------------------
cols <- c("black", "blue", "darkgreen", "orange", "red")

op <- par(no.readonly = TRUE)
par(mar = c(4, 4, 2, 1))

plot(1:ndays, Theta_bar[, 1], type = "l", lwd = 2, col = cols[1],
     xlab = "Day", ylab = expression(bar(theta)~"(0–z)"),
     ylim = range(Theta_bar, na.rm = TRUE),
     main = expression(paste("Depth-averaged ", bar(theta),
                             " up to selected depths")))

if (ncol(Theta_bar) > 1) {
  for (j in 2:ncol(Theta_bar)) {
    lines(1:ndays, Theta_bar[, j], lwd = 2, col = cols[j])
  }
}

legend("topleft", inset = 0.02, bg = "white",
       legend = colnames(Theta_bar),
       col = cols[seq_len(ncol(Theta_bar))], lwd = 2, cex = 0.9)

par(op)

# --- 10) Optional: plot rainfall forcing ------------------------------
plot(seq_len(ndays), rain_mm, type = "h", lwd = 2,
     xlab = "Day", ylab = "Rainfall [mm/day]",
     main = "Rainfall forcing")
abline(h = 0, lwd = 1)