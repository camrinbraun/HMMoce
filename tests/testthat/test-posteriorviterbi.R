# Reference cell-by-cell Viterbi implementation (pre-optimization logic) for exact parity verification
ref_posteriorviterbi <- function(pars, s, g, iniloc, dateVec, lims.lat, lims.lon, percentile = 99.9, approach = "separate") {
  if (approach == "joint") {
    if (dim(s)[1] > 1) {
      p.resid <- apply(s, c(1, 2), sum)[2, ]
      p.migrant <- 1 - p.resid
    } else {
      p.resid <- rep(NA, length.out = length(dateVec))
    }

    if (dim(s)[1] == 2) {
      p.1 <- apply(s[1, , , ], 1, sum)
      p.2 <- apply(s[2, , , ], 1, sum)
      s_norm <- array(NA, dim = dim(s)[c(3, 4, 2, 1)])
      for (i in 1:dim(s_norm)[3]) {
        s_norm[, , i, 1] <- (s[1, i, , ] / max(s[1, i, , ], na.rm = TRUE)) * p.1[i]
        s_norm[, , i, 2] <- (s[2, i, , ] / max(s[2, i, , ], na.rm = TRUE)) * p.2[i]
      }
      s_norm <- rowSums(s_norm, dims = 3, na.rm = TRUE)
      s_norm <- aperm(s_norm, c(3, 1, 2))
      s_viterbi <- array(0, dim = c(1, dim(s)[2], dim(s)[3], dim(s)[4]))
      for (tt in 1:dim(s)[2]) s_viterbi[1, tt, , ] <- s_norm[tt, , ]
      s <- s_viterbi
    }

    s_idx <- s
    for (i in 1:(dim(s)[2])) s_idx[1, i, , ] <- imager::renorm(s[1, i, , ], 0, 1)
    s[s_idx < (1 - percentile / 100)] <- 0
  } else {
    s_idx <- s
    for (i in 1:(dim(s)[2])) {
      s_idx[1, i, , ] <- imager::renorm(s[1, i, , ], 0, 1)
      if (dim(s)[1] == 2) s_idx[2, i, , ] <- imager::renorm(s[2, i, , ], 0, 1)
    }
    s[s_idx < (1 - percentile / 100)] <- 0
  }

  crop.lat.index <- which(g$lat[, 1] < lims.lat[1] | g$lat[, 1] > lims.lat[2])
  crop.lon.index <- which(g$lon[1, ] < lims.lon[1] | g$lon[1, ] > lims.lon[2])
  if (length(crop.lat.index) == 0 && length(crop.lon.index) == 0) {
    L.cropped <- s[, , , , drop = FALSE]
  } else if (length(crop.lat.index) == 0) {
    L.cropped <- s[, , -crop.lon.index, , drop = FALSE]
  } else if (length(crop.lon.index) == 0) {
    L.cropped <- s[, , , -crop.lat.index, drop = FALSE]
  } else {
    L.cropped <- s[, , -crop.lon.index, -crop.lat.index, drop = FALSE]
  }

  L.vit <- aperm(L.cropped, c(4, 3, 2, 1))
  L.vit <- aperm(L.vit, c(2, 1, 3, 4))
  L.vit <- array(
    HMMoce:::reverse(HMMoce:::reverse(HMMoce:::reverse(apply(L.vit, 2, c)))),
    dim = c(dim(L.vit)[2], dim(L.vit)[1], dim(L.vit)[3], dim(L.vit)[4])
  )

  longs.vit <- if (length(crop.lon.index) == 0) g$lon[1, ] else g$lon[1, -crop.lon.index]
  lats.vit <- if (length(crop.lat.index) == 0) g$lat[, 1] else g$lat[-crop.lat.index, 1]
  rows.vit <- dim(L.vit)[1]
  cols.vit <- dim(L.vit)[2]
  timesteps.vit <- dim(L.vit)[3]
  behavs.vit <- dim(L.vit)[4]
  ilo.vit <- which.min(abs(iniloc[1, "lon"] - longs.vit))
  ila.vit <- which.min(abs(iniloc[1, "lat"] - rev(lats.vit)))

  mult <- if (approach == "joint") 6 else 4
  sigmas <- if (length(pars) == 4) pars[1:2] else pars[1]
  sizes <- rep(ceiling(sigmas[1] * mult), 2)
  if (sizes[1] %% 2 == 0) sizes[1] <- sizes[1] + 1
  if (sizes[2] %% 2 == 0) sizes[2] <- sizes[2] + 1
  K1 <- HMMoce:::gausskern.pg(sizes[1], sigmas[1], muadv = 0)
  K <- if (length(pars) == 4) list(K1, HMMoce:::gausskern.pg(sizes[2], sigmas[2], muadv = 0)) else list(K1)

  if (approach == "separate") {
    if (dim(s)[1] == 2) {
      sv <- (-(apply(s[1, , , ], 1, sum) > apply(s[2, , , ], 1, sum)) + 1) + 1
      p.resid <- apply(s, c(1, 2), sum)[2, ]
    } else {
      sv <- NA
      p.resid <- rep(NA, length.out = length(dateVec))
    }
  }

  mpt.calc <- array(NA, dim = c(rows.vit, cols.vit, timesteps.vit, 2))
  mpt.calc[ila.vit, ilo.vit, 1, 1] <- which(row(mpt.calc[, , 1, 1]) == ila.vit & col(mpt.calc[, , 1, 1]) == ilo.vit)
  mpt.calc[ila.vit, ilo.vit, 1, 2] <- 1
  mpt.calc[, , 1, 2] <- log(mpt.calc[, , 1, 2])

  for (i in 2:timesteps.vit) {
    if (approach == "joint") {
      kern <- if (is.na(p.resid[1])) K[[1]] else K[[1]] * p.migrant[i] + K[[2]] * p.resid[i]
      L_day <- L.vit[, , i, 1]
    } else {
      kern <- if (is.na(sv[1]) || sv[i] == 1) K[[1]] else K[[2]]
      L_day <- if (behavs.vit == 1 || sv[i] == 1) L.vit[, , i, 1] else L.vit[, , i, 2]
    }

    kern.dim <- dim(kern)[1]
    kern.y.lower <- ceiling(kern.dim / 2)
    kern.y.upper <- rows.vit - floor(kern.dim / 2)
    kern.x.lower <- ceiling(kern.dim / 2)
    kern.x.upper <- cols.vit - floor(kern.dim / 2)
    kern.offset <- floor(kern.dim / 2)

    valid_prev <- which(!is.na(mpt.calc[, , i - 1, 1]), arr.ind = TRUE)
    if (nrow(valid_prev) > 0) {
      for (k_idx in seq_len(nrow(valid_prev))) {
        y <- valid_prev[k_idx, 1]
        x <- valid_prev[k_idx, 2]
        if (y > kern.y.lower && y < kern.y.upper) {
          ymin <- y - kern.offset; ymax <- y + kern.offset; kern.y <- 1:kern.dim
        } else if (y <= kern.y.lower) {
          ymin <- 1; y.offset <- kern.y.lower - y; kern.y <- (1:kern.dim)[-(1:y.offset)]; ymax <- y + kern.offset
        } else {
          ymax <- rows.vit; y.offset <- y - kern.y.upper; kern.y <- 1:(kern.dim - y.offset); ymin <- y - kern.offset
        }
        if (x > kern.x.lower && x < kern.x.upper) {
          xmin <- x - kern.offset; xmax <- x + kern.offset; kern.x <- 1:kern.dim
        } else if (x <= kern.x.lower) {
          xmin <- 1; x.offset <- kern.x.lower - x; kern.x <- (1:kern.dim)[-(1:x.offset)]; xmax <- x + kern.offset
        } else {
          xmax <- cols.vit; x.offset <- x - kern.x.upper; kern.x <- 1:(kern.dim - x.offset); xmin <- x - kern.offset
        }
        if (length(ymin:ymax) == kern.dim) kern.y <- 1:kern.dim
        if (length(xmin:xmax) == kern.dim) kern.x <- 1:kern.dim

        sub_upd <- L_day[ymin:ymax, xmin:xmax, drop = FALSE] * kern[kern.y, kern.x, drop = FALSE]
        pos_mask <- !is.na(sub_upd) & (sub_upd > 0)
        if (any(pos_mask)) {
          update.index <- y + (x - 1L) * rows.vit
          cand_prob <- log(sub_upd) + mpt.calc[y, x, i - 1, 2]
          cur_idx <- mpt.calc[ymin:ymax, xmin:xmax, i, 1]
          cur_prob <- mpt.calc[ymin:ymax, xmin:xmax, i, 2]
          better <- pos_mask & (is.na(cur_idx) | (cand_prob > cur_prob))
          if (any(better)) {
            cur_idx[better] <- update.index
            cur_prob[better] <- cand_prob[better]
            mpt.calc[ymin:ymax, xmin:xmax, i, 1] <- cur_idx
            mpt.calc[ymin:ymax, xmin:xmax, i, 2] <- cur_prob
          }
        }
      }
    }

    if (all(is.na(mpt.calc[, , i, 1]))) {
      best_prev_lin <- which.max(mpt.calc[, , i - 1, 2])
      if (length(best_prev_lin) > 0) {
        bp <- arrayInd(best_prev_lin[1], c(rows.vit, cols.vit))
        pos_l <- which(!is.na(L_day) & L_day > 0, arr.ind = TRUE)
        if (nrow(pos_l) == 0) pos_l <- bp
        d2 <- (pos_l[, 1] - bp[1, 1])^2 + (pos_l[, 2] - bp[1, 2])^2
        sig_eff <- max(sigmas[1], 1)
        lp <- log(pmax(L_day[pos_l], 1e-300)) - 0.5 * (d2 / (sig_eff^2)) + mpt.calc[bp[1, 1], bp[1, 2], i - 1, 2]
        mpt.calc[, , i, 1][pos_l] <- best_prev_lin[1]
        mpt.calc[, , i, 2][pos_l] <- lp
      }
    }
  }

  mpt <- matrix(NA, nrow = timesteps.vit, ncol = 2)
  end.loc <- arrayInd(which.max(mpt.calc[, , timesteps.vit, 2]), c(rows.vit, cols.vit))
  mpt[timesteps.vit, ] <- end.loc
  prev.loc <- mpt.calc[end.loc[, 1], end.loc[, 2], timesteps.vit, 1]
  for (i in (timesteps.vit - 1):1) {
    mpt[i, ] <- arrayInd(prev.loc, c(rows.vit, cols.vit))
    prev.loc <- mpt.calc[mpt[i, 1], mpt[i, 2], i, 1]
  }

  data.frame(
    "date" = dateVec,
    "lat" = lats.vit[(rows.vit + 1) - mpt[, 1]],
    "lon" = longs.vit[mpt[, 2]],
    "p.resid" = p.resid
  )
}

make_synthetic_grid <- function(ncol_g = 24, nrow_g = 22) {
  lon_v <- seq(-78, -66, length.out = ncol_g)
  lat_v <- seq(28, 39, length.out = nrow_g)
  list(
    lon = matrix(rep(lon_v, each = nrow_g), nrow = nrow_g),
    lat = matrix(rep(lat_v, times = ncol_g), nrow = nrow_g)
  )
}

test_that("hmm.posteriorviterbi matches reference identically for 2-state separate and joint", {
  set.seed(101)
  g <- make_synthetic_grid(24, 22)
  T_steps <- 12
  dates <- seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "1 day", length.out = T_steps)
  iniloc <- data.frame(day = c(1, 12), month = c(1, 1), year = c(2020, 2020), lon = c(-74, -69), lat = c(31, 36))

  s <- array(runif(2 * T_steps * ncol(g$lon) * nrow(g$lat))^4, dim = c(2, T_steps, ncol(g$lon), nrow(g$lat)))
  # Alternate state dominance across days
  w1 <- c(0.8, 0.7, 0.3, 0.2, 0.6, 0.9, 0.4, 0.1, 0.55, 0.75, 0.25, 0.65)
  for (t in seq_len(T_steps)) {
    s[1, t, , ] <- (s[1, t, , ] / sum(s[1, t, , ])) * w1[t]
    s[2, t, , ] <- (s[2, t, , ] / sum(s[2, t, , ])) * (1 - w1[t])
  }
  pars <- c(1.5, 0.5, 0.8, 0.7)

  # Cropped bounds
  lims.lat <- c(29.5, 37.5)
  lims.lon <- c(-76.5, -67.5)

  expect_equal(
    HMMoce:::hmm.posteriorviterbi(pars, s, g, iniloc, dates, lims.lat, lims.lon, percentile = 95, approach = "separate"),
    ref_posteriorviterbi(pars, s, g, iniloc, dates, lims.lat, lims.lon, percentile = 95, approach = "separate"),
    tolerance = 0
  )

  expect_equal(
    HMMoce:::hmm.posteriorviterbi(pars, s, g, iniloc, dates, lims.lat, lims.lon, percentile = 95, approach = "joint"),
    ref_posteriorviterbi(pars, s, g, iniloc, dates, lims.lat, lims.lon, percentile = 95, approach = "joint"),
    tolerance = 0
  )

  # Uncropped bounds
  expect_equal(
    HMMoce:::hmm.posteriorviterbi(pars, s, g, iniloc, dates, c(20, 50), c(-90, -50), approach = "separate"),
    ref_posteriorviterbi(pars, s, g, iniloc, dates, c(20, 50), c(-90, -50), approach = "separate"),
    tolerance = 0
  )
})

test_that("hmm.posteriorviterbi matches reference for 1-state models", {
  set.seed(202)
  g <- make_synthetic_grid(20, 20)
  T_steps <- 10
  dates <- seq(as.POSIXct("2020-02-01", tz = "UTC"), by = "1 day", length.out = T_steps)
  iniloc <- data.frame(day = c(1, 10), month = c(2, 2), year = c(2020, 2020), lon = c(-75, -70), lat = c(30, 35))

  s1 <- array(runif(T_steps * ncol(g$lon) * nrow(g$lat))^3, dim = c(1, T_steps, ncol(g$lon), nrow(g$lat)))
  for (t in seq_len(T_steps)) s1[1, t, , ] <- s1[1, t, , ] / sum(s1[1, t, , ])
  pars1 <- c(1.25)

  for (app in c("separate", "joint")) {
    res <- HMMoce:::hmm.posteriorviterbi(pars1, s1, g, iniloc, dates, c(29, 38), c(-77, -67), approach = app)
    ref <- ref_posteriorviterbi(pars1, s1, g, iniloc, dates, c(29, 38), c(-77, -67), approach = app)
    expect_equal(res, ref, tolerance = 0)
    expect_true(all(is.na(res$p.resid)))
  }
})

test_that("hmm.posteriorviterbi handles underflowed secondary kernels and disjoint-support fallback", {
  set.seed(303)
  g <- make_synthetic_grid(30, 30)
  T_steps <- 6
  dates <- seq(as.POSIXct("2020-03-01", tz = "UTC"), by = "1 day", length.out = T_steps)
  iniloc <- data.frame(day = c(1, 6), month = c(3, 3), year = c(2020, 2020), lon = c(-76, -68), lat = c(30, 37))

  # Large sigma1 vs tiny sigma2 causes outer cells of K2 to underflow to exact 0
  pars <- c(2.0, 0.08, 0.8, 0.8)
  s <- array(0, dim = c(2, T_steps, ncol(g$lon), nrow(g$lat)))
  # Day 1..3 localized near bottom-left; Day 4 jumps outside bounded kernel support to trigger fallback
  spots <- rbind(c(5, 5), c(6, 6), c(6, 7), c(25, 24), c(25, 25), c(26, 26))
  for (t in seq_len(T_steps)) {
    st_dom <- if (t %% 2 == 0) 2 else 1
    s[st_dom, t, spots[t, 1], spots[t, 2]] <- 0.85
    s[3 - st_dom, t, spots[t, 1], spots[t, 2]] <- 0.15
  }

  for (app in c("separate", "joint")) {
    expect_equal(
      HMMoce:::hmm.posteriorviterbi(pars, s, g, iniloc, dates, c(28, 39), c(-78, -66), approach = app),
      ref_posteriorviterbi(pars, s, g, iniloc, dates, c(28, 39), c(-78, -66), approach = app),
      tolerance = 0
    )
  }
})

test_that("hmm.posteriorviterbi recovers deterministic track from hmm.smoother output", {
  dates <- seq(as.POSIXct("2020-01-01 00:00:00", tz = "UTC"), by = "1 day", length.out = 5)
  sp.lim <- list(lonmin = -76, lonmax = -68, latmin = 29, latmax = 37)
  g <- setup.locs.grid(sp.lim, res = "one")
  col <- dim(g$lon)[2]
  row <- dim(g$lon)[1]

  L <- array(1e-8, dim = c(5, col, row))
  true_idx <- cbind(c(2, 3, 4, 5, 6), c(2, 3, 4, 5, 6))
  for (t in 1:5) {
    L[t, true_idx[t, 1], true_idx[t, 2]] <- 1
    L[t, , ] <- L[t, , ] / sum(L[t, , ])
  }

  K1 <- HMMoce:::gausskern.pg(siz = 5, sigma = 1.2, muadv = 0)
  K2 <- HMMoce:::gausskern.pg(siz = 5, sigma = 0.5, muadv = 0)
  P2 <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
  f2 <- hmm.filter(g = g, L = L, K = list(K1, K2), P = P2, m = 2)
  s2 <- hmm.smoother(f = f2, L = L, K = list(K1, K2), P = P2)

  iniloc <- data.frame(
    year = c(2020, 2020), month = c(1, 1), day = c(1, 5),
    lon = c(g$lon[1, 2], g$lon[1, 6]),
    lat = c(g$lat[2, 1], g$lat[6, 1])
  )

  tr <- HMMoce:::hmm.posteriorviterbi(
    pars = c(1.2, 0.5, 0.8, 0.8), s = s2, g = g, iniloc = iniloc,
    dateVec = dates, lims.lat = c(29, 37), lims.lon = c(-76, -68), approach = "separate"
  )

  expect_s3_class(tr, "data.frame")
  expect_named(tr, c("date", "lat", "lon", "p.resid"))
  expect_equal(tr$lon, g$lon[1, true_idx[, 1]])
  expect_equal(tr$lat, g$lat[true_idx[, 2], 1])
})
