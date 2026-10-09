test_that("hmm.filter, hmm.smoother, and calc.track run end-to-end on synthetic data", {
  dates <- as.POSIXct(c("2020-01-01 00:00:00", "2020-01-02 00:00:00", "2020-01-03 00:00:00"), tz = "UTC")
  
  # Setup grid
  sp.lim <- list(lonmin = -75, lonmax = -70, latmin = 30, latmax = 35)
  locs.grid <- setup.locs.grid(sp.lim, res = "one")
  g <- locs.grid
  
  # Setup L array of dimension (T, ncol, nrow)
  # Notice: in hmm.filter, dim(L)[1] is T, and spatial dimensions match g
  col <- dim(g$lon)[2]
  row <- dim(g$lon)[1]
  T_steps <- length(dates)
  
  L <- array(1e-5, dim = c(T_steps, col, row))
  # Set a probability peak at each time step
  L[1, 2, 2] <- 1
  L[2, 3, 3] <- 1
  L[3, 4, 4] <- 1
  
  # Normalize slices
  for (t in 1:T_steps) {
    L[t,,] <- L[t,,] / sum(L[t,,])
  }
  
  # Setup 1-state and 2-state kernels
  K1 <- HMMoce:::gausskern.pg(siz = 5, sigma = 1, muadv = 0)
  K2 <- HMMoce:::gausskern.pg(siz = 5, sigma = 2, muadv = 0)
  P2 <- matrix(c(0.8, 0.2, 0.2, 0.8), nrow = 2, byrow = TRUE)
  
  # 1. 2-State Filter & Smoother
  f2 <- hmm.filter(g = g, L = L, K = list(K1, K2), P = P2, m = 2)
  expect_named(f2, c("phi", "pred", "psi"))
  expect_equal(dim(f2$phi), c(2, T_steps, col, row))
  expect_equal(length(f2$psi), T_steps - 1)
  
  s2 <- hmm.smoother(f = f2, L = L, K = list(K1, K2), P = P2)
  expect_equal(dim(s2), c(2, T_steps, col, row))
  
  # 2. Track Calculation
  iniloc <- data.frame(
    year = c(2020, 2020),
    month = c(1, 1),
    day = c(1, 3),
    lon = c(-74, -71),
    lat = c(31, 34)
  )
  tr <- calc.track(distr = s2, g = g, dateVec = dates, iniloc = iniloc, method = "mean")
  expect_s3_class(tr, "data.frame")
  expect_true(all(c("lon", "lat") %in% names(tr)))
  expect_equal(nrow(tr), T_steps)
})
