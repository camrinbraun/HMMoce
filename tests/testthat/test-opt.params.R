test_that("opt.params estimates parameters using optim and neg.log.lik.fun", {
  dates <- as.POSIXct(c("2020-01-01 00:00:00", "2020-01-02 00:00:00"), tz = "UTC")
  
  sp.lim <- list(lonmin = -75, lonmax = -72, latmin = 30, latmax = 33)
  locs.grid <- setup.locs.grid(sp.lim, res = "one")
  g <- locs.grid
  
  col <- dim(g$lon)[2]
  row <- dim(g$lon)[1]
  T_steps <- length(dates)
  
  L <- array(1e-4, dim = c(T_steps, col, row))
  L[1, 2, 2] <- 1
  L[2, 3, 3] <- 1
  for (t in 1:T_steps) L[t,,] <- L[t,,] / sum(L[t,,])
  
  # Test neg.log.lik.fun directly
  pars <- c(2, 0.5, 0.8, 0.8)
  nll <- HMMoce:::neg.log.lik.fun(pars = pars, g = g, L = L)
  expect_true(is.numeric(nll))
  expect_false(is.nan(nll))
  
  # Test pos.log.lik.fun (for GA fitness)
  pos_ll <- HMMoce:::pos.log.lik.fun(pars = pars, g = g, L = L)
  expect_equal(pos_ll, -nll)
  
  # Test opt.params with 'optim'
  res_optim <- HMMoce:::opt.params(
    pars.init = c(2, 0.5, 0.8, 0.8),
    lower.bounds = c(0.5, 0.1, 0.5, 0.5),
    upper.bounds = c(4, 1.0, 0.9, 0.9),
    g = g,
    L = L,
    alg.opt = "optim",
    write.results = FALSE
  )
  
  expect_named(res_optim, c("par", "value", "convergence", "message"))
  expect_equal(length(res_optim$par), 4)
  expect_true(all(res_optim$par >= c(0.5, 0.1, 0.5, 0.5)))
  expect_true(all(res_optim$par <= c(4, 1.0, 0.9, 0.9)))
})
