test_that("make180 and make360 convert longitudes properly", {
  lons_neg <- c(-75, -180, 0, 180)
  expect_equal(make360(lons_neg), c(285, 180, 0, 180))
  
  lons_pos <- c(285, 180, 0, 360)
  expect_equal(make180(lons_pos), c(-75, 180, 0, 0))
  
  # Idempotence where appropriate
  expect_equal(make180(make360(-45)), -45)
})

test_that("reverse transposes and flips matrix rows", {
  m <- matrix(1:4, nrow = 2, ncol = 2)
  rev_m <- HMMoce:::reverse(m)
  expect_equal(dim(rev_m), c(2, 2))
  expect_equal(rev_m, t(m)[, 2:1])
})

test_that("repmat and meshgrid behave like matrix tiling", {
  a <- matrix(1:4, nrow = 2)
  rep_a <- repmat(a, 2, 3)
  expect_equal(dim(rep_a), c(4, 6))
  
  mg <- meshgrid(1:3, 1:4)
  expect_equal(dim(mg$X), c(4, 3))
  expect_equal(dim(mg$Y), c(4, 3))
})

test_that("findDateFormat correctly identifies supported date formats", {
  expect_equal(findDateFormat("2015-01-01 05:30:17"), "%Y-%m-%d %H:%M:%S")
  expect_equal(findDateFormat("05/15/2020"), "%m/%d/%Y")
  expect_equal(findDateFormat("05/15/20 14:30"), "%m/%d/%y %H:%M")
})

test_that("km.per.gridunit calculates distances correctly", {
  limits <- list(lonmin = -80, lonmax = -70, latmin = 25, latmax = 35)
  res <- km.per.gridunit(limits, res = 0.1)
  expect_equal(length(res), 2)
  expect_true(res[1] <= res[2])
  expect_true(all(res > 0))
})

test_that("track distance functions calculate distances accurately", {
  tr <- data.frame(
    lon = c(-75, -75, -74),
    lat = c(35, 36, 36)
  )
  dates <- as.POSIXct(c("2020-01-01", "2020-01-02", "2020-01-03"), tz = "UTC")
  
  step_d <- track.stepdist(tr)
  expect_equal(length(step_d), nrow(tr) - 1)
  expect_true(all(step_d > 0))
  
  cum_d <- track.cumdist(tr)
  expect_true(is.numeric(cum_d))
  expect_true(cum_d > 0)
  
  avg_d <- track.avgdist(tr, dates)
  expect_true(is.numeric(avg_d))
  expect_true(avg_d > 0)
})
