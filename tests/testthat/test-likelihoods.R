test_that("calc.bathy generates valid likelihood raster brick", {
  dates <- as.POSIXct(c("2020-01-01 00:00:00", "2020-01-02 00:00:00"), tz = "UTC")
  
  mmd <- data.frame(
    Date = dates,
    MaxDepth = c(50, 150)
  )
  
  bathy <- raster::raster(xmn = -75, xmx = -70, ymn = 30, ymx = 35, res = 0.5)
  raster::values(bathy) <- seq(10, 500, length.out = raster::ncell(bathy))
  
  l_bathy <- calc.bathy(mmd = mmd, bathy.grid = bathy, dateVec = dates, focalDim = 3)
  expect_s4_class(l_bathy, "RasterBrick")
  expect_equal(raster::nlayers(l_bathy), 2)
  expect_true(raster::cellStats(l_bathy[[1]], max) > 0)
})

test_that("setup.locs.grid and calc.lightloc generate longitude/elliptical likelihoods", {
  sp.lim <- list(lonmin = -75, lonmax = -70, latmin = 30, latmax = 35)
  locs.grid <- setup.locs.grid(sp.lim, res = "quarter")
  
  dates <- as.POSIXct(c("2020-01-01 00:00:00", "2020-01-02 00:00:00"), tz = "UTC")
  lightloc <- data.frame(
    Date = dates,
    Longitude = c(-73, -72),
    Error.Semi.minor.axis = c(50000, 50000)
  )
  
  l_light <- calc.lightloc(lightloc = lightloc, locs.grid = locs.grid, dateVec = dates, errEll = FALSE)
  expect_s4_class(l_light, "RasterBrick")
  expect_equal(raster::nlayers(l_light), 2)
  expect_true(raster::cellStats(l_light[[1]], max) > 0)
})

test_that("make.L combines rasters into normalized 3D array", {
  dates <- as.POSIXct(c("2020-01-01 00:00:00", "2020-01-02 00:00:00"), tz = "UTC")
  
  r1 <- raster::raster(xmn = -75, xmx = -70, ymn = 30, ymx = 35, res = 0.5)
  raster::values(r1) <- 1
  b1 <- raster::brick(r1, r1)
  
  bathy <- raster::raster(xmn = -75, xmx = -70, ymn = 30, ymx = 35, res = 0.5)
  raster::values(bathy) <- 200
  
  iniloc <- data.frame(
    year = c(2020, 2020),
    month = c(1, 1),
    day = c(1, 2),
    lon = c(-74, -71),
    lat = c(31, 34)
  )
  
  L <- make.L(
    ras.list = list(bathy = b1),
    iniloc = iniloc,
    dateVec = dates,
    maxDepth = c(50, 50),
    bathy = bathy
  )
  
  expect_true(is.array(L))
  expect_equal(length(dim(L)), 3)
  expect_equal(dim(L)[1], 2) # time dimension first
  expect_true(all(L >= 0))
})
