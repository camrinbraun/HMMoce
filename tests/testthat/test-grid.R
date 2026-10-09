test_that("setup.grid and setup.grid.raster construct consistent grid structures", {
  # 1. setup.grid
  locs <- data.frame(
    Longitude = c(-75, -70),
    Latitude = c(30, 35)
  )
  g_quarter <- setup.grid(locs, res = "quarter")
  expect_named(g_quarter, c("lon", "lat", "dlo", "dla"))
  expect_equal(dim(g_quarter$lon), dim(g_quarter$lat))
  expect_true(g_quarter$dlo > 0.2 && g_quarter$dlo < 0.3)
  expect_true(g_quarter$dla > 0.2 && g_quarter$dla < 0.3)
  
  # 2. setup.grid.raster
  r <- raster::raster(
    xmn = -80, xmx = -70, ymn = 30, ymx = 40,
    resolution = 0.5, crs = "+proj=longlat +datum=WGS84"
  )
  raster::values(r) <- 1
  g_ras <- setup.grid.raster(r)
  expect_named(g_ras, c("lon", "lat", "dlo", "dla"))
  expect_equal(g_ras$dlo, 0.5)
  expect_equal(g_ras$dla, 0.5)
  expect_equal(dim(g_ras$lon), c(raster::nrow(r), raster::ncol(r)))
})

test_that("resample.grid aligns list of rasters to reference grid", {
  r1 <- raster::raster(xmn = -75, xmx = -70, ymn = 30, ymx = 35, res = 0.25)
  raster::values(r1) <- runif(raster::ncell(r1))
  
  r2 <- raster::raster(xmn = -75, xmx = -70, ymn = 30, ymx = 35, res = 0.5)
  raster::values(r2) <- runif(raster::ncell(r2))
  
  res <- resample.grid(list(r1 = r1, r2 = r2), L.res = r1)
  expect_named(res, c("L.rasters", "L.mle.res", "g", "g.mle"))
  expect_equal(raster::res(res$L.rasters[[2]])[1], 0.25)
  expect_equal(dim(res$g$lon), dim(res$g$lat))
})
