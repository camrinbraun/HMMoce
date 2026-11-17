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
  r <- terra::rast(
    xmin = -80, xmax = -70, ymin = 30, ymax = 40,
    resolution = 0.5, crs = "+proj=longlat +datum=WGS84"
  )
  terra::values(r) <- 1
  g_ras <- setup.grid.raster(r)
  expect_named(g_ras, c("lon", "lat", "dlo", "dla"))
  expect_equal(g_ras$dlo, 0.5)
  expect_equal(g_ras$dla, 0.5)
  expect_equal(dim(g_ras$lon), c(terra::nrow(r), terra::ncol(r)))
})

test_that("resample.grid aligns list of rasters to reference grid", {
  r1 <- terra::rast(xmin = -75, xmax = -70, ymin = 30, ymax = 35, res = 0.25)
  terra::values(r1) <- runif(terra::ncell(r1))
  
  r2 <- terra::rast(xmin = -75, xmax = -70, ymin = 30, ymax = 35, res = 0.5)
  terra::values(r2) <- runif(terra::ncell(r2))
  
  res <- resample.grid(list(r1 = r1, r2 = r2), L.res = r1)
  expect_named(res, c("L.rasters", "L.mle.res", "g", "g.mle"))
  expect_equal(terra::res(res$L.rasters[[2]])[1], 0.25)
  expect_equal(dim(res$g$lon), dim(res$g$lat))
})
