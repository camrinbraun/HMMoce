test_that("bin_TempTS bins and summarizes depth-temperature data correctly", {
  dates <- as.POSIXct(c("2020-01-01 00:00:00", "2020-01-02 00:00:00"), tz = "UTC")
  
  # Day 1 has 3 distinct depth bins (>= 3 required by filter)
  # Day 2 has only 2 distinct depth bins (should be excluded)
  ts <- data.frame(
    Date = as.POSIXct(c(
      "2020-01-01 01:00:00",
      "2020-01-01 02:00:00",
      "2020-01-01 03:00:00",
      "2020-01-01 04:00:00",
      "2020-01-02 01:00:00",
      "2020-01-02 02:00:00"
    ), tz = "UTC"),
    Depth = c(0, 10, 20, 22, 10, 20),
    Temperature = c(20.5, 18.2, 15.0, 14.8, 17.0, 16.0)
  )
  
  res <- bin_TempTS(ts, out_dates = dates, bin_res = 10)
  
  expect_s3_class(res, "data.frame")
  expect_named(res, c("Depth", "nrecs", "MeanTemp", "MinTemp", "MaxTemp", "Date", "bin", "MeanPDT"))
  
  # Only Day 1 (time_idx = 1) should remain
  expect_equal(unique(res$Date), dates[1])
  expect_equal(nrow(res), 3)
  
  # Check binning with bin_res = 10 (0, 10, and 20 [20 and 22 rounded to 20])
  expect_equal(res$Depth, c(0, 10, 20))
  expect_equal(res$nrecs, c(1, 1, 2))
  expect_equal(res$MinTemp[res$Depth == 20], 14.8)
  expect_equal(res$MaxTemp[res$Depth == 20], 15.0)
  expect_equal(res$MeanPDT[res$Depth == 20], (14.8 + 15.0) / 2)
})

test_that("bin_TempTS drops NA values and handles missing depth/temperature", {
  dates <- as.POSIXct(c("2020-01-01 00:00:00"), tz = "UTC")
  ts <- data.frame(
    Date = as.POSIXct(rep("2020-01-01 01:00:00", 5), tz = "UTC"),
    Depth = c(0, 10, 20, NA, 30),
    Temperature = c(20, 18, 15, 12, NA)
  )
  
  res <- bin_TempTS(ts, out_dates = dates, bin_res = 10)
  expect_equal(res$Depth, c(0, 10, 20))
  expect_equal(nrow(res), 3)
})
