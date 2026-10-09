test_that("gausskern and variants generate valid probability kernels", {
  # 1. gausskern.pg
  k_pg <- HMMoce:::gausskern.pg(siz = 9, sigma = 2, muadv = 0)
  expect_equal(dim(k_pg), c(9, 9))
  expect_equal(sum(k_pg), 1, tolerance = 1e-6)
  expect_true(all(k_pg >= 0))
  
  # Peak should be at the center (5, 5)
  expect_equal(which(k_pg == max(k_pg), arr.ind = TRUE), matrix(c(5, 5), nrow = 1, dimnames = list(NULL, c("row", "col"))))

  # 2. gausskern.nostd
  k_nostd <- HMMoce:::gausskern.nostd(siz = 9, sigma = 2, muadv = 0)
  expect_equal(dim(k_nostd), c(9, 9))
  expect_true(sum(k_nostd) > 0)
  
  # 3. gausskern.isotrop
  k_iso <- HMMoce:::gausskern.isotrop(siz = 9, sigma = 2)
  expect_true(is.matrix(k_iso))
  expect_equal(dim(k_iso), c(9, 9))
  
  # 4. mask.K
  masked <- HMMoce:::mask.K(k_pg)
  expect_equal(dim(masked), dim(k_pg))
  expect_true(all(masked >= 0))
})
