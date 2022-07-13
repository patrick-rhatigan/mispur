test_that("mhat_fun works", {
  k1 <- k2 <- 5
  mhat_func <- gen_mhat_fun(k1 = k1, k2 = k2)
  expect_equal(mhat_fun(1:10, 5, 0.5), c(-8, -6, -4, -2, 0, -2, -4, -6, -8))
})
