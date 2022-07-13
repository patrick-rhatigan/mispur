test_that("pos works", {
  expect_equal(pos(-2:2), c(0,0,0,1,2))
})

test_that("neg works", {
  expect_equal(neg(-2:2), c(2,1,0,0,0))
})

test_that("pmax1 works", {
  expect_equal(pmax1(-2:2), c(1,1,1,1,2))
})

test_that("phi works", {
  expect_equal(phi(-2:2), c(0,0,0,0,Inf))
})

test_that("chi works", {
  expect_equal(chi(c(5, -6, 0, 1, -4), -2:2), c(-2, 6, 0, 0, 2))
})
