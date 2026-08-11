test_that("Survr validates its input", {
  expect_error(Survr(c(1, 1, 2), c(5, 3, 4), c(1, 0, 1)),
               "censored time")
  expect_error(Survr(c(1, 1), c(5, 3), c(0, 2)), "0-1")

  x <- Survr(c(1, 1, 2), c(5, 3, 4), c(1, 0, 0))
  expect_s3_class(x, "Survr")
  expect_equal(dim(x), c(3L, 3L))
})

test_that("survfitr rejects unknown estimators", {
  expect_error(survfitr(Survr(id, time, event) ~ 1, data = MMC,
                        type = "nonsense"))
})
