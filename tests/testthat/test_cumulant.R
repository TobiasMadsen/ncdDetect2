library(ncdDetect2)
context("cumulant")

test_that("cumulant one {0,1} variable",{
  res <- cumulantDerivatives(0.5, x = c(1,1), p = c(0.3,0.7), s = c(1,0))

  expect_equal(res[1], log(0.7+0.3*exp(0.5)))
  expect_equal(res[2], (0.3*exp(0.5)) / (0.7+0.3*exp(0.5)) )
  expect_equal(res[3], 0.3*0.7*exp(0.5) / (0.7+0.3*exp(0.5))**2 )
})

test_that("cumulant four {0,1} variable",{
  res <- cumulantDerivatives(0.5, x = rep(1:4,each = 2), p = rep(c(0.3,0.7), 4), s = rep(c(1,0), 4) )

  expect_equal(res[1], 4*log(0.7+0.3*exp(0.5)))
  expect_equal(res[2], 4*(0.3*exp(0.5)) / (0.7+0.3*exp(0.5)) )
  expect_equal(res[3], 4*0.3*0.7*exp(0.5) / (0.7+0.3*exp(0.5))**2 )
})

test_that("cumulant one {0,1,2,3} variable", {
  res <- cumulantDerivatives(0.5, x = rep(1,4), p = 1:4/10, s = 0:3)

  phi0 <- sum(1:4*exp(0.5*0:3)/10)
  phi1 <- sum(1:4*0:3*exp(0.5*0:3)/10)
  phi2 <- sum(1:4*0:3*0:3*exp(0.5*0:3)/10)
  expect_equal(res[1], log(phi0) )
  expect_equal(res[2], phi1 / phi0)
  expect_equal(res[3], (phi0*phi2 - phi1*phi1)/phi0**2 )
})
