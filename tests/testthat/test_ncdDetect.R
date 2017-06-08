context("ncdDetect")

test_that("IO",{

  predictions <- matrix(c(0.9,0.1), 100, 2, byrow = T)
  scores <- matrix(c(0,1), 100, 2, byrow = T)

  expect_type(ncdDetect(predictions, scores, 30, 0.1, N = 3), 'double')
})

test_that("Agrees with binomial",{
  predictions <- matrix(c(0.99,0.01), 1000, 2, byrow = T)
  scores <- matrix(c(0,1), 1000, 2, byrow = T)

  # N=1 should correspond to a binomial model
  ncd <- ncdDetect(predictions, scores, 30, 1e-15, N = 1)
  bin <- pbinom(29, size = 1000, prob = 0.01, lower.tail = F, log.p = T)

  expect_lt( abs(ncd - bin), log(2))
})

test_that("Agrees with betabinomial",{
  predictions <- matrix(c(0.99,0.01), 1000, 2, byrow = T)
  scores <- matrix(c(0,1), 1000, 2, byrow = T)

  od <- 0.1
  ncd <- ncdDetect(predictions, scores, 30, od, N = 11)
  bin <- log( sum( VGAM::dbetabinom(30:1000, size = 1000, prob = 0.01, rho = od^2 * 0.01 / (1-0.01)) ) )
  expect_lt( abs(ncd - bin), log(2))

  od <- 0.3
  ncd <- ncdDetect(predictions, scores, 50, od, N = 11)
  bin <- log( sum( VGAM::dbetabinom(50:1000, size = 1000, prob = 0.01, rho = od^2 * 0.01 / (1-0.01)) ) )
  expect_lt( abs(ncd -bin), log(2))
})
