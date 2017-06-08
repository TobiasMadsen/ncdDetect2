
context("estimate overdispersion")

test_that("IO",{
  od <- estimateOverDispersion(mutLoads$observed[1:1000], mutLoads$expected[1:1000], mutLoads$length[1:1000] * mutLoads$n_samples[1:1000])
  expect_is(od, "numeric")
  expect_equal(names(od), c("Estimate", "CI_5%", "CI_95%"))
})

test_that("Simulation 0.3",{
  set.seed(1)

  od_true <- 0.3
  length <- rpois(400, 1000)*2000
  prob <- rbeta(400, 1.5, 3e5)
  obs <- VGAM::rbetabinom(400, length, prob, od_true^2*prob/(1-prob))
  exp <- length * prob
  od <- estimateOverDispersion(obs, exp, length)

  expect_gt( od['CI_95%'], od_true)
  expect_lt( od['CI_5%'], od_true)
  expect_lt( abs(od['Estimate']-od_true), 0.05 )
})

test_that("Simulation 0.6",{
  set.seed(1)

  od_true <- 0.6
  length <- rpois(400, 1000)*2000
  prob <- rbeta(400, 1.5, 3e5)
  obs <- VGAM::rbetabinom(400, length, prob, od_true^2*prob/(1-prob))
  exp <- length * prob
  od <- estimateOverDispersion(obs, exp, length)

  expect_gt( od['CI_95%'], od_true)
  expect_lt( od['CI_5%'], od_true)
  expect_lt( abs(od['Estimate']-od_true), 0.05 )
})

