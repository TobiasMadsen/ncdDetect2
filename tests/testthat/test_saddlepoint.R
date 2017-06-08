library(ncdDetect2)
library(data.table)

context("saddlepoint")

test_that("Binomial distribution",{
  sp <- saddlepoint(51, dat = data.table(x = rep(1:100, each = 2),
                                         y = rep(c(1,0),100),
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 50, size = 100, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space

  sp <- saddlepoint(61, dat = data.table(x = rep(1:100, each = 2),
                                         y = rep(c(1,0),100),
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
})

test_that("Unsorted x",{
  set.seed(1)

  dat <- data.table(x = rep(1:100, each = 2),
                    y = rep(c(1,0),100),
                    probability = rep(c(0.4,0.6), 100))
  dat <- dat[sample(200),]

  sp <- saddlepoint(61, dat = dat)
  pb <- pbinom(q = 60, size = 100, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp-pb), log(1.05) ) # At most 5 % difference in probability space
})

test_that("Saddlepoint in 0",{
  # Saddlepoint approximation should return 1 when the score is 0
  sp <- saddlepoint(0, dat = data.table(x = rep(1:100, each = 2),
                                        y = rep(c(1,0),100),
                                        probability = rep(c(0.4,0.6), 100)), log = F)

  expect_gt(sp, 0.99999)

  sp <- saddlepoint(0, dat = data.table(x = rep(1:100, each = 2),
                                        y = rep(c(1,0),100),
                                        probability = rep(c(0.4,0.6), 100)), log = T)

  expect_gt(sp, -.00001)
})

test_that("Saddlepoint high p-values",{
  sp <- saddlepoint(11, dat = data.table(x = rep(1:100, each = 2),
                                         y = rep(c(1,0),100),
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 10, size = 100, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp-pb), log(1.05) )

  sp <- saddlepoint(21, dat = data.table(x = rep(1:100, each = 2),
                                         y = rep(c(1,0),100),
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 20, size = 100, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp-pb), log(1.05))

  sp <- saddlepoint(31, dat = data.table(x = rep(1:100, each = 2),
                                         y = rep(c(1,0),100),
                                         probability = rep(c(0.4,0.6), 100)))
  pb <- pbinom(q = 30, size = 100, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp-pb), log(1.05))
})

test_that("Saddlepoint outside range",{
  sp <- saddlepoint(101, dat = data.table(x = rep(1:100, each = 2),
                                         y = rep(c(1,0),100),
                                         probability = rep(c(0.4,0.6), 100)), log = F)

  expect_lt(sp, 1e-300)
})

test_that("Saddlepoint t == min_score + lattice",{
  sp <- saddlepoint(1, dat = data.table(x = rep(1:30, each = 2),
                                        y = rep(c(1,0),30),
                                        probability = rep(c(0.4,0.6), 30)))
  pb <- pbinom(q = 0, size = 30, prob = 0.4, lower.tail = F, log.p = T)

  expect_lt(abs(sp - pb), log(1.05))
})

test_that("Avoid too large first jumps in Newton Raphson",{
  dat = data.table(x = rep(1:300, each = 4),
                   y = c(rep(c(177,177,177,0), 100),
                         rep(c(170,170,170,0), 100),
                         rep(c(163,163,163,0), 100)),
                   probability = c(rep(c(2e-8,2e-8,2e-8,1-6e-8), 100),
                                   rep(c(4e-8,4e-8,4e-8,1-1.2e-7), 100),
                                   rep(c(8e-8,8e-8,8e-8,1-2.4e-7), 100)) )

  sp <- saddlepoint(281, dat = dat, log = F)
})
