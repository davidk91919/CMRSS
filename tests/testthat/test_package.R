################################################################################
# Tests for the CMRSS package
################################################################################
## see https://testthat.r-lib.org/
library("testthat")

context("CMRSS Functions")

test_that("xBal univariate desriptive means agree w/ lm",{
    set.seed(20160406)
    n <- 7
     dat <- data.frame(x1=rnorm(n), x2=rnorm(n),
                        s=rep(c("a", "b"), c(floor(n/2), ceiling(n/2)))
                        )
     dat = transform(dat, z=as.numeric( (x1+x2+rnorm(n))>0 ) )

     lm1 <- lm(x1~z, data=dat)
     xb1 <- xBalance(z~x1, strata = list(`Unstrat` = NULL, s = ~s), data=dat, report=c("adj.mean.diffs"))
     expect_equal(xb1$results["x1", "adj.diff", "Unstrat"], coef(lm1)["z"], check.attributes=F)

     lm2a <- lm(x1~z+s, data=dat)
     expect_equivalent(xb1$results["x1", "adj.diff", "s"], coef(lm2a)[["z"]])
})


