`EAMM` <-
function (numsim, group, repl, fixed = c(0, 1, 0), VI = seq(0.05, 
    0.95, 0.05), VS = seq(0.05, 0.5, 0.05), CoIS = 0, relIS = "cor") 
{
    FM <- fixed[[1]]
    FV <- fixed[[2]]
    FE <- fixed[[3]]
    vgi <- numeric(length(VI) * length(VS))
    vgs <- numeric(length(VI) * length(VS))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(VI) * length(VS))
    slpowCIlower <- numeric(length(VI) * length(VS))
    slpowCIupper <- numeric(length(VI) * length(VS))
    slpvalestimate <- numeric(length(VI) * length(VS))
    slpvalCIlower <- numeric(length(VI) * length(VS))
    slpvalCIupper <- numeric(length(VI) * length(VS))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(VI) * length(VS))
    intpowCIlower <- numeric(length(VI) * length(VS))
    intpowCIupper <- numeric(length(VI) * length(VS))
    intpvalestimate <- numeric(length(VI) * length(VS))
    intpvalCIlower <- numeric(length(VI) * length(VS))
    intpvalCIupper <- numeric(length(VI) * length(VS))
    kk <- 0
    for (k in VI) {
        for (r in VS) {
            N <- group * repl
            VR <- 1 - k
            if (VR >= 0) {
                for (i in 1:numsim) {
                  if (relIS == "cor") {
                    CovIS <- CoIS * sqrt(k) * sqrt(r)
                  }
                  if (relIS == "cov") {
                    CovIS <- CoIS
                  }
                  sigma <- matrix(c(k, CovIS, CovIS, r), ncol = 2)
                  er <- rnorm(N, 0, sqrt(VR))
                  X <- sort(rep(c(1:repl), group))
                  EF <- rnorm(N, FM, sqrt(FV))
                  db <- data.frame(ID = rep(1:group, repl), obs = 1:N, 
                    error = er, X = X, EF = EF)
                  x <- rmvnorm(group, c(0, 0), sigma, method = "chol")
                  db$rand.int <- rep(x[, 1], repl)
                  db$rand.sl <- rep(x[, 2], repl)
                  db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + 
                    db$error
                  m.lm <- lm(Y ~ EF, data = db)
                  m1.lmer <- lmer(Y ~ EF + (1 | ID), data = db, 
                    control = list(maxIter = 200, msMaxIter = 200))
                  pvint <- pchisq(-2 * (logLik(m.lm, REML = TRUE) - 
                    logLik(m1.lmer, REML = TRUE))[[1]], 1, lower.tail = FALSE)
                  powerint[i] <- pvint <= 0.05
                  pvalint[i] <- pvint
                  if (r > 0) {
                    m2.lmer <- lmer(Y ~ EF + (EF | ID), data = db, 
                      control = list(maxIter = 200, msMaxIter = 200))
                    anosl <- anova(m2.lmer, m1.lmer)
                    powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
                    pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
                  }
                  else {
                    powersl[i] <- 0
                    pvalsl[i] <- 1
                  }
                }
                kk <- kk + 1
                vgi[kk] <- k
                vgs[kk] <- r
                slCIpow <- ci(powersl)
                slpowestimate[kk] <- slCIpow["Estimate"]
                slpowCIlower[kk] <- slCIpow["CI lower"]
                slpowCIupper[kk] <- slCIpow["CI upper"]
                slCIpval <- ci(pvalsl)
                slpvalestimate[kk] <- slCIpval["Estimate"]
                slpvalCIlower[kk] <- slCIpval["CI lower"]
                slpvalCIupper[kk] <- slCIpval["CI upper"]
                intCIpow <- ci(powerint)
                intpowestimate[kk] <- intCIpow["Estimate"]
                intpowCIlower[kk] <- intCIpow["CI lower"]
                intpowCIupper[kk] <- intCIpow["CI upper"]
                intCIpval <- ci(pvalint)
                intpvalestimate[kk] <- intCIpval["Estimate"]
                intpvalCIlower[kk] <- intCIpval["CI lower"]
                intpvalCIupper[kk] <- intCIpval["CI upper"]
            }
            else {
                kk <- kk + 1
                vgi[kk] <- k
                vgs[kk] <- r
                slpowestimate[kk] <- NA
                slpowCIlower[kk] <- NA
                slpowCIupper[kk] <- NA
                slpvalestimate[kk] <- NA
                slpvalCIlower[kk] <- NA
                slpvalCIupper[kk] <- NA
                intpowestimate[kk] <- NA
                intpowCIlower[kk] <- NA
                intpowCIupper[kk] <- NA
                intpvalestimate[kk] <- NA
                intpvalCIlower[kk] <- NA
                intpvalCIupper[kk] <- NA
            }
        }
    }
    sim.sum <- data.frame(nb.ID = rep(group, length(VI) * length(VS)), 
        nb.repl = rep(repl, length(VI) * length(VS)), VI = vgi, 
        VS = vgs, int.pval = intpvalestimate, CIlow.ipv = intpvalCIlower, 
        CIup.ipv = intpvalCIupper, int.power = intpowestimate, 
        CIlow.ipo = intpowCIlower, CIup.ipo = intpowCIupper, 
        sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
        CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
        CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper)
    class(sim.sum) <- c("EAMM", "data.frame")
    sim.sum
}
