`PAMM` <-
function (numsim, group, repl, randompart, fixed = c(0, 1, 0),intercept=0) 
{
    VI <- as.numeric(randompart[[1]])
    VS <- as.numeric(randompart[[2]])
    VR <- as.numeric(randompart[[3]])
    if (length(randompart) <= 4) {
        CorIS <- as.numeric(randompart[[4]])
        CovIS <- CorIS * sqrt(VI) * sqrt(VS)
    }
    else {
        if (randompart[[5]] == "cor") {
            CorIS <- as.numeric(randompart[[4]])
            CovIS <- CorIS * sqrt(VI) * sqrt(VS)
        }
        if (randompart[[5]] == "cov") {
            CovIS <- as.numeric(randompart[[4]])
        }
    }
    sigma <- matrix(c(VI, CovIS, CovIS, VS), ncol = 2)
    FM <- fixed[[1]]
    FV <- fixed[[2]]
    FE <- fixed[[3]]
    iD <- numeric(length(repl) * length(group))
    rp <- numeric(length(repl) * length(group))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(repl) * length(group))
    slpowCIlower <- numeric(length(repl) * length(group))
    slpowCIupper <- numeric(length(repl) * length(group))
    slpvalestimate <- numeric(length(repl) * length(group))
    slpvalCIlower <- numeric(length(repl) * length(group))
    slpvalCIupper <- numeric(length(repl) * length(group))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(repl) * length(group))
    intpowCIlower <- numeric(length(repl) * length(group))
    intpowCIupper <- numeric(length(repl) * length(group))
    intpvalestimate <- numeric(length(repl) * length(group))
    intpvalCIlower <- numeric(length(repl) * length(group))
    intpvalCIupper <- numeric(length(repl) * length(group))
    kk <- 0
    for (k in group) {
        for (r in repl) {
            N <- k * r
            for (i in 1:numsim) {
                er <- rnorm(N, intercept, sqrt(VR))
                EF <- rnorm(N, FM, sqrt(FV))
                db <- data.frame(ID = rep(1:k, r), obs = 1:N, 
                  error = er, EF = EF)
                x <- rmvnorm(k, c(0, 0), sigma, method = "svd")
                db$rand.int <- rep(x[, 1], r)
                db$rand.sl <- rep(x[, 2], r)
                db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + 
                  db$error
                m.lm <- lm(Y ~ EF, data = db)
                m1.lmer <- lmer(Y ~ EF + (1 | ID), data = db, 
                  control = list(maxIter = 200, msMaxIter = 200))
                pvint <- pchisq(-2 * (logLik(m.lm, REML = TRUE) - 
                  logLik(m1.lmer, REML = TRUE))[[1]], 1, lower.tail = FALSE)
                powerint[i] <- pvint <= 0.05
                pvalint[i] <- pvint
                m2.lmer <- lmer(Y ~ EF + (EF | ID), data = db, 
                  control = list(maxIter = 200, msMaxIter = 200))
                anosl <- anova(m2.lmer, m1.lmer)
                powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
                pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
            }
            kk <- kk + 1
            iD[kk] <- k
            rp[kk] <- r
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
    }
    sim.sum <- data.frame(nb.ID = iD, nb.repl = rp, int.pval = intpvalestimate, 
        CIlow.ipv = intpvalCIlower, CIup.ipv = intpvalCIupper, 
        int.power = intpowestimate, CIlow.ipo = intpowCIlower, 
        CIup.ipo = intpowCIupper, sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
        CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
        CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper)
    class(sim.sum) = c("PAMM", "data.frame")
    sim.sum
}
