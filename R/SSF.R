`SSF` <-
function (numsim, tss, nbstep = 10, randompart, fixed = c(0, 
    1, 0), intercept=0 ,exgr = NA, exrepl = NA, heteroscedasticity=c("null")) 
{
    if (is.na(exgr)[[1]]) {
        grmin = 2
        grmax = tss/2
    }
    if (is.na(exrepl)[[1]]) {
        remin = 2
        remax = tss/2
    }
    if (!is.na(exgr)[[1]]) {
        grmin = exgr[[1]]
        grmax = exgr[[2]]
    }
    if (!is.na(exrepl)[[1]]) {
        remin = exrepl[[1]]
        remax = exrepl[[2]]
    }
    group = round(seq(grmin, grmax, I((grmax - grmin)/nbstep)))
    repl = ceiling(tss/group)
    mg.r0 = unique(matrix(c(group, repl), ncol = 2))
    mg.r1 = subset(mg.r0, mg.r0[, 2] >= remin)
    mg.r = subset(mg.r1, mg.r1[, 2] <= remax)
    stepvec = c(1:length(mg.r[, 1]))
    VI <- as.numeric(randompart[[1]])
    VS <- as.numeric(randompart[[2]])
    VR <- as.numeric(randompart[[3]])
    if (length(randompart <= 4)) {
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
    
    Hetero <- heteroscedasticity[[1]]
    het <- as.numeric(heteroscedasticity[-1])

    FM <- fixed[[1]]
    FV <- fixed[[2]]
    FE <- fixed[[3]]
    iD <- numeric(length(mg.r[, 1]))
    rp <- numeric(length(mg.r[, 1]))
    ss <- numeric(length(mg.r[, 1]))
    powersl <- numeric(numsim)
    pvalsl <- numeric(numsim)
    slpowestimate <- numeric(length(mg.r[, 1]))
    slpowCIlower <- numeric(length(mg.r[, 1]))
    slpowCIupper <- numeric(length(mg.r[, 1]))
    slpvalestimate <- numeric(length(mg.r[, 1]))
    slpvalCIlower <- numeric(length(mg.r[, 1]))
    slpvalCIupper <- numeric(length(mg.r[, 1]))
    powerint <- numeric(numsim)
    pvalint <- numeric(numsim)
    intpowestimate <- numeric(length(mg.r[, 1]))
    intpowCIlower <- numeric(length(mg.r[, 1]))
    intpowCIupper <- numeric(length(mg.r[, 1]))
    intpvalestimate <- numeric(length(mg.r[, 1]))
    intpvalCIlower <- numeric(length(mg.r[, 1]))
    intpvalCIupper <- numeric(length(mg.r[, 1]))
    kk <- 0
    for (k in stepvec) {
        N <- tss
        for (i in 1:numsim) {
                EF <- rnorm(N, FM, sqrt(FV))
                er <- numeric(length(N))
                if (Hetero=="null") (er <- rnorm(N, intercept, sqrt(VR)))
                if (Hetero=="power") (
                for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*(het[1]+abs(EF[n])^het[2])^2))} )
                if (Hetero=="exp")  (
                for (n in 1:N) {er[n] <- rnorm(1, intercept, sqrt(VR*exp(2*het[1]*EF[n])))} )

            db <- data.frame(ID = rep(1:mg.r[k, 1], mg.r[k, 2])[1:N], 
                obs = 1:N, error = er, EF = EF)
            x <- rmvnorm(mg.r[k, 1], c(0, 0), sigma, method = "svd")
            db$rand.int <- rep(x[, 1], mg.r[k, 2])[1:N]
            db$rand.sl <- rep(x[, 2], mg.r[k, 2])[1:N]
            db$Y <- db$rand.int + (db$rand.sl + FE) * db$EF + 
                db$error
            m.lm <- lm(Y ~ EF, data = db)
            m1.lmer <- lmer(Y ~ EF + (1 | ID), data = db)
		#, control = list(maxIter = 200, msMaxIter = 200))
            pvint <- pchisq(-2 * (logLik(m.lm, REML = TRUE) - 
                logLik(m1.lmer, REML = TRUE))[[1]], 1, lower.tail = FALSE)
            powerint[i] <- pvint <= 0.05
            pvalint[i] <- pvint
            m2.lmer <- lmer(Y ~ EF + (EF | ID), data = db)
		 # , control = list(maxIter = 200,msMaxIter = 200))
            anosl <- anova(m2.lmer, m1.lmer)
            powersl[i] <- anosl[2, "Pr(>Chisq)"] <= 0.05
            pvalsl[i] <- anosl[2, "Pr(>Chisq)"]
        }
        kk <- kk + 1
        iD[kk] <- mg.r[k, 1]
        rp[kk] <- round(N/mg.r[k, 1], digits = 2)
        ss[kk] <- N
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
    sim.sum <- data.frame(nb.ID = iD, nb.repl = rp, N = ss, int.pval = intpvalestimate, 
        CIlow.ipv = intpvalCIlower, CIup.ipv = intpvalCIupper, 
        int.power = intpowestimate, CIlow.ipo = intpowCIlower, 
        CIup.ipo = intpowCIupper, sl.pval = slpvalestimate, CIlow.slpv = slpvalCIlower, 
        CIup.slpv = slpvalCIupper, sl.power = slpowestimate, 
        CIlow.slpo = slpowCIlower, CIup.slpo = slpowCIupper)
    class(sim.sum) <- c("SSF", "data.frame")
    sim.sum
}
