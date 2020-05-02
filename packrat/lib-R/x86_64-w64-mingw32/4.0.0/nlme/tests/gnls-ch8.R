## Really just one section from ../inst/scripts/ch08.R :

library(nlme)

(fm1Lis <- nlsList(rate ~ SSasympOff(pressure, Asym, lrc, c0) | QB,
                   data = Dialyzer))
(iLis <- intervals(fm1Lis))
## TODO: stopifnot(all.equal(..)) ?

fm1Dial.gnls <- gnls(rate ~ SSasympOff(pressure, Asym, lrc, c0),
                     data = Dialyzer,
                     params = list(Asym + lrc ~ QB, c0 ~ 1),
                     start = c(53.6, 8.6, 0.51, -0.26, 0.225))
summary(fm1Dial.gnls)

## Modified Data (==> rename it !) ----------------------------------
DialyzerM <- Dialyzer
DialyzerM$QBcontr <- 2 * (Dialyzer$QB == 300) - 1
fm1Dial.nls <-
  nls(rate ~ SSasympOff(pressure, Asym.Int + Asym.QB * QBcontr,
                        lrc.Int + lrc.QB * QBcontr, c0),
      data = DialyzerM,
      start = c(Asym.Int = 53.6, Asym.QB = 8.6, lrc.Int = 0.51,
                lrc.QB = -0.26, c0 = 0.225))
summary(fm1Dial.nls)
logLik(fm1Dial.nls)
plot(fm1Dial.gnls, resid(.) ~ pressure, abline = 0)
fm2Dial.gnls <- update(fm1Dial.gnls,
                       weights = varPower(form = ~ pressure))
anova(fm1Dial.gnls, fm2Dial.gnls)
ACF(fm2Dial.gnls, form = ~ 1 | Subject)
plot(ACF(fm2Dial.gnls, form = ~ 1 | Subject), alpha = 0.05)
fm3Dial.gnls <-
    update(fm2Dial.gnls, corr = corAR1(0.716, form = ~ 1 | Subject))
fm3Dial.gnls
(im3 <- intervals(fm3Dial.gnls))
(a23 <- anova(fm2Dial.gnls, fm3Dial.gnls))

anoC <-
    cbind(fm2Dial.gnls =
              c(df=7, AIC= 748.4749, BIC= 769.0664, logLik=-367.2375,
                L.Ratio=NA, "p-value"=NA),
          fm3Dial.gnls=
              c(8, 661.0424, 684.5756, -322.5212, 89.4325, 3.e-21))
                                        # NB: exact p-value irrelevant
stopifnot(
    all.equal(anoC, t(data.matrix(a23)[,rownames(anoC)]), tol = 7e-7)
)
