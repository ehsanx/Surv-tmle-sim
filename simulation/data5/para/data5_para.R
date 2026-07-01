require(SuperLearner)

# Scenario 5: time-varying treatment effect (extends Scenario 2; the conditional
# treatment effect on the outcome hazard declines across follow-up intervals).
# See source/generateData5.R (default beta = c(1.5, 1.0, 0.5, 0.25), a waning effect).
data.fcn <- "data5"
source(paste0("../../../source/generateData5.R"), local = TRUE)
source(paste0("../../../source/SurvTMLE.R"), local = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!dir.exists("result")) {
  dir.create("result")
}

ite <- 1000 # 1000
n <- 10000 # 5000
seed_base <- 20250616   # per-iteration seed; identical across methods (paired design)

for (i in 1:ite) {
  set.seed(seed_base + i)
  FullData <- generateData5(n, ntime = 4)

  ## True survival ATE trajectory (t = 1..4)
  true1 <- with(FullData, (1 - colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1))))
  true0 <- with(FullData, (1 - colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0))))
  trueATE <- true1 - true0

  dat <- FullData

  res.para <- surv.TMLE(dat = dat,
                        Yvar = c("Y1", "Y2", "Y3", "Y4"),
                        Cvar = c("Censor1", "Censor2", "Censor3", "Censor4"),
                        Avar = "A",
                        Lvar = list("L1", "L2", "L3", "L4"),
                        L0var = "L",
                        lookback = 1,
                        Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                        gbound = 0.005, V = 5,
                        MSM.form = ~A+time, Print = FALSE)

  ## Save the full per-interval ATE trajectory (rows t = 1..4); surv.TMLE returns
  ## length-5 vectors (row 1 = baseline t=0), and ATE is negated to give treated-minus-control.
  result.para <- data.frame(t = 1:4,
                            true  = trueATE[1:4],
                            ATE   = -res.para$ATE[2:5],
                            SE    =  res.para$SE.ATE[2:5],
                            lower = -res.para$UL.ATE[2:5],
                            upper = -res.para$LL.ATE[2:5])
  saveRDS(result.para, file = paste0("result/", data.fcn, "_para_", i, ".Rds"))
}
