require(SuperLearner)

# Continuous-time generation, discretized at varying resolutions (AE6 / Reviewer 2).
# Event and censoring times are drawn from Weibull proportional-hazards models with
# time-fixed continuous confounders L1, L2 (see source/generateDataCont.R), then each
# data set is discretized onto K = 4, 8, 16, 30 intervals over [0, tau]. The estimator
# is run at every resolution and the ATE is read at the final interval (= horizon tau).
# The true continuous-time ATE is a fixed constant: true_ate_cont(2e6)["ATE"] ~ 0.20.
data.fcn <- "continuous"
source(paste0("../../../source/generateDataCont.R"), local = TRUE)
source(paste0("../../../source/SurvTMLE.R"), local = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!dir.exists("result")) {
  dir.create("result")
}

ite <- 250
n <- 2000
Ks <- c(4, 8, 16, 30)
seed_base <- 20250616   # per-iteration seed; identical across methods (paired design)

for (i in 1:ite) {
  set.seed(seed_base + i)
  fac <- gen_cont_factual(n)
  rows <- list()
  for (K in Ks) {
    d  <- discretize_cont(fac, K)
    Lv <- replicate(K, c("L1", "L2"), simplify = FALSE)
    res <- surv.TMLE(dat = d, Yvar = paste0("Y", 1:K), Cvar = paste0("Censor", 1:K),
                     Avar = "A", Lvar = Lv, L0var = NULL, lookback = 1,
                     Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                     gbound = 0.005, V = 5, MSM.form = ~A+as.factor(time), Print = FALSE)
    r <- K + 1   # interval K (= horizon tau); row 1 is baseline
    rows[[length(rows) + 1]] <- data.frame(K = K, ATE = -res$ATE[r], SE = res$SE.ATE[r],
                                           lower = -res$UL.ATE[r], upper = -res$LL.ATE[r])
  }
  saveRDS(do.call(rbind, rows), file = paste0("result/", data.fcn, "_para_", i, ".Rds"))
}
