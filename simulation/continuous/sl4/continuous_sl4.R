require(SuperLearner)

# Continuous-time generation, discretized at varying resolutions
# -- SL[4] = {SL.mean, SL.glm, SL.glm.interaction, SL.earth}. Slowest at K=30 with MARS.
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
seed_base <- 20250616

for (i in 1:ite) {
  set.seed(seed_base + i)
  fac <- gen_cont_factual(n)
  rows <- list()
  for (K in Ks) {
    d  <- discretize_cont(fac, K)
    Lv <- replicate(K, c("L1", "L2"), simplify = FALSE)
    res <- surv.TMLE(dat = d, Yvar = paste0("Y", 1:K), Cvar = paste0("Censor", 1:K),
                     Avar = "A", Lvar = Lv, L0var = NULL, lookback = 1,
                     Ymod = "SL", Cmod = "SL", Amod = "SL",
                     SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth"),
                     gbound = 0.005, V = NULL, MSM.form = ~A+time, Print = FALSE)
    r <- K + 1
    rows[[length(rows) + 1]] <- data.frame(K = K, ATE = -res$ATE[r], SE = res$SE.ATE[r],
                                           lower = -res$UL.ATE[r], upper = -res$LL.ATE[r])
  }
  saveRDS(do.call(rbind, rows), file = paste0("result/", data.fcn, "_sl4_", i, ".Rds"))
}
