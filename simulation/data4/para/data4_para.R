require(SuperLearner)

data.fcn <- "data4"
source(paste0("../../../source/generateData4.R"), local = TRUE)
source(paste0("../../../source/SurvTMLE.R"), local = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!dir.exists("result")) {
  dir.create("result")
}

# split <- 0
ite <- 1 # 1000
n <- 10000 # 5000

results <- data.frame(true = numeric(ite),
                      ATE = numeric(ite), 
                      SE = numeric(ite),
                      lower = numeric(ite),
                      upper = numeric(ite))
# set.seed(42)
for (i in 1:ite) {
  FullData <- generateData4(n, ntime = 4)
  
  ## True Survival ATE
  true1 <- with(FullData, (1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1))))
  true0 <- with(FullData, (1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0))))
  trueATE <- true1 - true0
  
  dat <- FullData
  
  res.para <- surv.TMLE(dat = dat, 
                      Yvar = c("Y1", "Y2", "Y3", "Y4"),
                      Cvar = c("Censor1", "Censor2", "Censor3", "Censor4"),
                      Avar = "A", 
                      Lvar = list(c("L1", "Z1"), 
                                  c("L2", "Z2"), 
                                  c("L3", "Z3"), 
                                  c("L4", "Z4")),
                      L0var = "L",
                      lookback = 2,
                      Ymod = "parametric", Cmod = "parametric", Amod = "parametric",
                      gbound = 0.005, V = 5, 
                      MSM.form = ~A+time, Print = FALSE)
  
  result.para <- c(trueATE[4], -res.para$ATE[5], res.para$SE.ATE[5], 
                   -res.para$UL.ATE[5], -res.para$LL.ATE[5])
  
  saveRDS(result.para, file = paste0("result/", data.fcn, "_para_", i, ".Rds"))
  
  results[i,] <- result.para
}

# saveRDS(results, file = paste0("results.Rds"))

