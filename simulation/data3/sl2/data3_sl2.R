require(SuperLearner)

data.fcn <- "data3"
source(paste0("../../../source/generateData3.R"), local = TRUE)
source(paste0("../../../source/SurvTMLE.R"), local = TRUE)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
if (!dir.exists("result")) {
  dir.create("result")
}

ite <- 200 # 1000
n <- 10000 # 5000

results <- data.frame(true = numeric(ite),
                      ATE = numeric(ite), 
                      SE = numeric(ite),
                      lower = numeric(ite),
                      upper = numeric(ite))

# set.seed(42)
for (i in 1:ite) {
  FullData <- generateData3(n, ntime = 4)
  
  ## True Survival ATE
  true1 <- with(FullData, (1-colMeans(cbind(Y1.1, Y2.1, Y3.1, Y4.1))))
  true0 <- with(FullData, (1-colMeans(cbind(Y1.0, Y2.0, Y3.0, Y4.0))))
  trueATE <- true1 - true0
  
  dat <- FullData
  
  res.sl <- surv.TMLE(dat = dat, 
                      Yvar = c("Y1", "Y2", "Y3", "Y4"),
                      Cvar = c("Censor1", "Censor2", "Censor3", "Censor4"),
                      Avar = "A", 
                      Lvar = list("L1", "L2", "L3", "L4"),
                      L0var = "L",
                      lookback = 2,
                      Ymod = "SL", Cmod = "SL", Amod = "SL",
                      SL.library = c("SL.glm", "SL.glm.interaction"),
                      gbound = 0.005, V = NULL, 
                      MSM.form = ~A+time, Print = FALSE)
  
  result.sl <- c(trueATE[4], -res.sl$ATE[5], res.sl$SE.ATE[5], 
                 -res.sl$UL.ATE[5], -res.sl$LL.ATE[5])
  
  saveRDS(result.sl, file = paste0("result/", data.fcn, "_sl2_", i, ".Rds"))
  
  results[i,] <- result.sl
}

# saveRDS(results, file = paste0("results.Rds"))

