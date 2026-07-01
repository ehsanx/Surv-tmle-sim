# Continuous-time data-generating mechanism + discretizer.
# Addresses AE6 / Reviewer 2: the real examples originate from CONTINUOUS survival
# time, so here survival and censoring times are generated from continuous-time
# (Weibull proportional-hazards) distributions and then DISCRETIZED at varying
# resolutions (K intervals over [0, tau]). The discrete-time TMLE estimators are then
# run on each discretized version; comparing their bias/coverage across K (against the
# true CONTINUOUS-time effect at tau) shows how discretization resolution matters.
#
# Confounders L1, L2 are time-fixed and continuous; their effect on the event hazard
# is linear (log-HR), so the parametric discrete-time model is well specified in the
# covariates and the ONLY structural error is discretization granularity.

cont_params <- list(
  tau     = 1.0,                 # fixed follow-up horizon (the estimand is the ATE at tau)
  alpha0  = -0.5, alphaL = c(0.5, 0.5),   # exposure logit:  P(A=1) = expit(a0 + aL%*%(L1,L2))
  betaA   = -0.7, betaL = c(0.5, 0.5),    # event log-HR:   A (protective), L1, L2
  lambdaT = 0.9,  rhoT = 1.4,             # event Weibull (scale, shape)
  gammaA  = 0.2,  gammaL = 0.2,           # censoring log-HR: A, L1
  lambdaC = 0.5,  rhoC = 1.0              # censoring Weibull (exponential when rhoC=1)
)

## Factual continuous-time data: covariates, exposure, latent event time T, censoring C.
## Weibull-PH inversion (Bender, Augustin & Schumacher 2005):
##   T = ( -log(U) / (lambda * exp(linpred)) )^(1/rho).
gen_cont_factual <- function(n, p = cont_params) {
  L1 <- rnorm(n); L2 <- rnorm(n)
  A  <- rbinom(n, 1, plogis(p$alpha0 + p$alphaL[1]*L1 + p$alphaL[2]*L2))
  U  <- runif(n); Vv <- runif(n)
  lpT <- p$betaA*A + p$betaL[1]*L1 + p$betaL[2]*L2
  Tt  <- (-log(U)  / (p$lambdaT*exp(lpT)))^(1/p$rhoT)
  lpC <- p$gammaA*A + p$gammaL*L1
  Ct  <- (-log(Vv) / (p$lambdaC*exp(lpC)))^(1/p$rhoC)
  data.frame(L1, L2, A, Tt, Ct)
}

## Discretize factual continuous data onto K equal intervals over [0, tau].
## Encoding matches the real-data analyses (pbc.R / rhc.R):
##   Y_t      = 1 if the event occurred and obs_time <= t*w           (cumulative event)
##   Censor_t = 1 if no observed event and obs_time <  t*w            (censored before t*w)
## obs_time = min(T, C, tau); event = 1 iff T <= C and T <= tau.
## Subjects with T > tau and C > tau are administratively censored at tau (Y=0, Censor=0).
discretize_cont <- function(df, K, tau = cont_params$tau) {
  w     <- tau / K
  obs   <- pmin(df$Tt, df$Ct, tau)
  event <- as.integer(df$Tt <= df$Ct & df$Tt <= tau)
  out   <- data.frame(id = seq_len(nrow(df)), L1 = df$L1, L2 = df$L2, A = df$A)
  for (t in 1:K) {
    th <- t*w
    out[[paste0("Y", t)]]      <- as.integer(event == 1 & obs <= th)
    out[[paste0("Censor", t)]] <- as.integer(event == 0 & obs <  th)
  }
  out
}

## High-precision TRUE continuous-time ATE on the survival scale at tau,
## S^1(tau) - S^0(tau), via common-random-number counterfactual event times.
true_ate_cont <- function(n = 5e6, tau = cont_params$tau, p = cont_params, seed = 7777) {
  set.seed(seed)
  L1 <- rnorm(n); L2 <- rnorm(n); U <- runif(n)
  base <- p$betaL[1]*L1 + p$betaL[2]*L2
  T1 <- (-log(U) / (p$lambdaT*exp(p$betaA*1 + base)))^(1/p$rhoT)
  T0 <- (-log(U) / (p$lambdaT*exp(p$betaA*0 + base)))^(1/p$rhoT)
  S1 <- mean(T1 > tau); S0 <- mean(T0 > tau)
  c(S1 = S1, S0 = S0, ATE = S1 - S0)
}

## Convenience: descriptive rates at tau for a factual draw (for calibration).
cont_rates <- function(df, tau = cont_params$tau) {
  event <- df$Tt <= df$Ct & df$Tt <= tau
  cens  <- df$Ct <  df$Tt & df$Ct <= tau
  admin <- df$Tt >  tau & df$Ct > tau
  c(event = mean(event), censored = mean(cens), admin = mean(admin),
    treated = mean(df$A))
}
