# Scenario 5 - Time-varying treatment effect
# Extends Scenario 2 (generateData2): same binary baseline covariate L and
# binary time-varying covariate L_t, but the conditional treatment effect on the
# discrete-time outcome hazard is allowed to VARY across follow-up intervals.
#
# In Scenarios 1-4 the treatment contribution to every interval's outcome-hazard
# linear predictor is the same fixed offset (-1 for treated, 0 for control), i.e.
# a time-INVARIANT conditional effect. Here the treated-arm offset at interval t
# is -beta[t], so the conditional log-odds treatment effect changes over time.
#
# Arguments
#   n     : sample size
#   ntime : number of discrete time intervals (default 4, matching the paper)
#   beta  : numeric vector of length ntime giving the protective treatment
#           log-odds shift at each interval. Treated-arm hazard uses intercept
#           shift -beta[t]; control-arm uses 0. The default is a WANING protective
#           effect (strong early, attenuating later). Setting beta = rep(1, ntime)
#           exactly recovers generateData2 (constant effect).
#
# The counterfactual outcome columns Y_t.1 / Y_t.0 are returned so the true
# (time-specific) ATE can be recomputed empirically, exactly as in Scenarios 1-4.

generateData5 <- function(n, ntime = 4, beta = c(1.5, 1.0, 0.5, 0.25)){
  expit <- plogis;
  if (length(beta) != ntime) stop("length(beta) must equal ntime");

  ## Generate baseline data
  L  <- rbinom(n, size = 1, prob = 0.5);
  L1 <- rbinom(n, size = 1, prob = 0.3 + 0.4*L);
  A  <- rbinom(n, size = 1, prob = plogis(-3 + 0.6*L + 0.6*L1));

  ## time 1 (treatment effect = -beta[1] for treated arm)
  py1.1 <- expit(-2 - beta[1] + 0.25*L + 0.25*L1);
  py1.0 <- expit(-2 + 0       + 0.25*L + 0.25*L1);
  Y1.1 <- rbinom(n, 1, py1.1);
  Y1.0 <- rbinom(n, 1, py1.0);
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L + 0.2*L1);
  Censor1 <- rbinom(n, 1, censor.prob1);

  for(j in 2:ntime){
    Lj_1   = get(paste0("L", j-1));
    Cj_1   = get(paste0("Censor", j-1));
    Yj_1   = get(paste0("Y", j-1));
    Yj_1.1 = get(paste0("Y", j-1, ".1"));   ## state under A = 1 at previous interval
    Yj_1.0 = get(paste0("Y", j-1, ".0"));   ## state under A = 0 at previous interval
    nameL     = paste0("L", j);
    namepy1   = paste0("py", j, ".1");
    namepy0   = paste0("py", j, ".0");
    nameY1    = paste0("Y", j, ".1");
    nameY0    = paste0("Y", j, ".0");
    nameY     = paste0("Y", j);
    nameCprob = paste0("censor.prob", j);
    nameCensor= paste0("Censor", j);

    assign(nameL, rbinom(n, size = 1, prob = 0.2 + 0.3*L + 0.3*Lj_1));
    Lj = get(nameL);
    ## time-varying treatment effect: treated offset -beta[j], control offset 0
    assign(namepy1, ifelse(Yj_1.1 == 0, expit(-2 - beta[j] + 0.25*L + 0.25*Lj), 1));
    assign(namepy0, ifelse(Yj_1.0 == 0, expit(-2 + 0        + 0.25*L + 0.25*Lj), 1));
    pyj.1 = get(namepy1);
    pyj.0 = get(namepy0);
    assign(nameY1, rbinom(n, 1, pyj.1));
    assign(nameY0, rbinom(n, 1, pyj.0));
    assign(nameCprob, ifelse(Cj_1 == 0 & Yj_1 == 0, expit(-5 + 0.2*A + 0.2*L + 0.2*Lj), 1));
    pCj = get(nameCprob);
    assign(nameCensor, rbinom(n, 1, pCj));
    assign(nameY, get(nameY1)*A + get(nameY0)*(1 - A));
  }

  id = 1:n;
  var.names = c("id", "L", "A", paste0("Censor", 1:ntime), paste0("Y", 1:ntime),
                paste0("L", 1:ntime), paste0("Y", 1:ntime, ".1"), paste0("Y", 1:ntime, ".0"));
  ds = matrix(NA, nrow = n, ncol = length(var.names));
  for(j in 1:length(var.names)){
    ds[,j] = get(var.names[j]);
  }
  ds = data.frame(ds);
  names(ds) = var.names;
  return(ds);
}
