# Example 1 - Single baseline binary covariate 
# Example from the paper

generateData1 <- function(n){
  expit <- plogis;
  
  ## Generate baseline data
  L <- rbinom(n, size = 1, prob = 0.5);
  A <- rbinom(n, size = 1, prob = plogis(-3 + 0.6*L));
  
  ## Generate counterfactual outcome
  #time 1
  py1.1 <- expit(-2 - 1 + 0.25*L);
  py1.0 <- expit(-2 + 0 + 0.25*L);
  Y1.1 <- rbinom(n, 1, py1.1);
  Y1.0 <- rbinom(n, 1, py1.0);
  
  #time 2
  Y2.1 <- Y2.0 <- rep(1, n);
  py2.1 <- expit(-2 - 1 + 0.25*L)[Y1.1 == 0];
  py2.0 <- expit(-2 + 0 + 0.25*L)[Y1.0 == 0];
  Y2.1[Y1.1 == 0] <- rbinom(n = length(py2.1), 1, py2.1);
  Y2.0[Y1.0 == 0] <- rbinom(n = length(py2.0), 1, py2.0);
  
  #time 3
  Y3.1 <- Y3.0 <- rep(1, n);
  py3.1 <- expit(-2 - 1 + 0.25*L)[Y2.1 == 0];
  py3.0 <- expit(-2 + 0 + 0.25*L)[Y2.0 == 0];
  Y3.1[Y2.1 == 0] <- rbinom(n = length(py3.1), 1, py3.1);
  Y3.0[Y2.0 == 0] <- rbinom(n = length(py3.0), 1, py3.0);
  
  #time 4
  Y4.1 <- Y4.0 <- rep(1, n);
  py4.1 <- expit(-2 - 1 + 0.25*L)[Y3.1 == 0];
  py4.0 <- expit(-2 + 0 + 0.25*L)[Y3.0 == 0];
  Y4.1[Y3.1 == 0] <- rbinom(n = length(py4.1), 1, py4.1);
  Y4.0[Y3.0 == 0] <- rbinom(n = length(py4.0), 1, py4.0);
  
  
  ## Generate censoring and observed outcome
  #time 1
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L);
  Censor1 <- rbinom(n, 1, censor.prob1);
  
  #time 2
  Censor2 <- rep(1, n);
  censor.prob2 <- expit(-5 + 0.2*A + 0.2*L)[Censor1 == 0 & Y1 == 0];
  Censor2[Censor1 == 0 & Y1 == 0] <- rbinom(n = length(censor.prob2), 1, censor.prob2);
  Y2 <- Y2.1*A + Y2.0*(1 - A);
  
  #time 3
  Censor3 <- rep(1, n);
  censor.prob3 <- expit(-5 + 0.2*A + 0.2*L)[Censor2 == 0 & Y2 == 0];
  Censor3[Censor2 == 0 & Y2 == 0] <- rbinom(n = length(censor.prob3), 1,
                                            censor.prob3);
  Y3 <- Y3.1*A + Y3.0*(1 - A);
  
  #time 4
  Censor4 <- rep(1, n);
  censor.prob4 <- expit(-5 + 0.2*A + 0.2*L)[Censor3==0 & Y3 == 0];
  Censor4[Censor3 == 0 & Y3 == 0] <- rbinom(n = length(censor.prob4), 1,
                                            censor.prob4);
  Y4 <- Y4.1*A + Y4.0*(1 - A);
  
  
  ## The observed outcome is missing if censored
  Y1[Censor1==1] <- NA; Y2[Censor2==1] <- NA; Y3[Censor3==1] <- NA; Y4[Censor4==1] <- NA;
  
  ## Once an event occurred, future Y = 1
  Y2[Y1 == 1] = 1;  Y3[Y2 == 1] = 1;  Y4[Y3 == 1] = 1; 
  
  ## return counterfactual and observed data
  data.frame(id = 1:n, L, A, Censor1, Censor2, Censor3, Censor4, Y1, Y2, Y3, Y4, Y1.1, Y2.1, Y3.1,
             Y4.1, Y1.0, Y2.0, Y3.0, Y4.0);
}