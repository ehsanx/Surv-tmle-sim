generateData2 <- function(n, ntime = 10){
  expit <- plogis;
  
  ## Generate baseline data
  L <- rbinom(n, size = 1, prob = 0.5);
  L1 <- rbinom(n, size = 1, prob = 0.3 + 0.4*L);
  A <- rbinom(n, size = 1, prob = plogis(-3 + 0.6*L + 0.6*L1));
  py1.1 <- expit(-2 - 1 + 0.25*L + 0.25*L1);
  py1.0 <- expit(-2 + 0 + 0.25*L + 0.25*L1);
  Y1.1 <- rbinom(n, 1, py1.1);
  Y1.0 <- rbinom(n, 1, py1.0);
  Y1 <- Y1.1*A + Y1.0*(1 - A);
  censor.prob1 <- expit(-5 + 0.2*A + 0.2*L + 0.2*L1);
  Censor1 <- rbinom(n, 1, censor.prob1);
  
  for(j in 2:ntime){
    Lj_1 = get(paste0("L", j-1));
    Cj_1 = get(paste0("Censor", j-1));
    Yj_1 = get(paste0("Y", j-1));
    Yj_1.1 = get(paste0("Y", j-1, ".1")); ##
    ## add this above one, Y_{j-1}.1 refers to the state at the previous timepoint, for exposure A = 1 
    Yj_1.0 = get(paste0("Y", j-1, ".0")); ##
    ## add this one, Y_{j-1}.0 refers to the state at the previous timepoint, for exposure A = 0
    nameL = paste0("L", j);
    namepy1 = paste0("py", j, ".1");
    namepy0 = paste0("py", j, ".0");
    nameY1 = paste0("Y", j, ".1");
    nameY0 = paste0("Y", j, ".0");
    nameY = paste0("Y", j);
    nameCprob = paste0("censor.prob", j);
    nameCensor = paste0("Censor", j); 
    
    assign(nameL, rbinom(n, size = 1, prob = 0.2 + 0.3*L + 0.3*Lj_1));
    Lj = get(nameL);
    assign(namepy1, ifelse(Yj_1.1 == 0, expit(-2 - 1 + 0.25*L + 0.25*Lj), 1)); 
    ## original version: assign(namepy1, ifelse(Yj_1 == 0, expit(-2 - 1 + 0.25*L + 0.25*Lj), 1));
    ## replaced Yj_1 with Yj_1.1, based on the definition by the paper: 
    ## Y^a_{t} = 1 if Y^a_{t-1} = 0
    ## for exposure A = 1
    assign(namepy0, ifelse(Yj_1.0 == 0, expit(-2 + 0 + 0.25*L + 0.25*Lj), 1));
    ## original version: assign(namepy0, ifelse(Yj_1 == 0, expit(-2 + 0 + 0.25*L + 0.25*Lj), 1));
    ## replaced Yj_1 with Yj_1.0, based on the definition by the paper: 
    ## Y^a_{t} = 1 if Y^a_{t-1} = 0
    ## for exposure A = 0
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