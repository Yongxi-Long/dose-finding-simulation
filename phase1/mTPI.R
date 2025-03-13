#' @title Modified toxicity probability interval (mTPI) design
#' @description
#' This function implements the mTPI design for dose-finding
#' Modified from R codes from https://biostatistics.mdanderson.org/SoftwareDownload/SingleSoftware/Index/72
#' @param data A data frame from the gen_DLT function containing all possible 
#' outcome conditions of patients at each dose
#' @param start_dose Starting dose level, if not specified would use the lowest 
#' dose
#' @param ndose Number of doses to investigate, should be consistent with the
#' number of columns of the data
#' @param target_DLT Target DLT rate
#' @param csize cohort size for each dose, default is 3
#' @param prior beta prior for the beta-binomial model, default is uniform, i.e., beta(1,1)
#' @param eps1,eps2 
#' (0,target_DLT-eps1) is the under dosing interval
#' (target_DLT-eps1, target_DLT+eps2) is the proper dosing interval
#' (target_DLT+eps2) is the overdosing interval
#' if not specified, eps1 = eps2 = 0.05
#' @param tox_bound the cutoff probability for excessive toxicity, default is 0.95
#' @returns
#' est_MTD: if MTD is one of the supplied doses, then it means a MTD has been found
#' if MTD = 0, it means starting dose is too toxic.
#' trial_data: a data frame containing trial data:
#'    p_id: patient id
#'    dose: assigned dose level
#'    DLT: dose-limiting toxicity yes(1)/no(0)
#' tox_dose: the dose that exceeds the toxicity bound, NA if all supplied doses are safe
#' tox_prob: DLT rate for the lowest toxic dose, NA if all supplied doses are safe
#' @export mTPI_eval
mTPI_eval <- function(data,
                     start_dose = 1,
                     ndose=NULL,
                     target_DLT,
                     csize=3,
                     prior=c(1,1),
                     eps1=0.05,eps2=0.05,
                     tox_bound=0.95
                     )
{
  # get sample size
  N <- nrow(data)
  #' if user does not give numer of dose levels, then take the colunm
  #' number of the DLT data frame
  if(is.null(ndose)) ndose <- ncol(data)
  # indicator if the trial will be stopped for excessive toxicity
  stopped <- 0
  # vector to store to dose exploration trajectory
  dose_levels <- start_dose
  # record any too toxic dose
  toxdose <- ndose+1
  # vector to store observed DLTs for each dose level
  DLT <- rep(0,ndose)
  # vector to store number of patients assigned to each dose
  n <- rep(0,ndose)
  # vector to store individual DLT outcome
  outcome <- numeric()
  # record actual cluster size
  csize_obs <- numeric()
  # vector to store posterior beta parameters for each dose
  a <- b <- rep(NA,ndose)
  mes <- "MTD found"
  # boundary for the indifference interval
  delta1 <- target_DLT - eps1
  delta2 <- target_DLT + eps2
  # initial cohort
  while(stopped==0)
  {
    current_dose <- tail(dose_levels,1)
    # update data
    # if the sample size is not a multiple of csize, the last cohort is < csize
    if(sum(n)+csize > N) {
      cohort <- data[(sum(n)+1):N,current_dose]
    } else
    {
      cohort <- data[(sum(n)+1):(sum(n)+csize),current_dose]
    }
    csize_obs <- c(csize_obs,length(cohort))
    outcome <- c(outcome,cohort)
    DLT[current_dose] <- DLT[current_dose]+sum(cohort)
    n[current_dose] <- n[current_dose]+length(cohort)
    
    # update posterior
    # posterior beta parameters (a,b) for this dose
    a[current_dose] <- prior[1]+DLT[current_dose]
    b[current_dose] <- prior[2]+n[current_dose]-DLT[current_dose]
    # calculate the unit probability mass corresponding to each of the 3 intervals
    upm1 <- pbeta(delta1,a[current_dose],b[current_dose])/delta1
    upm2 <- (pbeta(delta2,a[current_dose],b[current_dose])-pbeta(delta1,a[current_dose],b[current_dose]))/(eps1+eps2)
    upm3 <- pbeta(delta2,a[current_dose],b[current_dose],lower.tail = FALSE)/(1-delta2)
    # assess toxicity probability of the current dose P(d_j|target|data)
    tox_prob_current <- pbeta(target_DLT,a[current_dose],b[current_dose],lower.tail = FALSE)
    # assess toxicity probability of the next dose P(d_{j+1}|target|data)
    if(current_dose < ndose)
    {
      tox_prob_next <- pbeta(target_DLT,prior[1]+DLT[current_dose+1],
                             prior[2]+n[current_dose+1]-DLT[current_dose+1],
                             lower.tail = FALSE)
      # if next dose is too toxic, then cannot escalate (need to make sure that ump1 is the smallest)
      if(tox_prob_next > tox_bound) {upm1 <- 0;toxdose <- current_dose+1} # set ump1 to zero so cannot escalate to the next dose
    }
    
    if(current_dose==1)
    {
      # if the first dose is too toxic, the trial will be terminated for toxicity
      if(tox_prob_current > tox_bound)
      {
        MTD <- 0; stopped <- 1; mes <- "1st dose too toxic"
      } else
      {
        # escalate
        if(upm1==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose+1)}
        # stay
        if(upm2==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose)}
        # de-escalate indicates that 1st dose too toxic
        if(upm3==max(upm1,upm2,upm3)) {MTD <- 0; stopped <- 1; mes <- "1st dose too toxic"}
      }
    } else # for current dose > 1
    {
      if(current_dose==ndose)
      {
        # escalate indicates the final dose is not enough, but we stay as no more dose is available for escalation
        if(upm1==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose)}
        # stay
        if(upm2==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose)}
        # de-escalate 
        if(upm3==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose-1)}
      } else # for dose between 1 and maximum dose
      {
        # escalate
        if(upm1==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose+1)}
        # stay
        if(upm2==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose)}
        # de-escalate 
        if(upm3==max(upm1,upm2,upm3)) {dose_levels <- c(dose_levels,current_dose-1)}
      }
    }
    # if reach the maximum sample size, then do not go to the next dose
    if(sum(n)>=N) {stopped <- 1; dose_levels <- dose_levels[-length(dose_levels)]}
  }
  
  #' compute posterior mean from beta distribution
  #' MTD is selected based on isotonic estimates of the pj calculated using the
  #' pooled adjacent violator algorithm (adapted from MDAnderson website codes)
  if(mes == "MTD found")
  {
    # only consider explored and not too toxic doses
    tdose <- min(max(dose_levels),toxdose-1)
    mindose <- min(dose_levels)
    # estimate for toxicity prob for each dose
    pjs <-rep(-100,tdose)
    # variance of pj
    pjs.var<-rep(0, tdose)
    ##pp.var <- rep(0, D)
    
    for(i in 1:tdose){
      pjs[i] <- (DLT[i]+.005)/(n[i]+.01)
      pjs.var[i] <- betavar(DLT[i]+0.005, n[i]-DLT[i]+0.005) ### here adding 0.005 is to use beta(0.005, 0.005) prior for estimating the MTD, which is different from the dose-finding.
    }
    
    pp<- pava(pjs, wt=1/pjs.var) 
    
    for(i in 2:tdose){
      pp[i] <- pp[i] + i*1E-10 ## by adding an increasingly small number to tox prob at higher doses, it will break the ties and make the lower dose level the MTD if the ties were larger than pT or make the higher dose level the MTD if the ties are smaller than pT
    }
    MTD <-c(mindose:tdose)[sort(abs(pp[mindose:tdose]-target_DLT), index.return=T)$ix[1]]
    # final MTD is selected based on the order-transformed posterior means
  } 
  # prepare output
  trial_data <- data.frame(
    p_id = 1:sum(n),
    dose = rep(dose_levels,csize_obs),
    DLT = outcome
  )
  out <- list(
    est_MTD = MTD,
    trial_data = trial_data,
    message = mes,
    tox_dose=ifelse(toxdose==ndose+1,NA,toxdose),
    tox_prob = ifelse(toxdose==ndose+1,NA,tox_prob_next)
  )
  return(out)
}

#mTPI.sim(data = data,ndose=5,target_DLT = 0.2,eps1 = 0.05,eps2 = 0.05)
#' Other functions need for mTPI design
#' pava is the pool adjacent violator algorithm to perform isotonic transformation for the posterior means later

pava <- function (x, wt = rep(1, length(x))) 
{
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

#' betavar computes variances of beta distributions 
betavar<-function(a,b){
  resp <- a*b/((a+b)^2*(a+b+1))
  return(resp)
}