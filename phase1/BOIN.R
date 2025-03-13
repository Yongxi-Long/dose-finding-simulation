#' @title Bayesian Optimal Interval (BOIN) design
#' @description
#' This function implements the BOIN design for dose-finding
#' dependent on the BOIN package
#' @param data A data frame from the gen_DLT function containing all possible 
#' outcome conditions of patients at each dose
#' @param start_dose Starting dose level, if not specified would use the lowest 
#' dose
#' @param ndose number of doses in investigation
#' @param target_DLT Target DLT rate
#' @param csize cohort size for each dose, default is 3
#' #' BOIN design with an uninformative prior is equivalent to dose escalation with LRT
#' need to specify two limits:
#' @param p_L highest DLT probability that is deemed to be under-dosing 
#' such that dose escalation is required, default is 0.6*target_DLT
#' @param p_U lowest DLT probability that is deemed to be overdosing 
#' such that dose de-escalation is required, default is 1.4*target_DLT
#' @param prior beta prior for the beta-binomial model, default is uniform, i.e., beta(1,1)
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
#' @export BOIN_eval

BOIN_eval <- function(data, 
                     start_dose = 1, ndose=NULL,
                     target_DLT,
                     csize=3, p_L = NULL, p_U = NULL,
                     prior = c(1,1),
                     tox_bound = 0.95)
{
  # get sample size
  N <- nrow(data)
  #' if user does not give numer of dose levels, then take the colunm
  #' number of the DLT data frame
  if(is.null(ndose)) ndose <- ncol(data)
  # vector to store observed DLTs for each dose level
  DLT <- rep(0,ndose)
  # vector to store individual DLT outcome
  outcome <- numeric()
  # record actual cluster size
  csize_obs <- numeric()
  # vector to store number of patients assigned to each dose
  n <- rep(0,ndose)
  # vector to store posterior beta parameters for each dose
  a <- b <- rep(NA,ndose)
  # vector to store to dose exploration trajectory
  dose_levels <- start_dose
  # record any too toxic dose
  toxdose <- ndose+1
  
  # get the boundaries for escalation
  if(is.null(p_L)) {p_L <- 0.6*target_DLT}
  if(is.null(p_U)) {p_U <- 1.4*target_DLT}
  bound <- BOIN::get.boundary(target = target_DLT,
                              ncohort = 100,
                              cohortsize = csize,
                              p.saf = p_L,
                              p.tox = p_U,
                              cutoff.eli = tox_bound)
  lambda_e <- bound$lambda_e
  lambda_d <- bound$lambda_d
  
  stopped <- 0
  mes <- "MTD found"
  # dose-escalation
  while(stopped==0)
  {
    # get current dose, trim toxic doses
    current_dose <- tail(dose_levels,1)
    if(current_dose >= toxdose)
    {
      current_dose <- toxdose - 1
      dose_levels[length(dose_levels)] <- current_dose
    }
    # update data
    # if the sample size is not a multiple of csize, the last cohort is < csize
    if(sum(n)+csize > N) {
      cohort <- data[(sum(n)+1):N,current_dose]
    } else
    {
      cohort <- data[(sum(n)+1):(sum(n)+csize),current_dose]
    }
    outcome <- c(outcome,cohort)
    DLT[current_dose] <- DLT[current_dose]+sum(cohort)
    n[current_dose] <- n[current_dose]+length(cohort)
    csize_obs <- c(csize_obs,length(cohort))
    # observed DLT rate for current dose
    DLT.obs <- DLT[current_dose]/n[current_dose]
    
    # see if current dose is too toxic and should eliminate it and all higher doses
    # posterior beta parameters (a,b) for this dose
    a[current_dose] <- prior[1]+DLT[current_dose]
    b[current_dose] <- prior[2]+n[current_dose]-DLT[current_dose]
    tox_prob_current <- pbeta(target_DLT,a[current_dose],b[current_dose],lower.tail = FALSE)
    if(tox_prob_current > tox_bound)
    {
      toxdose <- current_dose # eliminate too toxic doses
    } 
    
    # compare to two boundaries to decide (de)-escalation
    if(current_dose==1)
    {
      # if the first dose is too toxic, the trial will be terminated for toxicity
      if(toxdose==1)
      {
        MTD <- 0; stopped <- 1; mes <- "1st dose too toxic"
      } else
      {
        # escalate
        if(DLT.obs <= lambda_e) 
        {
          dose_levels <- c(dose_levels,current_dose+1)
        } # de-escalate
        else if(DLT.obs >= lambda_d) 
        {
          MTD <- 0 ; stopped <- 1; mes <- "1st dose too toxic"
        } else # stay at current dose
        {
          dose_levels <- c(dose_levels,current_dose)
        }
      }
    } else # for current dose > 1
    {
      if(current_dose==ndose)
      {
        # escalate indicates the final dose is not enough, we stay st the current dose
        if(DLT.obs <= lambda_e) 
        {
          dose_levels <- c(dose_levels,current_dose)
        } # de-escalate
        else if(DLT.obs >= lambda_d) 
        {
          dose_levels <- c(dose_levels,current_dose - 1)
        } else # stay at current dose
        {
          dose_levels <- c(dose_levels,current_dose)
        }
        
      } else # for dose between 1 and maximum dose
      {
        # escalate
        if(DLT.obs <= lambda_e) 
        {
          dose_levels <- c(dose_levels,current_dose+1)
        } # de-escalate
        else if(DLT.obs >= lambda_d) 
        {
          dose_levels <- c(dose_levels,current_dose - 1)
        } else # stay at current dose
        {
          dose_levels <- c(dose_levels,current_dose)
        }
      }
    } # end else current dose > 1
    # if reach the maximum sample size, then do not go to the next dose
    if(sum(n)>=N) {stopped <- 1; dose_levels <- dose_levels[-length(dose_levels)]}
  }
  # if not stopped for toxicity
  if(mes == "MTD found")
  {
    # estimate MTD via isotonic regression
    mtd.obj <- BOIN::select.mtd(target = target_DLT,
                                npts = n,
                                ntox = DLT,
                                cutoff.eli = tox_bound,
                                boundMTD = TRUE,
                                p.tox = p_U)
    MTD <- mtd.obj$MTD
  }
  # prepare output
  # trial data
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
    tox_prob = ifelse(toxdose==ndose+1,NA,tox_prob_current)
  )
  return(out)
}
