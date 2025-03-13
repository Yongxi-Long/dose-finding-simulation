#' @title Continuous Reassessment Method (CRM) implementation
#' @description
#' This function implements the CRM for dose-finding, dependent on
#' the R package dfcrm
#' @param data A data frame from the gen_DLT function containing all possible 
#' outcome conditions of patients at each dose
#' @param start_dose Starting dose level, if not specified would use the lowest 
#' dose
#' @param ndose Number of doses to investigate, should be consistent with the
#' number of columns of the data
#' @param target_DLT Target DLT rate
#' @param csize cohort size for each dose, default is 1
#' @param skeleton  A vector of initial guesses of toxicity probabilities associated the doses.
#' if user does not supply a skeleton, the default one will be from the getprior
#' function, with middle dose (rounding down of the median if number of doses is even)
#' as the prior guess of MTD, a symmetric indifference interval of length 0.1 around 
#' the target DLT rate
#' @param model A character string specify the working model of CRM. 
#' Supported model type from the dfcrm package. The default
#' model is “empiric”. A one-parameter logistic model is specified by “logistic”.
#' @param incpt The intercept of the working logistic model. The default is 3. 
#' If model=“empiric”,this argument will be ignored.
#' @param conf_level credible interval level, default is 0.95
#' @import dfcrm
#' @returns 
#' est_MTD: if MTD is one of the supplied doses, then it means a MTD has been found
#' if MTD = 0, it means starting dose is too toxic.
#' trial_data: a data frame containing trial data:
#'    p_id: patient id
#'    dose: assigned dose level
#'    DLT: dose-limiting toxicity yes(1)/no(0)
#' @export CRM_eval
CRM_eval <- function(data,
                    start_dose=1,
                    ndose=NULL,
                    target_DLT,
                    csize=1,
                    skeleton=NULL,
                    model="empiric",
                    conf_level = 0.95,
                    incpt=3)
{
  #' if user does not give numer of dose levels, then take the colunm
  #' number of the DLT data frame
  if(is.null(ndose)) ndose <- ncol(data)
  if(is.null(skeleton))
  {
    skeleton <- dfcrm::getprior(halfwidth = 0.05,target = target_DLT,
                             nu = floor(median(1:ndose)),
                             nlevel = ndose,model = model)
  }
  #' get starting data from the 1st cohort of patients
  dose_levels <- start_dose
  count <- 0
  # record individual DLT
  DLT <- data[(count+1):(count+csize),start_dose]
  # record actual cluster size
  csize_obs <- csize
  # starting dose levels for the 1st cohort
  level <- rep(start_dose,csize)
  #' get sample size
  N <- nrow(data)
  stopped <- 0
  while(stopped==0)
  {
    count <- count + csize
    update <- dfcrm::crm(prior = skeleton,target = target_DLT,
                         tox = DLT,level = level,model = model,
                         conf.level = conf_level)
    #' updated estimate of the MTD
    update.mtd <- update$mtd
    #' assign next patient to the updated estimate of MTD and obtain his/her result
    dose_levels <- c(dose_levels,update.mtd)
    # if the sample size is not a multiple of csize, the last cohort is < csize
    if(count + csize >= N)
    {
      cohort <- data[(count+1):N,update.mtd]
      DLT <- c(DLT,cohort)
      
    } else
    {
      cohort <- data[(count+1):(count+csize),update.mtd]
      DLT <- c(DLT,cohort)
    }
    csize_obs <- c(csize_obs,length(cohort))
    level <- rep(dose_levels,csize_obs)
    if(length(DLT) >= N) stopped <- 1
  }
  #prepare output
  trial_data <- data.frame(
    p_id = 1:length(DLT),
    dose = level,
    DLT = DLT
  )
  out <- update[c("prior",
                  "prior.var",
                  "estimate",
                  "post.var",
                  "ptox",
                  "ptoxL",
                  "ptoxU",
                  "conf.level")
  ]
  out$est_MTD <- update.mtd
  out$trial_data <- trial_data
  return(out)
}
