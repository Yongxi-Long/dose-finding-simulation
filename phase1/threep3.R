#' @title 3+3 Design implementation
#' @description
#' This function implements the rule-based 3+3 design
#' @param data A data frame from the gen_DLT function or from user input that
#' contains all possible outcome conditions of patients at each dose
#' @param start_dose Starting dose level, if not specified would use the lowest 
#' dose
#' @param ndose number of doses in investigation
#' @returns 
#' est_MTD: if MTD is one of the supplied doses, then it means a MTD has been found
#' if MTD = 0, it means starting dose is too toxic.
#' trial_data: a data frame containing trial data:
#'    p_id: patient id
#'    dose: assigned dose level
#'    DLT: dose-limiting toxicity yes(1)/no(0)
#' @export threep3_eval
threep3_eval <- function(data,start_dose=1,ndose=NULL)
{
  #' if user does not give number of dose levels, then take the colunm
  #' number of the DLT data frame
  if(is.null(ndose)) ndose <- ncol(data)
  #' partition the data into cohorts of three patients
  cohort_labels <- (0:(nrow(data) - 1) %/%  3)+1
  num_cohorts <- length(unique(cohort_labels))
  dose_levels <- start_dose
  outcome <- numeric()
  mes <- "MTD found"
  for(i in unique(cohort_labels))
  {
    #number of total doses
    current_dose <- tail(dose_levels,1)
    previous_dose <- ifelse(length(dose_levels)>1,head(tail(dose_levels,2),1),NA)
    #' see if enough patients have been evaluated at this dose level
    #' if it is the first dose
    if(is.na(previous_dose))
    {
      # if current dose is the first dose, then there can only be three patients
      cohort <- data[which(cohort_labels==i),current_dose]
      outcome <- c(outcome,cohort)
      if(sum(cohort)>1)
      {
        # no MTD because already too toxic at the first dose
        mes <- "1st dose too toxic"
        MTD <- 0
        break
      } else if (sum(cohort)==0)
      {
        dose_levels <- c(dose_levels,current_dose+1)
      } else if (sum(cohort)==1)
      {
        # add three more patients to this level
        dose_levels <- c(dose_levels,current_dose)
      }
    } else
    {
      # dose levels > 1
      # if six patients have been evaluated at this dose
      six_pts_thisdose <- current_dose == previous_dose
      if(six_pts_thisdose)
      {
        cohort <- data[c(which(cohort_labels==i-1),
                         which(cohort_labels==i)),current_dose]
        outcome <- c(outcome,tail(cohort,3))
        if(sum(cohort)>1)
        {
          # if six patients at the 1st level and more than 1 DLT: first level too toxic
          if(current_dose==1)
          {
            MTD <- 0
            mes <- "1st dose too toxic"
          } else
          {
            MTD <- current_dose - 1
          }
          break
        } else
        {
          # escalate if it is not the highest dose and there are still patients,
          # otherwise current dose is recommended as the MTD
          if(current_dose==ndose | i==num_cohorts)
          {
            MTD <- ndose
            break
          } else
          {
            dose_levels <- c(dose_levels,current_dose+1)
          }
        }
      } else  # if only three patients at this dose
      {
        cohort <- data[which(cohort_labels==i),current_dose]
        outcome <- c(outcome,cohort)
        if(sum(cohort)>1)
        {
          MTD <- current_dose - 1
          break
        } else if (sum(cohort)==0)
        {
          # escalate if it is not the highest dose and there are still patients,
          # otherwise current dose is recommended as the MTD
          if(current_dose==ndose | i==num_cohorts)
          {
            MTD <- ndose
            break
          } else
          {
            dose_levels <- c(dose_levels,current_dose+1)
          }
        } else if (sum(cohort)==1)
        {
          # add three more patients to this level if there are more patients
          if(i == num_cohorts)
          {
            MTD <- ndose
            break
          } else
          {
            dose_levels <- c(dose_levels,current_dose)
          }
        }
      }
    } # end of else statement for the first dose
  } # end of for loop
  # prepare output
  # trial data
  trial_data <- data.frame(
    p_id = 1:(3*length(dose_levels)),
    dose = rep(dose_levels,each=3),
    DLT = outcome
  )
  out <- list(
    est_MTD = MTD,
    trial_data = trial_data,
    message = mes
  )
  return(out)
}
