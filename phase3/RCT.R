#' @title randomized controlled trial design for phase 3
#' @description Use the two-arm randomized controlled design. By default the outcome
#' is time-to-event 
#' @param alpha,beta type I and type II error rate. Power=1-beta 
#' @param p_trt probability of being randomized to the treatment arm
#' @param hazard0,hazard1 event hazard for the control and the treatment group respectively
#' @param t_FU time duration for study follow-up
#' @param t_accrual time duration for subject accrual
#' @param hazard_dropout drop out hazard (equal) for both groups, default is zero
#' @export RCT_design
#' @returns design parameters for RCT
#' N: total sample size
#' alpha: significance level of the test
require(gsDesign)
RCT_design <- function(alpha = 0.05, beta = 0.2,
                       p_trt = 0.5,
                       hazard0, hazard1,
                       t_FU, t_accrual,
                       hazard_dropout = 0)
{
  library(gsDesign)
  design_para <- nSurvival(alpha = alpha,
                           beta = beta,
                           lambda1 = hazard0,
                           lambda2 = hazard1,
                           Ts = t_FU,
                           Tr = t_accrual,
                           eta = hazard_dropout,
                           ratio = (1-p_trt)/p_trt)
  # get the sample size
  N <- round(design_para$n)
  return(list(N = N, alpha = alpha))
}
  
  
#' @title RCT design evaluation for phase 3
#' @description
#' Evaluate the performance under a randomized controlled trial design
#'  description
#' @param data A list of J data frames for the J doses containing the Phase 3 survival data to be evaluated,
#' must have the following four columns:
#' id: patient id
#' trt: treatment assignment
#' eventtime: time to event
#' status: event indicator
#' @param design_para A named list of design parameters, can be either the output of
#' RCT_design or user-input
#' @export RCT_eval
#' @returns Performance metrics under RCT design
#' reject_H0: A vector of length J indicating whether H0 is rejected for each dose j (j=1,..,J)
#' hazard_ratio: A vector of length J giving the hazard ratio of treatment vs. control
#'  for each dose
#' median_surv: A vector of length J giving the median survival for each dose
#' rmean_surv: A vector of length J giving the restricted mean survival until the maximum followup time (maxt)
#' for each dose
#' num_pts: A vector of length J giving the number of patients treated at each dose
require(survival)
RCT_eval <- function(data,
                   design_para)
{
  library(survival)
  # significance level for the test
  alpha <- design_para$alpha
  # number of patients
  res <- sapply(data, function(data_thisdose)
    {
    num_pts <- nrow(data_thisdose)
    fit_km <- survfit(Surv(eventtime, status) ~ trt, data = data_thisdose)
    mod_cox <- coxph(Surv(eventtime,status) ~ trt, data = data_thisdose,ties = "breslow")
    HR <- exp(mod_cox$coefficients["trt"]) |>
      unname()
    p_value <- summary(mod_cox)$logtest["pvalue"] |>
      unname()
    reject <- p_value < alpha
    rmean_surv <- summary(fit_km)$table[2,"rmean"]
    median_surv <- summary(fit_km)$table[2,"median"]
    return(c(reject_H0=reject, 
             hazard_ratio=HR, 
             median_surv = median_surv,
             rmean_surv = rmean_surv,
             num_pts = num_pts))
  })
  return(list(
    reject_H0 = res["reject_H0",],
    hazard_ratio=res["hazard_ratio",], 
    median_surv = res["median_surv",],
    rmean_surv = res["rmean_surv",],
    num_pts = res["num_pts",]
  ))
}
