#' @title Single-arm design for phase 2
#' @description Use the single arm design for a certain dose, use the ph2single function
#' @param p0,p1 response rate of the standard of care and the new intervention
#' @param alpha,beta type I and type II error rate. Power=1-beta
#' @export single_arm_design
#' @returns design parameters for single arm design
#' r: number of responses needed to exceeded at the end of trial
#' reject H0 only if responses > r at N
#' N: sample size
require(clinfun)
single_arm_design <- function(p0=NULL,
                              p1=NULL,
                              alpha=0.05,
                              beta=0.2)
{
  design_para <- tryCatch(ph2single(pu=p0,
                                    pa=p1,
                                    ep1=alpha,
                                    ep2=beta,
                                    nsoln = 1),
                          error = function(e)
                            {
                            e$message
                          })
  if(class(design_para)!="data.frame")
  {
    stop(design_para)
  } else
  {
    r <- design_para$r
    N <- design_para$n
  }
return(list(N=N,r=r)) 
}


#' @title Single-arm design evaluation for phase 2
#' @description Evaluate the performance under Sinle-arm design
#' @param data A data frame containing Phase 2 trial data to be evaluated
#' @param design_para A named list of design parameters, can be either the output of single_arm_design, 
#' or user-input
#' @export single_arm_eval
#' @returns Performance metrics under single-arm design
#' reject_H0: A vector of length J indicating whether H0 is rejected for each dose j (j=1,..,J)
#' rsp_rate: A vector of length J giving the response rate for each dose
#' num_pts: A vector of length J giving the number of patients treated at each dose
single_arm_eval <- function(data,design_para)
{
  # check if design parameters are complete
  check <- sapply(c("r","N"), function(i) is.null(design_para[[i]]))
  if(any(check))
  {
    stop("Design parameters for single_arm are not complete!")
  }
  
  # extract number of doses based on supplied data columns
  ndose <- ncol(data)
  # extract design parameters
  r <- design_para$r; N <- design_para$N
  # evaluate trial under each dose
  # if success at the end
  reject_H0 <- cumsum(data)[N,] > r 
  # record num of patients per dose
  num_pts <- rep(N,ndose)
  # record num of responses per dose
  num_rsp <-  as.vector(unlist(cumsum(data)[N,]))
  # response rate per dose
  rsp_rate <- num_rsp/num_pts 
  return(list(
    reject_H0 = reject_H0,
    rsp_rate = rsp_rate,
    num_pts = num_pts
  ))
}
  

