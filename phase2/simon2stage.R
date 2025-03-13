#' @title Simon's two stage design for phase 2
#' @description Use the Simon's two stage design for a certain dose, use the ph2simon function
#' @param p0,p1 assumed response rate of the standard of care and the new intervention
#' @param alpha,beta type I and type II error rate. Power=1-beta
#' @param method method for the Simon's two stage design, default is minimax
#' available options are "Minimax" or "Optimal"
#' @param nmax maximum sample size to solve for the Simon's design, default
#' is 200, maximum possible is 1000. If no solution is found with the current
#' supplied nmax, an error will occur and consider increase the nmax until
#' a solution is found.
#' @export simon2stage_design
#' @returns design parameters for Simon's twostage design
#' r1: number of responses needed to exceeded in first stage, continue to second stage only if responses > r1 at n1
#' n1: number of subjects treated in first stage
#' r2: number of responses needed to exceeded at the end of trial, reject H0 only if responses > r2 at N
#' N: sample size
require(clinfun)
simon2stage_design <- function(p0 = NULL,
                               p1 = NULL,
                               alpha = 0.05, beta = 0.2,
                               method = "Minimax",
                               nmax = 200)
{
  # get design parameters
  design_para <- tryCatch(ph2simon(pu=p0,
                    pa=p1,
                    ep1=alpha,
                    ep2=beta,
                    nmax = nmax),
           error = function(e) {
             e$message
           })
  if(class(design_para)!="ph2simon")
  {
    stop(design_para)
  } else
  {
    #summary(ph2design)
    # use the corresponding minimax or optimal design specified
    # by the method parameter
    r1 <- twostage.admissible(design_para)[method,"r1"]
    n1 <- twostage.admissible(design_para)[method,"n1"]
    r2 <- twostage.admissible(design_para)[method,"r"]
    N <- twostage.admissible(design_para)[method,"n"]
  }
  return(list(r1=r1,n1=n1,r2=r2,N=N))
}
#' @title Simon's two stage design evaluation for phase 2
#' @description Evaluate the performance under Simon's two stage design
#' @param data A data frame containing Phase 2 trial data to be evaluated
#' @param design_para A named list of design parameters, can be either the output of simon2stage_design, 
#' or user-input
#' @export simon2stage_eval
#' @returns Performance metrics under Simon's twostage design
#' reject_H0: A vector of length J indicating whether H0 is rejected for each dose j (j=1,..,J)
#' rsp_rate: A vector of length J giving the response rate for each dose
#' num_pts: A vector of length J giving the number of patients treated at each dose
#' stop_futility: A vector of length J indicating whether the trial is stopped for futility at n1 for each dose
simon2stage_eval <- function(data,
                            design_para)
{
  # check if design parameters are complete
  check <- sapply(c("r1","r2","n1","N"), function(i) is.null(design_para[[i]]))
  if(any(check))
  {
    stop("Design parameters for simon2stage are not complete!")
  }
  # extract number of doses based on supplied data columns
  ndose <- ncol(data)
  # extract design parameters
  r1 <- design_para$r1; r2 <- design_para$r2
  n1 <- design_para$n1; N <- design_para$N
  # evaluate trial under each dose
  # if stopped for futility after 1st stage
  stop_futility <- cumsum(data)[n1,] <= r1
  # if make it into stage 2 and success at the end
  reject_H0 <- cumsum(data)[N,] > r2 & !stop_futility
  # record num of patients per dose
  num_pts <- n1*stop_futility + N*(!stop_futility)
  # record num of responses per dose
  num_rsp <- as.vector(unlist(cumsum(data)[n1,]*stop_futility + cumsum(data)[N,]*!stop_futility))
  # response rate per dose
  rsp_rate <- num_rsp/num_pts 
  return(list(
    reject_H0 = reject_H0,
    rsp_rate = rsp_rate,
    num_pts = num_pts,
    stop_futility = stop_futility
  ))
}
