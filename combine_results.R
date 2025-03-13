combine_results <- function(res_ph1=NULL,res_ph2=NULL,res_ph3=NULL)
{
  if(is.null(res_ph1) | is.null(res_ph2))
  {
    stop("Please at least supply the simulation results from phase 1 and phase 2")
  }
  # calculate phase 1 & 2 combined power
  # match results based on dose
  res_ph12 <- dplyr::inner_join(res_ph1$`summary table`,
                                res_ph2$`summary table`,
                                by = "dose",
                                suffix = c("_ph1","_ph2"))
  res_ph12 <- res_ph12 |>
    dplyr::mutate(avg_num_pts_ph12 = avg_num_pts_ph1 + avg_num_pts_ph2*prop_seleted_MTD)
  # combined power 
  power_ph12 <- sum(res_ph12$prop_seleted_MTD * res_ph12$avg_reject_H0)
  if(is.null(res_ph3))
  {
    # prepare phase 1 2 output
    out <- list()
    out$`Doses levels` <- res_ph12$dose
    out$`Design for each phase` <- c("phase 1"= res_ph1$design,
                                     "phase 2"=res_ph2$design)
    out$`Sample size for each phase` <- c("phase 1"= res_ph1$`sample size`,
                                          "phase 2"=res_ph2$`sample size`)
    out$`Combined power` <- c("phase 1&2"=power_ph12)
    out$`Summary table` <- res_ph12 |>
      dplyr::select(dose,avg_num_pts_ph12)
  } else
  {
    # calculate phase 1 & 2 combined power
    # match results based on dose
    res_ph123 <- dplyr::inner_join(res_ph12,res_ph3$`summary table`,
                                   by = "dose",
                                   suffix = c("_ph12","_ph3"))
    res_ph123 <- res_ph123 |>
      dplyr::mutate(avg_num_pts_ph123 = avg_num_pts_ph12 
                    + avg_num_pts*prop_seleted_MTD*avg_reject_H0_ph12)
    # combined power 
    power_ph123 <- sum(res_ph123$prop_seleted_MTD * res_ph123$avg_reject_H0_ph12
                       *res_ph123$avg_reject_H0_ph3)
    
    # prepare phase 1 2 output
    out <- list()
    out$`Doses levels` <- res_ph12$dose
    out$`Design for each phase` <- c("phase 1"= res_ph1$design,
                                     "phase 2"=res_ph2$design,
                                     "phase 3" = res_ph3$design)
    out$`Sample size for each phase` <- c("phase 1"= res_ph1$`sample size`,
                                          "phase 2"=res_ph2$`sample size`,
                                          "phase 3"=res_ph3$`sample size`)
    out$`Combined power` <- c("phase 1&2"=power_ph12, "phase 1&2&3" = power_ph123)
    out$`Summary table` <- res_ph123 |>
      dplyr::select(dose,avg_num_pts_ph12,avg_num_pts_ph123)
  }
  return(out)
}
