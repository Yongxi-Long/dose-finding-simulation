#' @title Simulate phase 1 trial
#' @description
#' The function does MC simulation for phase 1 trials given input design parameters
#' @name simulate_ph1
#' @param nsim number of MC simulations to run
#' @param N sample size for phase 1 trial
#' @param data_list a list of data frames that user can input instead of relying on
#' provided internal data generation. The length of the list is taken as the number of 
#' MC iterations
#' @param method dose-finding method
#' @param tox_profile toxicity profile, should be a function that takes a dose
#' value as input and outputs the DLT rate
#' @param dose_levels a vector of dose levels to be investigated
#' @param DLT_rates users can directly specify a vector of DLT rates for the corresponding
#' dose levels instead of a full function (tox_profile) and a vector of dose levels (dose_levels)
#' @param seed random seed for simulation
#' @param ... other parameters specific for the chosen dose-finding method
#' @returns summary statistics across all simulations
#' tab_MTD: 
require(dplyr)
simulate_ph1 <- function(nsim = NULL,
                         N = NULL,
                         data_list = NULL,
                         design,
                         tox_profile=NULL,
                         dose_levels=NULL,
                         DLT_rates=NULL,
                         seed=NULL,
                         ...)
{
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr is required but not installed.")
  }
  # get all input arguments
  call <- match.call(expand.dots = TRUE)
  if (!is.null(seed)) set.seed(seed)
  # check if user has supplied simulated subject data, if not, call gen_DLT
  if(is.null(data_list))
  {
    data_list <- gen_DLT(nsim = nsim,
                         drug_profile=tox_profile,
                         dose_levels = dose_levels,
                         DLT_rates = DLT_rates,
                         N = N)
  }
  # generate the phase 1 data 
  source("data_generating_functions.R") # based on phase 
  
  # source the corresponding design method
  source(paste0("phase1/",design,".R"))
  default_args <- formals(paste0(design,"_eval"))
  call_args <- as.list(call[-1]) # first one is the function name, not needed
  # change default values to user-input ones, if there is any
  final_args <- default_args
  final_args[names(call_args)] <- call_args
  # keep arguments relevant for threep3.sim only
  final_args <- final_args[names(default_args)]
  res <- sapply(data_list,function(data)
    {
    final_args$data <- data
    do.call(paste0(design,"_eval"),final_args)
  },simplify = FALSE)
  
  # return some design parameters
  out <- list()
  out$phase <- "phase1"
  out$design <- design
  out$`sample size` <- nrow(data_list[[1]])
  out$`number of iterations` <- length(data_list)
  # summary table
  MTD_seq <- sapply(res, function(i) i$est_MTD)
  MTD_tab <- data.frame(table(factor(MTD_seq, levels = 0:ncol(data_list[[1]]))))
  colnames(MTD_tab) <- c("dose","prop_seleted_MTD")
  MTD_tab$prop_seleted_MTD <- MTD_tab$prop_seleted_MTD/length(data_list)
  
  tmp=sapply(res, function(i) {
    dat <- i$trial_data |>
      dplyr::group_by(dose) |>
      dplyr::summarise(num_pts=dplyr::n(),DLT_rate=mean(DLT))
  },simplify = FALSE)
  
  patient_tab <- dplyr::bind_rows(tmp) |>
    dplyr::group_by(dose) |>
    dplyr::summarise(avg_DLT_rate=sum(DLT_rate)/length(data_list),
              avg_num_pts = sum(num_pts)/length(data_list))
  patient_tab$dose <- factor(patient_tab$dose,levels = 1:ncol(data_list[[1]]))
  summary_tab <- dplyr::left_join(MTD_tab,patient_tab,by="dose")
  out$`summary table` <- summary_tab
  return(out)
}


#' @title Simulate phase 2 trial
#' @description
#' The function does MC simulation for phase 2 trials given input design parameters
#' @name simulate_ph2
#' @param nsim number of MC simulations to run
#' @param design  design for phase 2, default is Simon's two-stage
#' @param alpha type I error rate, default is 0.05
#' @param beta type II error rate, default is 0.2 (1-beta) is the target power
#' @param rsp_profile dose-response profile, should be a function that takes a dose
#' value as input and outputs the objective response rate (ORR)
#' @param dose_levels a vector of dose levels to be investigated
#' @param rsp_rates users can directly specify a vector of ORR for the corresponding
#' dose levels instead of a full function (rsp_profile) and a vector of dose levels (dose_levels)
#' @param ph1_res results from phase 1 simulation. Optional, if supplied, will yield additional
#' summary results for phase 1&2 combined.
#' @param seed random seed for simulation
#' @param ... other parameters specific for the chosen phase 2 design
#' @returns summary statistics across all simulations
#' tab_MTD: 
simulate_ph2 <- function(nsim,
                         data_list = NULL,
                         design="simon2stage",
                         design_para = NULL,
                         rsp_profile=NULL,
                         dose_levels=NULL,
                         rsp_rates=NULL,
                         ph1_res = NULL,
                         seed=NULL,
                         ...)
{
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr is required but not installed.")
  }
  # get all input arguments
  call <- match.call(expand.dots = TRUE)
  if (!is.null(seed)) set.seed(seed)
  
  # source the corresponding design method
  source(paste0("phase2/",design,".R"))
  # if user does not supply design parameters for the design method
  # then internally calls the design function
  if(is.null(design_para))
  {
    # get sample size/other design parameters from specified design
    default_args <- formals(paste0(design,"_design"))
    call_args <- as.list(call[-1]) # first one is the function name, not needed
    # change default values to user-input ones, if there is any
    final_args <- default_args
    final_args[names(call_args)] <- call_args
    # keep arguments relevant only
    final_args <- final_args[names(default_args)]
    design_para <- do.call(paste0(design,"_design"),final_args)
  } else
  {
    # check if the design parameters are complete
    names_design_para <- extract_return_names(paste0(design,"_design"))
    if(!identical(sort(names_design_para),sort(names(design_para))))
    {
      stop(paste0("Design parameters for ",design," are not complete!"))
    }
  }
  # get sample size
  # sample size
  N <- design_para$N
  
  # if data is supplied externally, check if the sample size N matches
  # the number of observations in the data frame
  if(!is.null(data_list))
  {
    nsim <- length(data_list)
    if(N > nrow(data_list[[1]]))
    {
      stop("The supplied data set has fewer observations than the design sample size, please input a larger data set!")
    }
    if(N < nrow(data_list[[1]]))
    {
      data_list <- sapply(data_list, function(data)
        {
        data[1:N,]
      },simplify = FALSE)
      warning(paste0("The supplied data set has more observations than the design sample size, only the first ",N, " observations will be used!"))
    }
  } else {
    # generate list of data frames based on calculated sample size
    data_list <- gen_rsp(nsim = nsim,
                         drug_profile=rsp_profile,
                         dose_levels = dose_levels,
                         rsp_rates = rsp_rates,
                         N = N)
  }
  
  res <- sapply(data_list,function(data)
  {
    do.call(paste0(design,"_eval"),list(data = data, design_para = design_para))
  },simplify = FALSE)
  # prepare output
  # return some design features
  out <- list()
  out$phase <- "phase2"
  out$design <- design
  out$`sample size` <- design_para$N
  out$`number of iterations` <- nsim
  # summary table
  # see how many categories are returned by the design function
  res_names <- names(res[[1]])
  summary_tab <- data.frame(matrix(NA,
                                   ncol = length(res_names)+1,
                                   nrow = ncol(data_list[[1]])))
  colnames(summary_tab) <- c("dose",paste0("avg_",res_names))
  summary_tab$dose <- factor(1:ncol(data_list[[1]]))
  summary_tab[,paste0("avg_",res_names)]  = sapply(res_names,function(x) rowMeans(sapply(res,function(i) i[[x]])))
  out$`summary table` <- summary_tab
  return(out)
}


#' @title Simulate phase 3 trial
#' @description
#' The function does MC simulation for phase 3 trials given input design parameters
#' @name simulate_ph3
#' @param nsim number of MC simulations to run
#' @param N sample size for phase 2 trial
#' @param design  design for phase 3, default is two-arm RCT
#' @param alpha type I error rate, default is 0.05
#' @param median_surv0: median survival for the control group. Then internally
#' the function will convert it to hazard assuming the survival probability follows
#' an exponential distribution
#' @param median_surv1: a vector of median survivals for the treatment group under each
#' dose. 
#' @param maxt maximum follow-up time (will apply administrative censoring for survival
#' times larger than this)
#' @param use_simsurv: default is FALSE. If users wish to generate more complicated
#' survival data, they can set this to be TRUE and supply parameters for the simsurv()
#' function. 
#' Note: if users wish to simulate survival data for multiple doses using simsurv() function,
#' then all the corresponding parameters should be lists, with the j-th element
#' corresponding to parameters of the j-th dose
#' @param ph1_res results from phase 1 simulation. Optional, if supplied, will yield additional
#' summary results for phase 1&2 combined.
#' @param ph2_res results from phase 2 simulation. Optional, if supplied, will yield additional
#' summary results for phase 2&3 combined.
#' @param seed random seed for simulation
#' @param ... other parameters specific for the chosen phase 3 design
#' @returns summary statistics across all simulations
simulate_ph3 <- function(nsim,
                         data_list = NULL,
                         design = "RCT",
                         design_para = NULL,
                         median_surv0 = NULL,
                         median_surv1 = NULL,
                         p_trt = 0.5,
                         use_simsurv = FALSE,
                         maxt,
                         ph1_res = NULL,
                         ph2_res = NULL,
                         seed=NULL,
                         ...)
{
  # get all input arguments
  call <- match.call(expand.dots = TRUE)
  if (!is.null(seed)) set.seed(seed)
  
  source("data_generating_functions.R")
  # source the corresponding design method
  source(paste0("phase3/",design,".R"))
  
  # if user does not supply design parameters for the design method
  # then internally calls the design function
  if(is.null(design_para))
  {
    # get sample size/other design parameters from specified design
    default_args <- formals(paste0(design,"_design"))
    call_args <- as.list(call[-1]) # first one is the function name, not needed
    # change default values to user-input ones, if there is any
    final_args <- default_args
    final_args[names(call_args)] <- call_args
    # keep arguments relevant only
    final_args <- final_args[names(default_args)]
    design_para <- do.call(paste0(design,"_design"),final_args)
  } else
  {
    # check if the design parameters are complete
    names_design_para <- extract_return_names(paste0(design,"_design"))
    if(!identical(sort(names_design_para),sort(names(design_para))))
    {
      stop(paste0("Design parameters for ",design," are not complete!"))
    }
  }
  # get sample size
  # sample size
  N <- design_para$N
  
  # if data is supplied externally, check if the sample size N matches
  # the number of observations in the data frame
  if(!is.null(data_list))
  {
    nsim <- length(data_list)
    if(N > nrow(data_list[[1]][[1]]))
    {
      stop("The supplied data set has fewer observations than the design sample size, please input a larger data set!")
    }
    if(N < nrow(data_list[[1]][[1]]))
    {
      data_list <- sapply(data_list, function(data)
      {
        sapply(data, function(data_thisdose)
          {
          data_thisdose[1:N,]
        },simplify = FALSE)
      },simplify = FALSE)
      warning(paste0("The supplied data set has more observations than the design sample size, only the first ",N, " observations will be used!"))
    }
  } else {
  # generate list of data frames
  default_args <- formals(gen_surv)
  call_args <- sapply(call, function(i) eval(i,parent.frame()))[-1]
  final_args <- default_args
  final_args[names(call_args)] <- call_args
  final_args <- final_args[names(default_args)]
  final_args$N <- N
  data_list <- do.call(gen_surv,final_args)
  }
  
  # evaluate by the corresponding design
  # evaluate for each of the data frame
  # data_list is a list of length nsim 
  # every element is again a list of length(dose_levels)
  # data_list[[i]][[j]] stores the data frame for individual survival data
  # for the i-th iteration of the j-th dose level
  res <- sapply(data_list,function(data)
  {
    do.call(paste0(design,"_eval"),list(data = data, design_para = design_para))
  },simplify = FALSE)
  
  # prepare output
  # return some design features
  out <- list()
  out$phase <- "phase3"
  out$design <- design
  out$`sample size` <- design_para$N
  out$`number of iterations` <- nsim
  # summary table
  # see how many categories are returned by the design function
  res_names <- names(res[[1]])
  summary_tab <- data.frame(matrix(NA,
                                   ncol = length(res_names)+1,
                                   nrow = ncol(data_list[[1]][[1]])))
  colnames(summary_tab) <- c("dose",paste0("avg_",res_names))
  summary_tab$dose <- factor(1:ncol(data_list[[1]][[1]]))
  summary_tab[,paste0("avg_",res_names)]  = sapply(res_names,function(x) rowMeans(sapply(res,function(i) i[[x]])))
  out$`summary table` <- summary_tab
  return(out)
}

# Extract return variable names from a function
extract_return_names <- function(f) {
  body_expr <- body(f)  # Get function body
  return_lines <- body_expr[[length(body_expr)]]  # Extract the return() statement
  if (is.call(return_lines) && return_lines[[1]] == "return") {
    return_names <- names(eval(return_lines[[2]]))  # Extract list names
    return(return_names)
  }
  return(NULL)
}


