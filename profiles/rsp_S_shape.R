logistic_fun <- function(d,alpha,beta)
{
  L = alpha+beta*d
  return(exp(L)/(1+exp(L)) - exp(alpha)/(1+exp(alpha)))
}

rsp_S_shape <- function(d)
{
  logistic_fun(d,alpha=-3,beta=0.01)/1.7
}

attr(rsp_S_shape,"label") <- "rsp_S_shape"

