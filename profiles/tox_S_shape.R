logistic_fun <- function(d,alpha,beta)
{
  L = alpha+beta*d
  return(exp(L)/(1+exp(L)) - exp(alpha)/(1+exp(alpha)))
}

tox_S_shape <- function(d)
{
  logistic_fun(d,alpha = -4, beta = 0.008)
}
attr(tox_S_shape,"label") <- "tox_S_shape"
