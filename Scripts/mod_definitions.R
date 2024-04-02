get_mod_specs <- function(mod = 1) {
  ### simulation specs
  
  n = c(1000)
  pk = 1
  nu_xy = 1
  nu_conf = c(1)
  treat_p = c(0.2, 0.5, 0.8)
  treat_beta = c(0.1, 0.5, 1)
  complexity = "linear"
  
  if(mod == 1) {
    type="onevar_linear"
    
  } else if (mod == 2) {
    type="multivar_linear"
    pk = c(100)
    nu_xy = c(0.1, 1)
    nu_conf = c(0.1, 1)
    
  } else if (mod == 3) {
    type="onevar_nonlinear"
    complexity="sin"
    
  } else if (mod == 4) {
    type= "multivar_nonlinear"
    pk = c(100)
    nu_xy = c(0.1, 1)
    nu_conf = c(0.1, 1)
    treat_p = c(0.2, 0.5, 0.8)
    complexity = "sin"
  }
  
  mod <- list(mod=mod, type=type, n=n, pk=pk, nu_xy=nu_xy, nu_conf=nu_conf, treat_p=treat_p, treat_beta=treat_beta, complexity=complexity)
  
  return(mod)
}

