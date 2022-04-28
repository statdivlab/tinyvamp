simpl_opt_linesearch_fnnls <- function(fn, #objective function to be minimized
                                       x, #starting values of simplex-constrained parameter
                                       xhess, # hessian matrix of objective function at x
                                       xgrad, # gradient of objective function at x
                                       lambda = 0, #trust penalty
                                       maxit = 100,
                                       constraint_tolerance = 1e-10
){

  x0 <- x
  curr_fn_value <- fn(x0)


  new_x <-  simpl_auglag_fnnls(x = x,
                               fn = fn,
                               xhess = xhess,
                               xgrad = xgrad,
                               lambda = lambda,
                               constraint_tolerance = constraint_tolerance)

  new_fn_value <- fn(new_x)

  step_direction <-  new_x - x0

  stepsize <- 1

  fn_decrease <- 10

  while((fn_decrease >  0)&(stepsize>1e-2)){


    prop_x <- x0 + stepsize*step_direction

    new_fn_value <- fn(prop_x)
    if(is.nan(new_fn_value)){
      new_fn_value <- 1e100
    }

    if(min(prop_x)<0){
      new_fn_value <- 1e100
    }
    #
    fn_decrease <- new_fn_value - curr_fn_value
    stepsize <- stepsize/2

  }

  if(fn_decrease <= 0){
    return(prop_x)
  } else{
    return(x0)
  }


}
