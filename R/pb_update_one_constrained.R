#
# pb_update_one_constrained <-
#   function(constr_gradient_term, # alpha root-n grad l-star - grad l
#            H2dot_transp_theta_hat,
#            #H2dot_n^Ttheta^hat %*%
#            #rbind(curr_theta_const - theta_const_mle,
#            #-theta_var_mle)
#            # where theta_var is (simplex-constrained) parameter to be optimized
#            # over, and theta_const is all other parameters (held constant in
#            # this step)
#            H22n, #submatrix of criterion hessian corresp. to theta_var
#            curr_theta_var,
#            curr_theta_const,
#            theta_hat_var,
#            theta_hat_const
#   ){
#
#     #directly use auglag - no simpl_auglag_fnnls
#     # prox_criterion <- function(theta_var){
#     #   theta_diff <- ( rbind(curr_theta_const,
#     #                   x) -
#     #               rbind(theta_hat_const,
#     #                     theta_hat_var))
#     #   gr_term <- constr_gradient_term %*%theta_diff
#     #   hess_term <- t(theta_diff)%*%H22n%*%theta_diff
#     #   return(as.numeric(gr_term + 0.5*hess_term))
#     # }
#     #
#     #
#     # simpl_auglag_fnnls(x = curr_theta_var,
#     #                    fn = prox_criterion, #function of x to optimize
#     #                    xhess = H22n, #hessian at x
#     #                    xgrad = , #gradient at x
#     #                                lambda, #penalty parameters
#     #                                nu = 1, #starting lagrangian penalty
#     #                                mu = 1, #starting augmented lagrangian penalty
#     #                                constraint_tolerance = 1e-10, #sum-to-one constraint tolerance
#     #                                maxit = 100 # maximum number of iterations (outer loop)
#     #
#     # )
#   }
