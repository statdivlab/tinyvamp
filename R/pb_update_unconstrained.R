#
#
#
# pb_update_unconstrained <-
#   function(hessian_premultiplier, #H22n_inv%**H2dotn
#            theta_hat_unconstr, #mle
#            premult_unconstr_gradient_term, #H22n_inv %*% alpha root-n diff unconstr. gradients
#            theta_hat_constr, #mle for unconstrained variables
#            theta_hat_unconstr, #mle for constrained variables
#            curr_theta_constr #current value of constrained pars (in pb update)
#            ){
#   return(
#     -hessian_premultiplier%*% rbind(curr_theta_constr - theta_hat_constr,
#                                    -theta_hat_unconstr) -
#       premult_unconstr_gradient_term
#   )
# }
