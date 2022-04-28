#
# sum_of_logs <- function(x){ #calculates log(sum(exp(x)))
#   if(length(x) ==1){
#     return(x)
#   } else{
#     x <- x[order(x, decreasing = TRUE)]
#     return(x[1] + log(1 + sum(exp(x[2:length(x)] - x[1]))))}
# }
#
# expit <- function(x){
#   xmaxed <- x-max(x)
#   return(exp(xmaxed) - sum_of_logs(c(-max(x),xmaxed)))
# }
#
# abspit <- function(x){
#   xmaxed <- x/max(x)
#   return(abs(xmaxed)/(sum(abs(xmaxed))))
# }
#
#
# diff_of_logs <- function(larger, smaller){
#   return(larger + log( 1 - exp(smaller - larger) ) )
# }
#
# signed_log_sum <- function(log_magnitudes,
#                            signs){
#   # if no nonzero summands, return zero
#   if(sum(signs != 0) == 0){
#     return(c(0,0))
#   }
#
#   # if no negative summands, return sum of positive summands
#   if(sum(signs == -1) == 0){
#     return(c(sum_of_logs(log_magnitudes[signs == 1]),1))
#   }
#
#   # if no positive summands, return sum of negative summands (with sign -1)
#   if(sum(signs == 1) == 0){
#     return(c(sum_of_logs(log_magnitudes[signs == -1]),-1))
#   }
#
#   # otherwise sum positive and negative magnitudes separately
#   positives <- log_magnitudes[signs == 1]
#   negatives <- log_magnitudes[signs == -1]
#
#   pos_logsum <- sum_of_logs(positives)
#   neg_logsum <- sum_of_logs(negatives)
#
#   # if positive magnitude equals negative magnitude, return zero
#   if(pos_logsum == neg_logsum){
#     return(c(0,0))
#   }
#
#   # otherwise return difference of positive and negative terms with appropriate sign
#   return(c(diff_of_logs(max(pos_logsum,neg_logsum),
#                         min(pos_logsum, neg_logsum)),
#            ifelse(pos_logsum > neg_logsum, 1, -1)))
# }
#
# mixture_dist <- data.frame(support = 1:5,
#                            weights = rep(1/5,5))
#
#
# ### Tests for this
# # - point mixtures
# # - simple mixtures (2 comps)
# poisson_lagrange <- function(y,mixture_dist,
#                              return_sum = TRUE,
#                              absolute = FALSE){
#   n <- length(y)
#   S <- nrow(mixture_dist)
#
#   lls <- matrix(nrow = n, ncol = S)
#
#   for(i in 1:n){
#     for(s in 1:S){
#       lls[i,s] <- dpois(y[i], lambda = exp(mixture_dist$support[s]), log = T)+ log(abs(mixture_dist$weights[s]))
#     }
#   }
#
#   weight_signs <- sign(mixture_dist$weights)
#   marg_lls <- numeric(n)
#   for(i in 1:n){
#     to_add <- signed_log_sum(lls[i,],weight_signs)
#     if(to_add[2] != 1){
#       return(1e10)
#     }
#     marg_lls[i] <- to_add[1]
#   }
#   if(return_sum){
#     to_return <- -(1/n)*sum(marg_lls) + ifelse(absolute,sum(abs(mixture_dist$weights)),
#                                                sum(mixture_dist$weights))
#     if(is.infinite(to_return)){
#       return(1e10)
#     } else{
#       return(to_return)
#     }
#   } else{
#     return(marg_lls)
#   }
# }
#
#
# ### Tests
# # - vs numerical differentiation
#
# poisson_lagrange_dir_deriv <- function(y,
#                                        mixture_dist,
#                                        new_supp = NULL,
#                                        absolute = FALSE){
#
#   n <- length(y)
#   if(is.null(new_supp)){
#     new_supp <- mixture_dist$support
#     use_mix_supp <- TRUE
#   } else{
#     use_mix_supp <- FALSE
#   }
#   S <- length(new_supp)
#
#   marg_lls <- poisson_lagrange(y,mixture_dist,
#                                return_sum = FALSE)
#
#   cond_lls <- matrix(nrow = n, ncol = S)
#
#
#   for(i in 1:n){
#     for(s in 1:S){
#       cond_lls[i,s] <- -exp(dpois(y[i],exp(new_supp[s]),log  =TRUE) - marg_lls[i] - log(n))
#     }
#   }
#
#   if(use_mix_supp&absolute){
#   constraint_term <- sign(mixture_dist$weights)
#   } else{
#     constraint_term <- 1
#   }
#
#     grad <- apply(cond_lls, 2, sum) + constraint_term
#
#   return(grad)
# }
# set.seed(0)
# y <- rpois(100, lambda = 400)
#
# mixture_dist <- data.frame("support"  = c(log(300),log(400)),
#                            "weights" = c(0.5,0.5))
#
# # weights in (0,1)
# numDeriv::grad(
#                function(x) poisson_lagrange(y, mixture_dist = data.frame("support" = mixture_dist$support,
#                                               "weights" = x)),
#                mixture_dist$weights)
# poisson_lagrange_dir_deriv(y,mixture_dist)
#
# # positive and negative weights
#
# mixture_dist$weights <- c(-1,2)
#
# numDeriv::grad(
#   function(x) poisson_lagrange(y, mixture_dist = data.frame("support" = mixture_dist$support,
#                                                             "weights" = x)),
#   mixture_dist$weights)
# poisson_lagrange_dir_deriv(y,mixture_dist)
#
#
#
# # big ol' weights
# mixture_dist$weights <- c(-1e10,1e10 + 1)
#
# numDeriv::grad(
#   function(x) poisson_lagrange(y, mixture_dist = data.frame("support" = mixture_dist$support,
#                                                             "weights" = x)),
#   mixture_dist$weights)
# poisson_lagrange_dir_deriv(y,mixture_dist)
#
#
#
#
#
#
# poisson_ll_dir_deriv_log <- function(y,
#                                        mixture_dist,
#                                        new_supp = NULL){
#
#   n <- length(y)
#   if(is.null(new_supp)){
#     new_supp <- mixture_dist$support
#     use_mix_supp <- TRUE
#   } else{
#     use_mix_supp <- FALSE
#   }
#   S <- length(new_supp)
#
#   marg_lls <- poisson_lagrange(y,mixture_dist,
#                                return_sum = FALSE)
#
#   cond_lls <- matrix(nrow = n, ncol = S)
#
#
#   for(i in 1:n){
#     for(s in 1:S){
#       cond_lls[i,s] <- dpois(y[i],exp(new_supp[s]),log  =TRUE) - marg_lls[i]
#     }
#   }
#
#
#     grad <- apply(cond_lls, 2, sum_of_logs)
#
#     return(grad)
#
# }
#
# poisson_ll_dir_deriv_log_dsupp <-  function(y,
#                                             mixture_dist,
#                                             new_supp){
#   n <- length(y)
#   S <- length(new_supp)
#
#   cond_lls <- matrix(nrow = n, ncol = S)
#
#   marg_lls <- poisson_lagrange(y, mixture_dist, return_sum = F)
#
#   for(i in 1:n){
#     for(s in 1:S){
#       cond_lls[i,s] <- dpois(y[i], exp(new_supp[s]), log = TRUE) - marg_lls[i]
#     }
#   }
#
#   grads <- numeric(S)
#
#   log_derivs <- poisson_ll_dir_deriv_log(y,
#                                          mixture_dist,
#                                          new_supp)
#
#   for(s in 1:S){
#     diffs <- y - exp(new_supp[s])
#
#     grads[s] <- sum(exp(cond_lls[,s] - log_derivs[s])*diffs)
#
#   }
#
#   return(grads)
#
# }
#
# # numgrad <- numDeriv::grad(function(x) poisson_ll_dir_deriv_log(y,mixture_dist,
# #                                                                new_supp = x),
# #                           x = 5:7)
# #
# # poisson_ll_dir_deriv_log_dsupp(y, mixture_dist,
# #                                         5:7)
#
# # poisson_lagrange_supp_deriv <- function(y,
# #                                         mixture_dist){
# #
# #   n <- length(y)
# #   new_supp <- mixture_dist$support
# #
# #   S <- length(new_supp)
# #
# #   marg_lls <- poisson_lagrange(y,mixture_dist,
# #                                return_sum = FALSE)
# #
# #   cond_lls <- matrix(nrow = n, ncol = S)
# #   ll_derivs <- matrix(nrow = n, ncol = S)
# #
# #
# #   for(i in 1:n){
# #     for(s in 1:S){
# #       cond_lls[i,s] <- -exp(dpois(y[i],exp(new_supp[s]),log  =TRUE) - marg_lls[i] - log(n))*mixture_dist$weights[s]
# #       ll_derivs[i,s] <- y[i] - exp(new_supp[s])
# #     }
# #   }
# #
# #
# #   grad <- apply(cond_lls*ll_derivs, 2, sum)
# #
# #   return(grad)
# #
# # }
#
# # numgrad <- grad(function(t) poisson_lagrange(y, mixture_dist = data.frame(weights = t, support = 1:5)),
# #                 x = rep(1/5,5))
# #
# # analgrad <- poisson_lagrange_dir_deriv(y, mixture_dist = data.frame(weights = rep(1/5,5), support = 1:5))$grad
# #
# # numgrad <- grad(function(t) poisson_lagrange(y, mixture_dist = data.frame(weights = 1/5, support = t)),
# #                 x = 1:5)
# #
# # analgrad <- poisson_lagrange_supp_deriv(y, mixture_dist = data.frame(weights = rep(1/5,5), support = 1:5))
# #
# #
# # numhess <- hessian(function(t) poisson_lagrange(y, mixture_dist = data.frame(weights = t, support = 1:5)),
# #                    x = rep(1/5,5))
# #
# # analhess <-  poisson_lagrange_dir_deriv(y, mixture_dist = data.frame(weights = rep(1/5,5), support = 1:5),
# #                                         new_supp = 1:5,
# #                                         hessian = TRUE)
#
#
# par(mfrow = c(1,2))
# nsim <- 3
# for_sims <- vector(nsim, mode = "list")
# sims <- vector(3,mode = "list")
#
# for(k in 1:3){
#   sims[[k]] <- for_sims
# }
#
# set.seed(0)
# ords <- 1:3
# for(ord in ords){
# for(sim in 1:nsim){
#   print(sim)
#   supp <- seq(1,10, by = 1)
#   mixture_dist <- data.frame(weights = sqrt(1/length(supp)),
#                              support = supp)
#
#   # y <- sapply(1:100,
#   #             function(k){ rpois(1,lambda = ifelse(sample(c(T,F),1),
#   #                                                  exp(rnorm(1,sample(c(5,10),1))),
#   #                                                  exp(7.5)))})
#
#   y <- replicate(10^ord,sample(500:1000,1))
#   y_stor <- y
#
# # keep_going <- TRUE
# #   while(keep_going){
#
# curr_derivs <- 10
# while(sum(curr_derivs^2)> 0.0001){
#
#   absolute <- sum(curr_derivs)> 10
#
# ### Looks like this optimization is unstable for some reason -- INVESTIGATE!
#   attempt <- optim(par = mixture_dist$weights,
#                    function(x) poisson_lagrange(y, mixture_dist = data.frame(weights = x,
#                                                                              support = mixture_dist$support),
#                                                 absolute = absolute),
#                    gr = function(x) poisson_lagrange_dir_deriv(y,
#                                                                mixture_dist = data.frame(weights = x,
#                                                                                          support = mixture_dist$support),
#                                                                absolute = absolute),
#                    method = "BFGS",
#                    control = list(
#                      # maxit = 5,
#                                   trace = 6)
#   )
#
#   # print(attempt$value)
#
#   mixture_dist$weights <-   attempt$par
#
#   mixture_dist <- mixture_dist[mixture_dist$weights>0,]
#
#   # mixture_dist$weights <- mixture_dist$weights/sum(mixture_dist$weights)
#
#
#
#     mixture_dist <- mixture_dist[order(mixture_dist$support),]
#     S <- nrow(mixture_dist)
#
#     candidate_supp <- seq(-10,10, by = .01)
#
#     grs <- poisson_lagrange_dir_deriv(y,mixture_dist,
#                                       new_supp = candidate_supp)
#
#
#     plot(candidate_supp, grs, type = "l")
#
#     candidate_supp <- candidate_supp[grs <0]
#     grs <- grs[grs<0]
#
#     accepted_supp <- numeric(0)
#     if(length(candidate_supp)>0){
#     for(s in 1:(S + 1)){
#       lb <- ifelse(s >1, mixture_dist$support[s - 1],-Inf)
#       ub <- ifelse(s < S + 1, mixture_dist$support[s], Inf)
#
#       supp_filt <- (candidate_supp > lb)&(candidate_supp < ub)
#       to_add <- (candidate_supp[supp_filt])[which.min(grs[supp_filt])]
#       accepted_supp <- c(accepted_supp, to_add)
#
#     }
#
#     }
#
#     mixture_dist <- rbind(mixture_dist,
#                           data.frame("support" = accepted_supp,
#                                      "weights" = rep(1e-3, length(accepted_supp))))
#     mixture_dist$weights <- mixture_dist$weights/sum(mixture_dist$weights)
#
#
#     # lbs <- c(mixture_dist$support[1] - 1,
#     #          mixture_dist$support)
#     #
#     # ubs <- c(mixture_dist$support,
#     #          mixture_dist$support[S] + 1)
#     # starting <- apply(cbind(lbs,ubs),1,mean)
#     # all_together <- optim(par = starting,
#     #                       function(x) -sum(poisson_ll_dir_deriv_log(y,
#     #                                                                 mixture_dist,
#     #                                                                 new_supp = x)),
#     #                       gr = function(x) -poisson_ll_dir_deriv_log_dsupp(y,
#     #                                                           mixture_dist,
#     #                                                           x),
#     #                       method = "L-BFGS-B",
#     #                       lower = lbs,
#     #                       upper = ubs)
#     #
#     # supp_candidates <- all_together$par
#     #
#     # grads <- poisson_lagrange_dir_deriv(y,
#     #                            mixture_dist,
#     #                            new_supp = supp_candidates)
#     #
#     #
#     # gr_condition<- grads < 0
#     #
#     # supp_condition <- sapply(supp_candidates,
#     #                               function(x) min(abs(x - mixture_dist$support))>1e-8)
#     #
#     # supp_candidates <- supp_candidates[gr_condition&supp_condition]
#     #
#     # if(length(supp_candidates)>0){
#     # mixture_dist <- rbind(mixture_dist,
#     #                       data.frame("support"= supp_candidates,
#     #                                  "weights" = 1e-3))
#     # mixture_dist$weights <- mixture_dist$weights/sum(mixture_dist$weights)
#     # } else{
#     #   keep_going <- FALSE
#     # }
#
#     curr_derivs <- poisson_lagrange_dir_deriv(y,
#                                mixture_dist)
#
#     print(sum(curr_derivs^2))
#
#     t <- seq(6,7,by = .001)
#     derivs <- poisson_lagrange_dir_deriv(y, mixture_dist,
#                                          t)
#     par(mfrow = c(1,2))
#     plot(t, derivs, type = "l")
#     points(mixture_dist$support, rep(0, length(mixture_dist$support)),
#            pch = "|")
#     abline(a = 0, b = 0, lty = 2)
#
#     mixture_dist <- mixture_dist[order(mixture_dist$support),]
#     S <- nrow(mixture_dist)
#     plot(stepfun(x = mixture_dist$support, y= c(0,sapply(1:S, function(k) sum(mixture_dist$weights[1:k])))))
#
#
#
#   }
#
#   sims[[ord]][[sim]] <- mixture_dist
#
#   print(mixture_dist)
# }
# }
#
# t <- seq(6,7,by = .001)
# derivs <- poisson_lagrange_dir_deriv(y, mixture_dist,
#                            t)$grad
# par(mfrow = c(1,2))
# plot(t, derivs, type = "l")
# points(mixture_dist$support, rep(0, length(mixture_dist$support)),
#        pch = "|")
# abline(a = 0, b = 0, lty = 2)
#
# mixture_dist <- mixture_dist[order(mixture_dist$support),]
# S <- nrow(mixture_dist)
# plot(stepfun(x = mixture_dist$support, y= c(0,sapply(1:S, function(k) sum(mixture_dist$weights[1:k])))))
