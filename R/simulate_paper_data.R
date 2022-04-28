

simulate_paper_data <- function(n,
                                J,
                                B_multiplier,
                                distrib,
                                seed){

if(!(J %in% c(5,20))){
  stop("Number of taxa J must be equal to 5 or 20 in this simulation")
}

if(!(n %in% c(1,3))){
  stop("Number of technical replicates n must be 1 or 3 in this simulation")
}

  ### generate Z_tilde (alpha_tilde_k = 0 for k = 1, 2)
Z_tilde <- do.call(rbind,lapply(1:n,
                                function(x) matrix(rep(c(1,9,81,729),
                                                       4),
                                                   ncol = 1)))

 ### generate Z
Z <- do.call(rbind,lapply(1:4,
                    function(x) do.call(rbind,
                                        lapply(1:(4*n),function(k) matrix(
                                          as.numeric(x == 1:4),nrow = 1
                                        )))))

X <- matrix(1,nrow = n*16,ncol = 1)

B_star <- matrix(c(
  rep(c(3,-1,1,-3),ceiling(J/4))[1:(J - 1)],
                   0),nrow = 1,ncol = J)

B <- B_star*B_multiplier

P1 <- 2^(seq(0,4,length.out = J))
P1 <- P1/sum(P1)
P2 <- P1[J:1]
P3 <- c(rep(0,0.4*J),
        10^seq(0,2,length.out = 0.6*J))
P3 <- P3/sum(P3)
P4 <- P3[J:1]

P <- rbind(P1,P2,P3,P4)

P_tilde <- matrix(1/J,ncol = J, nrow = 1)
X_tilde <- matrix(1,ncol = 1, nrow = 1)
gamma_tilde <- -3.7

dilutions <- log(as.numeric(Z_tilde))/log(9)

gamma_means <- sapply(dilutions, function(d) min(13.5 - 1.5*d,12) )
set.seed(seed)
gammas <- sapply(gamma_means, function(x) rnorm(1,mean = x, sd = sqrt(0.05)))


means <- meaninate(gammas = gammas,
                   B = B,
                   X= X,
                   Z = Z,
                   P = P,
                   X_tilde = X_tilde,
                   Z_tilde = Z_tilde,
                   Z_tilde_gamma_cols = ncol(Z_tilde),
                   P_tilde = P_tilde,
                   gamma_tilde = gamma_tilde)


if(distrib== "Poisson"){
  W <- apply(means, c(1,2), function(x) rpois(1,lambda = x))
}

if(distrib == "NB"){
  W <- apply(means, c(1,2), function(x) rnbinom(1,mu = x, size = 13))
}

return(W)

}
