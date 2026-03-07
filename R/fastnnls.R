fast_nnls <- function(
  ZTx,
  ZTZ,
  active = NULL,
  unconstrained = NULL,
  tolerance = 1e-10
) {
  # store number of covariates
  p <- nrow(ZTZ)

  if (length(intersect(active, unconstrained)) > 0) {
    stop(
      "A variable cannot be assigned to the unconstrained set and to the active set."
    )
  }
  ### getting stuck in this initial loop for some reason! check it out!
  if (!is.null(active)) {
    #set d = 0 to start
    d <- matrix(0, nrow = p, ncol = 1)
    #check appropriateness of initial active set
    s <- matrix(0, nrow = p, ncol = 1)
    not_active <- !(1:p %in% active)
    s_passive <- qr.solve(
      ZTZ[not_active, not_active, drop = F],
      ZTx[not_active, drop = F],
      tol = 1e-14
    )
    s[not_active, ] <- s_passive

    ### redefine not_active to excluded unconstrained variables (if included)
    if (!is.null(unconstrained)) {
      not_active <- !(1:p %in% union(active, unconstrained))
      s_passive <- s[not_active]
    }
    #Enter inner loop
    if (length(s_passive) > 0) {
      while (ifelse(length(s_passive) > 0, min(s_passive) <= 0, FALSE)) {
        alphas <- (-(d[not_active]) / (d[not_active] - s_passive))
        alphas <- alphas[s_passive < 0]
        alpha <- min(-alphas)

        d_archive <- d
        d <- d + alpha * (s - d)
        stopifnot(sum(is.na(d)) == 0)

        if (is.null(unconstrained)) {
          active <- (1:p)[d == 0]
          not_active <- !(1:p %in% active)
        }
        if (!is.null(unconstrained)) {
          active <- (1:p)[(d == 0) & (!((1:p) %in% unconstrained))]
          not_active <- !(1:p %in% active)
        }

        s <- matrix(0, nrow = p)
        if (length(active) > 0) {
          s_passive <- qr.solve(
            ZTZ[not_active, not_active],
            ZTx[not_active, ],
            tol = 1e-14
          )
          s[not_active, ] <- s_passive
        } else {
          s <- s_passive <- qr.solve(ZTZ, ZTx, tol = 1e-14)
        }

        if (!is.null(unconstrained)) {
          not_active <- !(1:p %in% union(active, unconstrained))
          s_passive <- s[not_active]
        }
      }
    }
    d <- s
    active <- (1:p)[d == 0]

    w <- ZTx - ZTZ %*% d
  }
  if (is.null(active)) {
    if (is.null(unconstrained)) {
      active <- 1:p
      d <- matrix(0, nrow = p, ncol = 1)
      w <- ZTx
    } else {
      active <- (1:p)[!(1:p %in% unconstrained)]
      d <- matrix(0, nrow = p, ncol = 1)
      d[unconstrained] <-
        qr.solve(
          ZTZ[unconstrained, unconstrained],
          ZTx[unconstrained, ],
          tol = 1e-14
        )
      w <- (ZTx - ZTZ %*% d)
    }
  }

  # counter <- 0
  #Enter main loop - ???? are loop conditions correct ???
  while ((ifelse(length(active) > 0, max(w[active]), -1) > tolerance)) {
    # counter <- counter + 1
    # print(counter)

    to_remove <- which.max(w[active])
    # print(to_remove)
    active <- active[-to_remove]

    if (length(active) > 0) {
      s <- matrix(0, nrow = p)
      not_active <- !(1:p %in% active)

      s_passive <- qr.solve(
        ZTZ[not_active, not_active],
        ZTx[not_active, ],
        tol = 1e-14
      )
      s[not_active, ] <- s_passive
    } else {
      s <- s_passive <- qr.solve(ZTZ, ZTx, tol = 1e-14)
    }

    ### redefine not_active to excluded unconstrained variables (if included)
    if (!is.null(unconstrained)) {
      not_active <- !(1:p %in% union(active, unconstrained))
      s_passive <- s[not_active]
    }

    #Enter inner loop
    if (length(s_passive) > 0) {
      while (ifelse(length(s_passive) > 0, min(s_passive) <= 0, FALSE)) {
        alphas <- (-(d[not_active]) / (d[not_active] - s_passive))
        alphas <- alphas[s_passive < 0]
        alpha <- min(-alphas)

        d_archive <- d
        d <- d + alpha * (s - d)
        stopifnot(sum(is.na(d)) == 0)

        if (is.null(unconstrained)) {
          active <- (1:p)[d == 0]
          not_active <- !(1:p %in% active)
        }
        if (!is.null(unconstrained)) {
          active <- (1:p)[(d == 0) & (!((1:p) %in% unconstrained))]
          not_active <- !(1:p %in% active)
        }

        s <- matrix(0, nrow = p)
        if (length(active) > 0) {
          s_passive <- qr.solve(
            ZTZ[not_active, not_active],
            ZTx[not_active, ],
            tol = 1e-14
          )
          s[not_active, ] <- s_passive
        } else {
          s <- s_passive <- qr.solve(ZTZ, ZTx, tol = 1e-14)
        }

        ### redefine not_active to excluded unconstrained variables (if included)
        if (!is.null(unconstrained)) {
          not_active <- !(1:p %in% union(active, unconstrained))
          s_passive <- s[not_active]
        }
      }
    }

    d <- s

    stopifnot(sum(is.na(d)) == 0)

    w <- (ZTx - ZTZ %*% d)

    # counter <- counter + 1

    # print(counter)
    # stopifnot(counter<13)
  }

  return(d)
}

fast_nnls_Matrix <- function(
  ZTx,
  ZTZ,
  active = NULL,
  unconstrained = NULL,
  tolerance = 1e-10
) {
  # store number of covariates
  p <- nrow(ZTZ)

  if (length(intersect(active, unconstrained)) > 0) {
    stop(
      "A variable cannot be assigned to the unconstrained set and to the active set."
    )
  }
  ### getting stuck in this initial loop for some reason! check it out!
  if (!is.null(active)) {
    #set d = 0 to start
    d <- matrix(0, nrow = p, ncol = 1)
    #check appropriateness of initial active set
    s <- matrix(0, nrow = p, ncol = 1)
    not_active <- !(1:p %in% active)
    s_passive <- qr.solve(
      ZTZ[not_active, not_active, drop = F],
      ZTx[not_active, drop = F],
      tol = 1e-14
    )
    s[not_active, ] <- s_passive

    ### redefine not_active to excluded unconstrained variables (if included)
    if (!is.null(unconstrained)) {
      not_active <- !(1:p %in% union(active, unconstrained))
      s_passive <- s[not_active]
    }
    #Enter inner loop
    if (length(s_passive) > 0) {
      while (ifelse(length(s_passive) > 0, min(s_passive) <= 0, FALSE)) {
        alphas <- (-(d[not_active]) / (d[not_active] - s_passive))
        alphas <- alphas[s_passive < 0]
        alpha <- min(-alphas)

        d_archive <- d
        d <- d + alpha * (s - d)
        stopifnot(sum(is.na(d)) == 0)

        if (is.null(unconstrained)) {
          active <- (1:p)[d == 0]
          not_active <- !(1:p %in% active)
        }
        if (!is.null(unconstrained)) {
          active <- (1:p)[(d == 0) & (!((1:p) %in% unconstrained))]
          not_active <- !(1:p %in% active)
        }

        s <- matrix(0, nrow = p)
        if (length(active) > 0) {
          s_passive <- qr.solve(
            ZTZ[not_active, not_active],
            ZTx[not_active, ],
            tol = 1e-14
          )
          s[not_active, ] <- s_passive
        } else {
          s <- s_passive <- qr.solve(ZTZ, ZTx, tol = 1e-14)
        }

        if (!is.null(unconstrained)) {
          not_active <- !(1:p %in% union(active, unconstrained))
          s_passive <- s[not_active]
        }
      }
    }
    d <- s
    active <- (1:p)[d == 0]

    w <- ZTx - ZTZ %*% d
  }
  if (is.null(active)) {
    if (is.null(unconstrained)) {
      active <- 1:p
      d <- matrix(0, nrow = p, ncol = 1)
      w <- ZTx
    } else {
      active <- (1:p)[!(1:p %in% unconstrained)]
      d <- matrix(0, nrow = p, ncol = 1)
      d[unconstrained] <-
        qr.solve(
          ZTZ[unconstrained, unconstrained],
          ZTx[unconstrained, ],
          tol = 1e-14
        )
      w <- (ZTx - ZTZ %*% d)
    }
  }

  # counter <- 0
  #Enter main loop - ???? are loop conditions correct ???
  while ((ifelse(length(active) > 0, max(w[active]), -1) > tolerance)) {
    # counter <- counter + 1
    # print(counter)

    to_remove <- which.max(w[active])
    # print(to_remove)
    active <- active[-to_remove]

    if (length(active) > 0) {
      s <- matrix(0, nrow = p)
      not_active <- !(1:p %in% active)

      s_passive <- qr.solve(
        ZTZ[not_active, not_active],
        ZTx[not_active, ],
        tol = 1e-14
      )
      s[not_active, ] <- s_passive
    } else {
      s <- s_passive <- qr.solve(ZTZ, ZTx, tol = 1e-14)
    }

    ### redefine not_active to excluded unconstrained variables (if included)
    if (!is.null(unconstrained)) {
      not_active <- !(1:p %in% union(active, unconstrained))
      s_passive <- s[not_active]
    }

    #Enter inner loop
    if (length(s_passive) > 0) {
      while (ifelse(length(s_passive) > 0, min(s_passive) <= 0, FALSE)) {
        alphas <- (-(d[not_active]) / (d[not_active] - s_passive))
        alphas <- alphas[s_passive < 0]
        alpha <- min(-alphas)

        d_archive <- d
        d <- d + alpha * (s - d)
        stopifnot(sum(is.na(d)) == 0)

        if (is.null(unconstrained)) {
          active <- (1:p)[d == 0]
          not_active <- !(1:p %in% active)
        }
        if (!is.null(unconstrained)) {
          active <- (1:p)[(d == 0) & (!((1:p) %in% unconstrained))]
          not_active <- !(1:p %in% active)
        }

        s <- matrix(0, nrow = p)
        if (length(active) > 0) {
          s_passive <- qr.solve(
            ZTZ[not_active, not_active],
            ZTx[not_active, ],
            tol = 1e-14
          )
          s[not_active, ] <- s_passive
        } else {
          s <- s_passive <- qr.solve(ZTZ, ZTx, tol = 1e-14)
        }

        ### redefine not_active to excluded unconstrained variables (if included)
        if (!is.null(unconstrained)) {
          not_active <- !(1:p %in% union(active, unconstrained))
          s_passive <- s[not_active]
        }
      }
    }

    d <- s

    stopifnot(sum(is.na(d)) == 0)

    w <- (ZTx - ZTZ %*% d)

    # counter <- counter + 1

    # print(counter)
    # stopifnot(counter<13)
  }

  return(d)
}
