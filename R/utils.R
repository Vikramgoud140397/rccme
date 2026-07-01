#' Validate inputs for rccme_calib_me
#'
#' @inheritParams rccme_calib_me
#' @return Returns cleaned obects if all checks pass.
#' @keywords internal
.validate_calib_me <- function(
    w_mat, rel_vec = NULL, w_se_vec = NULL, w_se_mat = NULL, z_mat) {
  if (!is.matrix(w_mat)) {
    w_mat <- tryCatch(as.matrix(w_mat), error = function(e) NULL)
  }
  if (is.null(w_mat)) {
    stop("Unable to transform `w_mat` into a matrix.")
  }

  if (!is.matrix(z_mat)) {
    z_mat <- tryCatch(as.matrix(z_mat), error = function(e) NULL)
  }
  if (is.null(z_mat)) {
    stop("Unable to transform `z_mat` into a matrix.")
  }

  if (!is.numeric(w_mat)) {
    stop("`w_mat` must be a matrix of numbers.")
  }
  if (!is.numeric(z_mat)) {
    stop("`z_mat` must be a matrix of numbers.")
  }

  is_null_rel <- is.null(rel_vec)
  is_null_wse <- is.null(w_se_vec)
  is_null_wma <- is.null(w_se_mat)
  if (is_null_rel && is_null_wse && is_null_wma) {
    stop("Provide one of `rel_vec` `w_se_vec` or `w_se_mat`.")
  }
  if ((
    as.integer(!is_null_rel) + as.integer(!is_null_wse) +
      as.integer(!is_null_wma)
  ) > 1) {
    stop("Provide only one of `rel_vec`, `w_se_vec` or `w_se_mat`.")
  }

  n <- nrow(w_mat)
  p <- ncol(w_mat)
  q <- ncol(z_mat)

  if (!is_null_rel) {
    if (!is.numeric(rel_vec)) {
      stop("`rel_vec` must be a vector of numbers.")
    }

    if (length(rel_vec) != p) {
      stop(paste(
        "`rel_vec` have the same number of elements are there",
        "are columns in w_mat.",
        sep = " "
      ))
    }

    if (any(!is.finite(rel_vec))) {
      stop("`rel_vec` contains values that are not valid numbers.")
    }

    if (any(rel_vec <= 0) || any(rel_vec >= 1)) {
      stop("All values in `rel_vec` must be between 0 and 1.")
    }

    # No missing allowed in w_mat
    if (anyNA(w_mat)) {
      stop(paste(
        "`w_mat` cannot contain missing values under",
        "homoskedastic measurement error.",
        sep = " "
      ))
    }
  }

  if (!is_null_wse) {
    if (!is.numeric(w_se_vec)) {
      stop("`w_se_vec` must be a vector of numbers.")
    }

    if (length(w_se_vec) != p) {
      stop(paste(
        "`w_se_vec` have the same number of elements are there",
        "are columns in w_mat.",
        sep = " "
      ))
    }

    if (any(!is.finite(w_se_vec))) {
      stop("`w_se_vec` contains values that are not valid numbers.")
    }

    if (any(w_se_vec <= 0)) {
      stop("All values in `w_se_vec` must be positive.")
    }

    # No missing allowed in w_mat
    if (anyNA(w_mat)) {
      stop(paste(
        "`w_mat` cannot contain missing values under",
        "homoskedastic measurement error.",
        sep = " "
      ))
    }
  }

  if (!is_null_wma) {
    if (!is.matrix(w_se_mat)) {
      w_se_mat <- tryCatch(as.matrix(w_se_mat),
        error = function(e) NULL
      )
    }
    if (is_null_wma || !is.numeric(w_se_mat)) {
      stop("`w_se_mat` must be a matrix of numbers.")
    }

    if (!all(dim(w_se_mat) == dim(w_mat))) {
      stop("`w_se_mat` must have same dimensions as `w_mat`.")
    }

    na_w <- is.na(w_mat)
    na_se <- is.na(w_se_mat)

    if (!identical(unname(na_w), unname(na_se))) {
      stop("Missingness in `w_mat` and `w_se_mat` must match.")

      
          }

    if (any(na_w)) {
      w_mat[na_w] <- 0
      w_se_mat[na_se] <- 1
    }

    if (any(!is.finite(w_se_mat))) {
      stop("`w_se_mat` contains values that are not valid numbers.")
    }

    if (any(w_se_mat <= 0)) {
      stop("All entries of `w_se_mat` must be positive.")
    }

    if (any(!is.finite(w_mat))) {
      stop("`w_mat` contains values that are not valid numbers.")
    }
  }

  if (any(is.infinite(z_mat))) {
    stop("`z_mat` contains infinite values.")
  }

  if (nrow(z_mat) != n) {
    stop("`z_mat` must have the same number of rows as `w_mat`.")
  }

  # Sample size check
  if (n < p) {
    stop(paste(
      "There must be more cases than variables to compute in `w_mat`",
      "to compute valid covariance matrices.",
      sep = " "
    ))
  }
  if (n < q) {
    stop(paste(
      "There must be more cases than variables to compute in `z_mat`",
      "to compute valid covariance matrices.",
      sep = " "
    ))
  }

  # Constant column checks (can cause singular covariance)
  if (any(apply(w_mat, 2, stats::var) == 0)) {
    stop(paste(
      "`w_mat` has constant columns;",
      "Unable to compute valid covariance matrices.",
      sep = " "
    ))
  }

  if (ncol(z_mat) > 0 && any(apply(z_mat, 2, stats::var) == 0)) {
    stop(paste(
      "`z_mat` has constant columns;",
      "Unable to compute valid covariance matrices.",
      sep = " "
    ))
  }

  list(
    w_mat = w_mat,
    rel_vec = rel_vec,
    w_se_vec = w_se_vec,
    w_se_mat = w_se_mat,
    z_mat = z_mat
  )
}

#' Prepare sum-score, EAP or regression factor scores for calibration
#'
#' @inheritParams rccme_calib_me
#' @return Returns rescaled scores.
#' @keywords internal
.prep_score <- function(
    w_mat,
    rel_vec = NULL, w_se_vec = NULL, w_se_mat = NULL, standard = FALSE) {
  if (!is.null(rel_vec)) {
    sd_s <- apply(w_mat, 2, stats::sd)
    if (isTRUE(standard)) {
      w_mat <- sweep(w_mat, 2, sd_s * sqrt(rel_vec), "/")
      sd_s <- 1.0 / sqrt(rel_vec)
    }
    rel_ests <- rep(1.0, length(rel_vec))
    w_se_vec <- sd_s * sqrt(1 - rel_vec)
  } else if (!is.null(w_se_vec)) {
    rel_ests <- 1 - w_se_vec^2
    w_se_vec <- w_se_vec / sqrt(rel_ests)
  } else {
    rel_ests <- colMeans(1 - w_se_mat^2)
    w_se_mat <- sweep(w_se_mat, 2, sqrt(rel_ests), "/")
  }
  w_mat <- sweep(w_mat, 2, rel_ests, "/")

  list(
    w_mat = w_mat,
    w_se_vec = w_se_vec,
    w_se_mat = w_se_mat
  )
}

#' Get conditional expectation under homoskedastic measurement error
#'
#' @inheritParams rccme_calib_me
#' @return Returns cleaned obects if all checks pass.
#' @keywords internal
.get_e_z_zw_hom <- function(
    w_mat, w_se_vec, z_mat) {
  n <- nrow(w_mat)
  p <- ncol(w_mat)
  q <- ncol(z_mat)

  sigma_ww <- stats::cov(w_mat)
  mu_w <- colMeans(w_mat)

  mu_z <- colMeans(z_mat)
  sigma_zz <- stats::cov(z_mat)
  sigma_xz <- stats::cov(w_mat, z_mat)
  if (anyNA(mu_z)) {
    ret_cov <- .optim_cov_fiml(z_mat)
    mu_z <- ret_cov$mu
    sigma_zz <- ret_cov$c_mat
    sigma_xz <- .cov_xy_fiml(w_mat, z_mat, mu_w, mu_z, sigma_zz)
  }

  exp_x <- matrix(nrow = n, ncol = p)
  mat_stable <- matrix(nrow = p + q, ncol = p + q)
  if (q > 0) {
    mat_stable[(p + 1):(p + q), seq_len(p)] <- t(sigma_xz)
    mat_stable[seq_len(p), (p + 1):(p + q)] <- sigma_xz
    mat_stable[(p + 1):(p + q), (p + 1):(p + q)] <- sigma_zz
  }
  sigma_uu <- diag(w_se_vec^2, p)
  sigma_xx <- sigma_ww - sigma_uu
  mat_stable[seq_len(p), seq_len(p)] <- sigma_ww
  inv_mat <- chol2inv(chol(mat_stable))

  mat_xxxz <- cbind(sigma_xx, sigma_xz)

  exp_x <- mu_w +
    mat_xxxz %*% inv_mat %*% rbind(t(w_mat) - mu_w, t(z_mat) - mu_z)
  t(exp_x)
}

#' Get conditional expectation under heteroskedastic measurement error
#'
#' @inheritParams rccme_calib_me
#' @return Returns cleaned obects if all checks pass.
#' @keywords internal
.get_e_z_zw_het <- function(
    w_mat, w_se_mat, z_mat) {
  n <- nrow(w_mat)
  p <- ncol(w_mat)
  q <- ncol(z_mat)

  sigma_xx <- stats::cov(w_mat)
  diag(sigma_xx) <- diag(sigma_xx) - colMeans(w_se_mat^2)

  mu_w <- colMeans(w_mat)

  mu_z <- colMeans(z_mat)
  sigma_zz <- stats::cov(z_mat)
  sigma_xz <- stats::cov(w_mat, z_mat)
  if (anyNA(mu_z)) {
    ret_cov <- .optim_cov_fiml(z_mat)
    mu_z <- ret_cov$mu
    sigma_zz <- ret_cov$c_mat
    sigma_xz <- .cov_xy_fiml(w_mat, z_mat, mu_w, mu_z, sigma_zz)
  }

  exp_x <- matrix(nrow = n, ncol = p)
  mat_stable <- matrix(nrow = p + q, ncol = p + q)
  if (q > 0) {
    mat_stable[(p + 1):(p + q), seq_len(p)] <- t(sigma_xz)
    mat_stable[seq_len(p), (p + 1):(p + q)] <- sigma_xz
    mat_stable[(p + 1):(p + q), (p + 1):(p + q)] <- sigma_zz
  }
  mat_xxxz <- cbind(sigma_xx, sigma_xz)

  for (i in 1:n) {
    sigma_ww_i <- sigma_xx + diag(w_se_mat[i, ]^2, p)
    inv_mat <- mat_stable
    inv_mat[seq_len(p), seq_len(p)] <- sigma_ww_i
    inv_mat <- chol2inv(chol(inv_mat))
    exp_x[i, ] <- mu_w +
      mat_xxxz %*% inv_mat %*% c(w_mat[i, ] - mu_w, z_mat[i, ] - mu_z)
  }

  exp_x
}

#' Negative log-likelihood for MVN
#'
#' @inheritParams .score_mvn
#' @return Negative log-likelihood for MVN
#' @keywords internal
.negll_fiml <- function(theta, x_mat) {
  n <- nrow(x_mat)
  p <- ncol(x_mat)

  mu <- theta[1:p]
  c_chol <- matrix(0, p, p)
  c_chol[lower.tri(c_chol, diag = TRUE)] <- theta[(p + 1):length(theta)]
  diag(c_chol) <- exp(diag(c_chol)) + 1e-6
  c_mat <- tcrossprod(c_chol)

  x_miss <- !is.na(x_mat)
  miss_pats <- apply(x_miss, 1, function(x) {
    paste0(as.integer(x), collapse = "")
  })
  groups <- split(seq_len(n), miss_pats)

  nll <- 0

  for (idx_s in groups) {
    curr_idxs <- which(x_miss[idx_s[1], ])
    if (length(curr_idxs) == 0) next

    x_i <- x_mat[idx_s, curr_idxs, drop = FALSE]
    mu_i <- mu[curr_idxs]
    c_mat_i <- c_mat[curr_idxs, curr_idxs, drop = FALSE]
    diag(c_mat_i) <- diag(c_mat_i) + 1e-8
    c_chol_i <- tryCatch(chol(c_mat_i), error = function(e) NULL)
    if (is.null(c_chol_i)) {
      return(1e10)
    }

    x_centered <- sweep(x_i, 2, mu_i)
    z <- forwardsolve(t(c_chol_i), t(x_centered))
    ln_det <- 2 * sum(log(diag(c_chol_i)))

    nll <- nll + 0.5 * (length(idx_s) * ln_det + sum(z^2))
  }

  nll
}

#' Score function for FIML of MVN
#'
#' @param theta (vector) parameter candidates
#' @param x_mat (matrix) matrix with missing elements
#' @return Gradients of mean and covariance
#' @keywords internal
.score_mvn <- function(theta, x_mat) {
  # https://people.csail.mit.edu/jrennie/writing/multivariateNormal.pdf
  n <- nrow(x_mat)
  p <- ncol(x_mat)

  mu <- theta[1:p]
  c_chol <- matrix(0, p, p)
  c_chol[lower.tri(c_chol, diag = TRUE)] <- theta[(p + 1):length(theta)]
  diag(c_chol) <- exp(diag(c_chol)) + 1e-6
  c_mat <- tcrossprod(c_chol)

  grad_mu <- rep(0, p)
  grad_c_mat <- matrix(0, p, p)

  x_miss <- !is.na(x_mat)
  miss_pats <- apply(x_miss, 1, function(x) {
    paste0(as.integer(x), collapse = "")
  })
  groups <- split(seq_len(n), miss_pats)

  for (idx_s in groups) {
    curr_idxs <- which(x_miss[idx_s[1], ])
    if (length(curr_idxs) == 0) next

    x_i <- x_mat[idx_s, curr_idxs, drop = FALSE]
    mu_i <- mu[curr_idxs]
    c_mat_i <- c_mat[curr_idxs, curr_idxs, drop = FALSE]

    c_inv <- chol2inv(chol(c_mat_i))
    x_centered <- sweep(x_i, 2, mu_i)

    grad_mu_i <- -colSums(x_centered %*% c_inv)
    grad_mu[curr_idxs] <- grad_mu[curr_idxs] + grad_mu_i

    x_t_x <- crossprod(x_centered)
    p1 <- c_inv %*% x_t_x %*% c_inv
    p2 <- length(idx_s) * c_inv

    grad_pat <- -.5 * (p1 - p2)
    grad_c_mat[curr_idxs, curr_idxs] <- grad_c_mat[curr_idxs, curr_idxs] +
      grad_pat
  }

  grad_c_mat <- .5 * (grad_c_mat + t(grad_c_mat))
  grad_l <- 2 * grad_c_mat %*% c_chol
  diag(grad_l) <- diag(grad_l) * diag(c_chol)
  grad_l_pars <- grad_l[lower.tri(grad_l, diag = TRUE)]

  c(grad_mu, grad_l_pars)
}

#' Optimisation function for MVN FIML
#'
#' @inheritParams .score_mvn
#' @return Result from optim
#' @keywords internal
.optim_cov_fiml <- function(x_mat) {
  p <- ncol(x_mat)
  x_mat_s <- scale(x_mat)
  init_mu <- attr(x_mat_s, "scaled:center")
  d_mat <- diag(p)
  diag(d_mat) <- attr(x_mat_s, "scaled:scale")
  init_cov <- stats::cov(x_mat_s, use = "p")
  init_cov[is.na(init_cov)] <- 0
  diag(init_cov) <- diag(init_cov) + 1e-2
  init_chol <- chol(init_cov)
  diag(init_chol) <- log(diag(init_chol))
  init_pars <- c(
    rep(0, p),
    init_chol[lower.tri(init_chol, diag = TRUE)]
  )
  ret <- stats::optim(
    init_pars,
    function(theta) rccme:::negll_fiml_cpp(theta, x_mat_s),
    function(theta) rccme:::score_mvn_cpp(theta, x_mat_s),
    method = "BFGS",
    control = list(maxit = 1e3, pgtol = 1e-6)
  )
  mu <- ret$par[seq_len(p)]
  c_mat <- diag(p)
  c_mat[lower.tri(c_mat, diag = TRUE)] <- ret$par[-seq_len(p)]
  diag(c_mat) <- exp(diag(c_mat))
  c_mat <- tcrossprod(c_mat)
  list(mu = mu + init_mu, c_mat = d_mat %*% c_mat %*% d_mat)
}

#' Cross-covariance when y has missing values
#'
#' @param x_mat (matrix) complete matrix
#' @param y_mat (matrix) matrix with missing elements
#' @param mu_x (vector) known mean vector of x_mat
#' @param mu_y (vector) known mean vector of y_mat
#' @param cov_y (matrix) known covariance of y_mat
#' @return Complete cross-covariance
#' @keywords internal
.cov_xy_fiml <- function(x_mat, y_mat, mu_x, mu_y, cov_y) {
  n <- nrow(x_mat)
  p <- ncol(x_mat)
  q <- ncol(y_mat)

  x_centered <- sweep(x_mat, 2, mu_x)

  y_miss <- !is.na(y_mat)
  miss_pats <- apply(y_miss, 1, function(x) {
    paste0(as.integer(x), collapse = "")
  })
  groups <- split(seq_len(n), miss_pats)

  cov_xy <- matrix(0, p, q)

  for (idx_s in groups) {
    curr_idxs <- which(y_miss[idx_s[1], ])
    miss_idxs <- which(is.na(y_mat[idx_s[1], ]))

    y_i <- y_mat[idx_s, curr_idxs, drop = FALSE]
    x_i <- x_centered[idx_s, , drop = FALSE]

    mu_c <- mu_y[curr_idxs]

    y_centered <- sweep(y_i, 2, mu_c)
    y_comp <- matrix(0, nrow(y_i), q)
    y_comp[, curr_idxs] <- y_centered

    if (length(miss_idxs) > 0) {
      cov_y_c <- cov_y[curr_idxs, curr_idxs, drop = FALSE]
      cov_y_m <- cov_y[miss_idxs, curr_idxs, drop = FALSE]

      c_inv <- chol2inv(chol(cov_y_c))
      y_comp[, miss_idxs] <- t(cov_y_m %*% c_inv %*% t(y_centered))
    }

    # Cross-cov contribution
    cov_xy <- cov_xy + crossprod(x_i, y_comp)
  }

  cov_xy / n
}


