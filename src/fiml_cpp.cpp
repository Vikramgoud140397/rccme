// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ------------------------------------------------------------
// Helper: build Sigma from the packed lower-Cholesky parameter
// vector (same parameterisation as the R code).
//
// theta = [mu_1 ... mu_p | L_11 L_21 L_22 L_31 ... L_pp]
// diagonal elements of L are stored on the log scale, so the
// true diagonal is exp(L_ii).
// ------------------------------------------------------------
static mat chol_to_sigma(const vec& theta, int p) {
  mat L(p, p, fill::zeros);
  // fill lower triangle from theta[p .. end]
  int idx = p;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      L(i, j) = theta(idx++);
    }
    L(j, j) = std::exp(L(j, j)) + 1e-6;   // back-transform diagonal
  }
  return L * L.t();   // Sigma = L L'
}

// same thing but also return L itself (needed for the gradient)
static mat chol_to_L(const vec& theta, int p) {
  mat L(p, p, fill::zeros);
  int idx = p;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      L(i, j) = theta(idx++);
    }
    L(j, j) = std::exp(L(j, j)) + 1e-6;
  }
  return L;
}

// ------------------------------------------------------------
// Build the missing-data pattern table once and reuse it for
// both the NLL and the gradient.  Returns a list of index
// vectors (0-based) — one entry per unique pattern.
// ------------------------------------------------------------
static std::vector<uvec> make_groups(const mat& X) {
  int n = X.n_rows;
  // encode each row as a string "0101..."
  std::map<std::string, std::vector<int>> grp_map;
  for (int i = 0; i < n; ++i) {
    std::string key(X.n_cols, '0');
    for (int j = 0; j < (int)X.n_cols; ++j)
      if (!std::isnan(X(i, j))) key[j] = '1';
    grp_map[key].push_back(i);
  }
  std::vector<uvec> groups;
  groups.reserve(grp_map.size());
  for (auto& kv : grp_map) {
    uvec idx(kv.second.size());
    for (int i = 0; i < (int)kv.second.size(); ++i)
      idx(i) = kv.second[i];
    groups.push_back(idx);
  }
  return groups;
}

// observed column indices for the first row of a group
static uvec obs_cols(const mat& X, int row) {
  std::vector<uword> cols;
  for (int j = 0; j < (int)X.n_cols; ++j)
    if (!std::isnan(X(row, j))) cols.push_back(j);
  uvec out(cols.size());
  for (int i = 0; i < (int)cols.size(); ++i) out(i) = cols[i];
  return out;
}

// ------------------------------------------------------------
// Negative log-likelihood (FIML, MVN)
// Exported so optim() can call it from R.
// ------------------------------------------------------------
// [[Rcpp::export]]
double negll_fiml_cpp(const arma::vec& theta, const arma::mat& X) {
  int n = X.n_rows;
  int p = X.n_cols;

  vec mu  = theta.head(p);
  mat Sigma = chol_to_sigma(theta, p);

  auto groups = make_groups(X);

  double nll = 0.0;

  for (auto& idx : groups) {
    uvec cols = obs_cols(X, idx(0));
    if (cols.is_empty()) continue;

    mat  X_i     = X.submat(idx, cols);           // n_i x q_i
    vec  mu_i    = mu(cols);
    mat  Sigma_i = Sigma.submat(cols, cols);
    Sigma_i.diag() += 1e-8;

    // Cholesky of Sigma_i
    mat  C_i;
    bool ok = chol(C_i, Sigma_i, "upper");        // C_i' * C_i = Sigma_i
    if (!ok) return 1e10;

    int  ni = idx.n_elem;
    int  qi = cols.n_elem;

    // centred data
    mat  Xc = X_i.each_row() - mu_i.t();          // n_i x q_i

    // solve C_i' z = Xc'  =>  z = (C_i')^{-1} Xc'
    mat  Z = solve(trimatu(C_i).t(), Xc.t());     // qi x ni

    double ln_det = 2.0 * sum(log(C_i.diag()));

    nll += 0.5 * (ni * ln_det + accu(Z % Z));
  }

  return nll;
}

// ------------------------------------------------------------
// Gradient of the negative log-likelihood w.r.t. theta
// (same parameterisation).
// ------------------------------------------------------------
// [[Rcpp::export]]
arma::vec score_mvn_cpp(const arma::vec& theta, const arma::mat& X) {
  int n = X.n_rows;
  int p = X.n_cols;

  vec mu = theta.head(p);
  mat L  = chol_to_L(theta, p);
  mat Sigma = L * L.t();

  vec  grad_mu(p,    fill::zeros);
  mat  grad_Sigma(p, p, fill::zeros);

  auto groups = make_groups(X);

  for (auto& idx : groups) {
    uvec cols = obs_cols(X, idx(0));
    if (cols.is_empty()) continue;

    mat  X_i      = X.submat(idx, cols);
    vec  mu_i     = mu(cols);
    mat  Sigma_i  = Sigma.submat(cols, cols);

    mat  C_inv = inv_sympd(Sigma_i);

    mat  Xc = X_i.each_row() - mu_i.t();          // n_i x q_i
    int  ni = idx.n_elem;

    // gradient w.r.t. mu (observed subset)
    vec gmu_i = -C_inv * sum(Xc, 0).t();           // q_i-vector
    for (int k = 0; k < (int)cols.n_elem; ++k)
      grad_mu(cols(k)) += gmu_i(k);

    // gradient w.r.t. Sigma (observed subset)
    mat XtX   = Xc.t() * Xc;                       // q_i x q_i
    mat gS_i  = -0.5 * (C_inv * XtX * C_inv - ni * C_inv);
    for (int r = 0; r < (int)cols.n_elem; ++r)
      for (int c2 = 0; c2 < (int)cols.n_elem; ++c2)
        grad_Sigma(cols(r), cols(c2)) += gS_i(r, c2);
  }

  // symmetrise
  grad_Sigma = 0.5 * (grad_Sigma + grad_Sigma.t());

  // chain rule: Sigma = L L',  dNLL/dL = 2 * (dNLL/dSigma) * L
  mat grad_L = 2.0 * grad_Sigma * L;
  // diagonal lives on log scale: multiply by exp(l_ii) = L_ii (before +1e-6)
  // The stored diagonal value gave  L_ii = exp(theta_ii)+1e-6,
  // so d/d(theta_ii) = L_ii * exp(theta_ii) / (exp(theta_ii)+1e-6)
  // ≈ L_ii for large values; exact chain rule:
  for (int j = 0; j < p; ++j)
    grad_L(j, j) *= (L(j, j) - 1e-6);   // exp(theta_jj) = L_jj - 1e-6

  // pack lower triangle (column-major, same order as R)
  int n_pars = p + p * (p + 1) / 2;
  vec grad(n_pars, fill::zeros);
  grad.head(p) = grad_mu;
  int idx2 = p;
  for (int j = 0; j < p; ++j)
    for (int i = j; i < p; ++i)
      grad(idx2++) = grad_L(i, j);

  return grad;
}
