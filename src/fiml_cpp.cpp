// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// theta = [mu_1 ... mu_p | L_11 L_21 L_22 L_31 ... L_pp]
// diagonal of L is on log scale
static mat chol_to_sigma(const vec& theta, int p) {
  mat L(p, p, fill::zeros);
  int idx = p;
  for (int j = 0; j < p; ++j) {
    for (int i = j; i < p; ++i) {
      L(i, j) = theta(idx++);
    }
    L(j, j) = std::exp(L(j, j)) + 1e-6;
  }
  return L * L.t();
}

// returns L rather than L L'
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

static std::vector<uvec> make_groups(const mat& X) {
  int n = X.n_rows;
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

static uvec obs_cols(const mat& X, int row) {
  std::vector<uword> cols;
  for (int j = 0; j < (int)X.n_cols; ++j)
    if (!std::isnan(X(row, j))) cols.push_back(j);
  uvec out(cols.size());
  for (int i = 0; i < (int)cols.size(); ++i) out(i) = cols[i];
  return out;
}

// [[Rcpp::export]]
double negll_fiml_cpp(const arma::vec& theta, const arma::mat& X) {
  int p = X.n_cols;
  vec mu = theta.head(p);
  mat sigma = chol_to_sigma(theta, p);
  auto groups = make_groups(X);
  double nll = 0.0;
  for (auto& idx : groups) {
    uvec cols = obs_cols(X, idx(0));
    if (cols.is_empty()) continue;
    mat x_i = X.submat(idx, cols);
    vec mu_i = mu(cols);
    mat sigma_i = sigma.submat(cols, cols);
    sigma_i.diag() += 1e-8;
    mat c_i;
    bool ok = chol(c_i, sigma_i, "upper");
    if (!ok) return 1e10;
    int ni = idx.n_elem;
    mat xc = x_i.each_row() - mu_i.t();
    mat z = solve(trimatu(c_i).t(), xc.t());
    double ln_det = 2.0 * sum(log(c_i.diag()));
    nll += 0.5 * (ni * ln_det + accu(z % z));
  }
  return nll;
}

// [[Rcpp::export]]
arma::vec score_mvn_cpp(const arma::vec& theta, const arma::mat& X) {
  int p = X.n_cols;
  vec mu = theta.head(p);
  mat L = chol_to_L(theta, p);
  mat sigma = L * L.t();
  vec grad_mu(p, fill::zeros);
  mat grad_sigma(p, p, fill::zeros);
  auto groups = make_groups(X);
  for (auto& idx : groups) {
    uvec cols = obs_cols(X, idx(0));
    if (cols.is_empty()) continue;
    mat x_i = X.submat(idx, cols);
    vec mu_i = mu(cols);
    mat sigma_i = sigma.submat(cols, cols);
    mat c_inv = inv_sympd(sigma_i);
    mat xc = x_i.each_row() - mu_i.t();
    int ni = idx.n_elem;
    // gradient w.r.t. mu
    vec gmu_i = -c_inv * sum(xc, 0).t();
    for (int k = 0; k < (int)cols.n_elem; ++k)
      grad_mu(cols(k)) += gmu_i(k);
    // gradient w.r.t. sigma
    mat xtx = xc.t() * xc;
    mat gs_i = -0.5 * (c_inv * xtx * c_inv - ni * c_inv);
    for (int r = 0; r < (int)cols.n_elem; ++r)
      for (int c2 = 0; c2 < (int)cols.n_elem; ++c2)
        grad_sigma(cols(r), cols(c2)) += gs_i(r, c2);
  }
  // symmetrise
  grad_sigma = 0.5 * (grad_sigma + grad_sigma.t());
  // chain rule: sigma = L L'
  mat grad_l = 2.0 * grad_sigma * L;
  for (int j = 0; j < p; ++j)
    grad_l(j, j) *= (L(j, j) - 1e-6);
  // pack lower triangle, column-major
  int n_pars = p + p * (p + 1) / 2;
  vec grad(n_pars, fill::zeros);
  grad.head(p) = grad_mu;
  int idx2 = p;
  for (int j = 0; j < p; ++j)
    for (int i = j; i < p; ++i)
      grad(idx2++) = grad_l(i, j);
  return grad;
}
