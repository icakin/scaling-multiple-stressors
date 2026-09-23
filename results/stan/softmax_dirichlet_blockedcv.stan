
functions {
  vector segment_vector(vector x, int start, int L) {
    vector[L] y;
    for (m in 1:L) y[m] = x[start + m - 1];
    return y;
  }
}
data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> J;
  int<lower=1> S;
  int<lower=1> start_idx[N];
  int<lower=1> len[N];
  vector[J] p_obs;
  vector[J] g_z;
  int<lower=1, upper=K> tax_of[J];
  int<lower=1, upper=S> stress_id[N];
  int<lower=1> Rdim;
  int<lower=1, upper=Rdim> rich_id[N];
  int<lower=0, upper=1> use_taxon_biases;
}
parameters {
  real mu_kappa;
  real<lower=0> sigma_kappa;
  matrix[S, Rdim] kappa_raw;
  vector[K-1] delta_raw;
  real<lower=0> sigma_delta;
  real mu_phi;
  real<lower=0> sigma_phi;
  vector[S] log_phi_raw;
}
transformed parameters {
  matrix[S, Rdim] kappa;
  vector[K] delta;
  vector[S] log_phi;
  vector[S] phi;

  kappa = mu_kappa + sigma_kappa * kappa_raw;

  delta[1:(K-1)] = use_taxon_biases * sigma_delta * delta_raw;
  delta[K]       = 0;

  log_phi = mu_phi + sigma_phi * log_phi_raw;
  for (s in 1:S) phi[s] = exp(log_phi[s]);
}
model {
  mu_kappa    ~ normal(0, 1);
  sigma_kappa ~ normal(0, 0.5);
  to_vector(kappa_raw) ~ normal(0, 1);

  sigma_delta ~ normal(0, 1);
  delta_raw   ~ normal(0, 1);

  mu_phi    ~ normal(log(50), 1);
  sigma_phi ~ normal(0, 1);
  log_phi_raw ~ normal(0, 1);

  for (n in 1:N) {
    int a = start_idx[n];
    int L = len[n];
    vector[L] eta;
    vector[L] pseg = segment_vector(p_obs, a, L);

    pseg = pseg + 1e-12;
    pseg = pseg / sum(pseg);

    for (m in 1:L) {
      int j = a + m - 1;
      int t = tax_of[j];
      eta[m] = kappa[stress_id[n], rich_id[n]] * g_z[j] + delta[t];
    }

    target += dirichlet_lpdf(pseg | phi[stress_id[n]] * softmax(eta));
  }
}


