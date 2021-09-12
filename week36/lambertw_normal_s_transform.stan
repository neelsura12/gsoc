data {
    int N;
    vector[N] y;
    real mu;
    real<lower=0> sigma;
}
parameters {
    real<lower=0> gamma;
}
model {
    vector[N] z = (y - mu) / sigma;
    vector[N] w_gamma_z_pr = lambert_w0(gamma * z);
    vector[N] w_gamma_z_nonpr = lambert_wm1(gamma * z);
    vector[N] w_gamma_z_pr_normal = (w_gamma_z_pr  / gamma) * sigma + mu;
    vector[N] w_gamma_z_nonpr_normal = (w_gamma_z_nonpr  / gamma) * sigma + mu;

    gamma ~ exponential(1);
    
    real lower_bound = -sigma/(gamma * e()) + mu;
    for (i in 1:N) {
        if (y[i] >= lower_bound) {
            w_gamma_z_nonpr_normal[i] ~ normal(mu, sigma);
            target += log(w_gamma_z_pr[i]) - log(gamma * z[i]) - log1p(w_gamma_z_pr[i]);
            if (y[i] <= mu) {
                target += N * log(sigma) + 0.5 / square(gamma) * w_gamma_z_nonpr[i];
                target += -log(w_gamma_z_nonpr[i]) + log(gamma * z[i]) + log1p(w_gamma_z_nonpr[i]);
            }
        }
    }
}