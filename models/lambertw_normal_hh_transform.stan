data {
    int N;
    vector[N] y;
}
parameters {
    real mu;
    real<lower=0> sigma;
    real<lower=0> delta_left;
    real<lower=0> delta_right;
}
model {
    mu ~ normal(0, 1);
    sigma ~ normal(0, sqrt(pi()/2));
    delta_left ~ exponential(1);
    delta_right ~ exponential(1);
    
    real z, w_delta_z_sq;
    for (i in 1:N) {
        z = (y[i] - mu) / sigma;
        if (z <= 0) {
            w_delta_z_sq = lambert_w0(delta_left * square(z));
            sigma * sqrt(w_delta_z_sq / delta_left) ~ normal(0, sigma);
            target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta_left) - log(fabs(z)) - log1p(w_delta_z_sq);    
        } else {
            w_delta_z_sq = lambert_w0(delta_right * square(z));
            sigma * sqrt(w_delta_z_sq / delta_right) ~ normal(0, sigma);
            target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta_right) - log(fabs(z)) - log1p(w_delta_z_sq);    
        }
    }
}