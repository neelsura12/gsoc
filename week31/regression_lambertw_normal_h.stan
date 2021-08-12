data {
    int N;
    vector[N] y;
    vector[N] x;
}
parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
    real<lower=0> delta;
}
transformed parameters {
    vector[N] mu = alpha + beta * x;
}
model {
    vector[N] z = (y - mu) / sigma;
    vector[N] w_delta_z_sq = lambert_w0(delta * square(z));
    
    alpha ~ normal(0.5, 1);
    beta ~ normal(3, 1);
    sigma ~ normal(0, sqrt(pi()/2));
    delta ~ exponential(1);

    target += -N * log(sigma) - 0.5 * square(z)' * exp(-w_delta_z_sq);
    target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta) - log(fabs(z)) - log1p(w_delta_z_sq);
}
generated quantities {
    // new datapoint of (X, Y)
    real new_x = poisson_rng(3);
    real mu_hat = alpha + beta * new_x;
    real u = normal_rng(0, 1);
    real new_y = u * exp(delta/2 * square(u))*sigma + mu_hat;
}