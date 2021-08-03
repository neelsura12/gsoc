data {
    int N;
    vector[N] y;
    vector[N] s;
}
parameters {
    real<lower=0> sigma;
    real<lower=0> delta;
    real beta;
    real alpha;
}
transformed parameters
{
    vector[N] mu = alpha + beta * s;
}
model {
    vector[N] z = (y - mu) / sigma;
    vector[N] w_delta_z_sq = lambert_w0(delta * square(z));
    
    sigma ~ normal(0, sqrt(pi()/2));
    delta ~ exponential(1);
    alpha ~ normal(0.5, 1);
    beta ~ normal(3, 1);

    target += -N * log(sigma) - 0.5 * square(z)' * exp(-w_delta_z_sq);
    target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta) - log(fabs(z)) -log1p(w_delta_z_sq);
}
generated quantities {
    // new datapoint of covariate s
    real new_s = poisson_rng(3);
    real mu_hat = new_s * beta + alpha;
    real new_x = normal_rng(mu_hat, sigma);
    real u = (new_x-mu_hat)/sigma;
    real new_y = u * exp(delta/2 * square(u))*sigma + mu_hat;
}