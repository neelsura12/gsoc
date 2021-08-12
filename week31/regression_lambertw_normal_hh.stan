data {
    int N;
    vector[N] y;
    vector[N] x;
}
parameters {
    real alpha;
    real beta;
    real<lower=0> sigma;
    real<lower=0> delta_left;
    real<lower=0> delta_right;
}
transformed parameters {
    vector[N] mu = alpha + beta * x;
}
model {
    alpha ~ normal(0.5, 1);
    beta ~ normal(3, 1);
    sigma ~ normal(0, sqrt(pi()/2));
    delta_left ~ exponential(1);
    delta_right ~ exponential(1);   
    
    real z, w_delta_z_sq;
    for (i in 1:N) {
        z = (y[i] - mu[i]) / sigma;
        if (z <= 0) {
            w_delta_z_sq = lambert_w0(delta_left * square(z));
            sigma * sqrt(w_delta_z_sq / delta_left) ~ normal(0, sigma);
            // jacobian adjustment for transformed-normal by func. u exp(d/2 u^2)
            target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta_left) - log(fabs(z)) - log1p(w_delta_z_sq);    
        } else {
            w_delta_z_sq = lambert_w0(delta_right * square(z));
            sigma * sqrt(w_delta_z_sq / delta_right) ~ normal(0, sigma);
            // jacobian adjustment for transformed-normal by func. u exp(d/2 u^2)
            target += 0.5 * log(w_delta_z_sq) - 0.5 * log(delta_right) - log(fabs(z)) - log1p(w_delta_z_sq);    
        }
    }
}
generated quantities {
    // new datapoints of (X, Y)
    real new_x = poisson_rng(3);
    real mu_hat = new_x * beta + alpha;
    real u = normal_rng(0, 1);
    real delta = u <= 0 ? delta_left : delta_right;
    real new_y = u * exp(delta/2 * square(u))*sigma + mu_hat;
}