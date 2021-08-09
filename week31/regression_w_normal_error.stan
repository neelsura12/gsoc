data {
    int N;
    vector[N] y;
    vector[N] x;
}
parameters {
    real<lower=0> sigma;
    real beta;
    real alpha;
}
transformed parameters {
    vector[N] mu = alpha + beta * x;
}
model {
    sigma ~ normal(0, sqrt(pi()/2));
    alpha ~ normal(0.5, 1);
    beta ~ normal(3, 1);
    y ~ normal(mu, sigma);
}
generated quantities {
    real x_new = poisson_rng(3);
    real mu_new = alpha + beta * x_new;
    real y_new = normal_rng(mu_new, sigma);
}