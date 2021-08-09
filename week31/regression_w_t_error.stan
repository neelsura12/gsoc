data {
    int N;
    vector[N] y;
    vector[N] x;
}
parameters {
    real<lower=2> nu;
    real beta;
    real alpha;
}
transformed parameters {
    real<lower=0> sigma = sqrt(nu/(nu - 2));
    vector[N] mu = alpha + beta * x;
}
model {
    sigma ~ normal(0, sqrt(pi()/2));
    nu ~ exponential(1);
    alpha ~ normal(0.5, 1);
    beta ~ normal(3, 1);
    y ~ student_t(nu, mu);
}
generated quantities {
    real x_new = poisson_rng(3);
    real mu_new = alpha + beta * x_new;
    real y_new = student_t_rng(mu_new, sigma);
}