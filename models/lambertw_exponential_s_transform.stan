data {
    int N;
    vector[N] y;
}
parameters {
    real<lower=0> lambda;
    real<lower=0> gamma;
}
model {
    vector[N] z = lambda * y;
    vector[N] w_gamma_z = lambert_w0(gamma * z);
    
    lambda ~ normal(0, sqrt(pi()/2));
    gamma ~ exponential(1);

    (w_gamma_z/(gamma * lambda)) ~ exponential(lambda);
    target += log(w_gamma_z) - log(gamma * z) - log1p(w_gamma_z);
}