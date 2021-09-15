functions {
    real lambertw_normal_h_transform_lpdf(vector y, real delta, real mu, real sigma)
    {
        int N = num_elements(y);
        vector[N] z = (y - mu) / sigma;
        vector[N] w_delta_z_sq = lambert_w0(delta * square(z));
        return (-N * log(sigma) - 0.5 * square(z)' * exp(-w_delta_z_sq)
        + sum(0.5 * log(w_delta_z_sq) - 0.5 * log(delta) - log(fabs(z)) - log1p(w_delta_z_sq)));
    }
}
data {
    int N;
    vector[N] y;
}
parameters {
    real mu;
    real<lower=0> sigma;
    real<lower=0> delta;
}
model {
    mu ~ normal(0, 1);
    sigma ~ normal(0, sqrt(pi()/2));
    delta ~ exponential(1);
    y ~ lambertw_normal_h_transform(delta, mu, sigma);
}