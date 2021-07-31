data {
    real mu; // offset
    real<lower=0> sigma; // multiplier
    real gamma; // skewness parameter
}
generated quantities {
    real x = normal_rng(mu, sigma);;
    real u = (x - mu)/sigma;
    real z = u * exp(gamma * u);
    real y = z * sigma + mu;
}
