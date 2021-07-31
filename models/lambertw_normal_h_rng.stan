data {
    real mu; // offset
    real<lower=0> sigma; // multiplier
    real<lower=0> delta; // symmetric-tail index
}
generated quantities {
    real x = normal_rng(mu, sigma);
    real u = (x - mu)/sigma;
    real z = u * exp(delta/2 * square(u));
    real y = z * sigma + mu;
}
