data {
    real mu; // offset
    real<lower=0> sigma; // multiplier
    real<lower=0> delta_left; // left-tail parameter
    real<lower=0> delta_right; // right-tail parameter
}
generated quantities {
    real x = normal_rng(mu, sigma);
    real u = (x - mu)/sigma;
    real z;
    if (u <= 0)
        z = u * exp(delta_left/2 * square(u));
    else
        z = u * exp(delta_right/2 * square(u));
    real y = z * sigma + mu;
}
