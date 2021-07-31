data {
    int N;
    real mu; // offset
    real<lower=0> sigma; // multiplier
    real<lower=0> delta_left; // left-tail index
    real<lower=0> delta_right; // right-tail index
}
generated quantities {
    real x;
    real y;
    real u; 
    real z;
    for ( i in 1:N ) {
        x = normal_rng(mu, sigma);
        u = (x - mu)/sigma;
        if (u <= 0)
            z = u * exp(delta_left/2 * square(u));
        else
            z = u * exp(delta_right/2 * square(u));
        y = z*sigma + mu;
    }
}
