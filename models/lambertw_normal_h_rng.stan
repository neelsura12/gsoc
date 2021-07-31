data {
    int N;
    real mu; // offset
    real<lower=0> sigma; // multiplier
    real<lower=0> delta; // symmetric-tail index
}
generated quantities {
    real x;
    real y;
    real u; 
    real z;
    for ( i in 1:N ) {
        x = normal_rng(mu, sigma);
        u = (x - mu)/sigma;
        z = u * exp(delta/2 * square(u));
        y = z*sigma + mu;
    }
}
