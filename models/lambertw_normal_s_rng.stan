data {
    int N;
    real mu; // offset
    real<lower=0> sigma; // multiplier
    real gamma; // skewness parameter
}
generated quantities {
    real x;
    real y;
    real u; 
    real z;
    for ( i in 1:N ) {
        x = normal_rng(mu, sigma);
        u = (x - mu)/sigma;
        z = u * exp(gamma * u);
        y = z*sigma + mu;
    }
}
