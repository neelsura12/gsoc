data
{
    real<lower=0> delta;
    real mu;
    real<lower=0> sigma;
}
generated quantities
{
    real x = normal_rng(mu, sigma);
    real y = x * exp(delta/2*x^2);
}
