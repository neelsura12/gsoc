functions
{
    int signum(real y)
    {
        return y > 0 ? 1 : -1;
    }
}
data
{
    int N;
    real y[N];
}
parameters
{
    real<lower=0> d;
}
transformed parameters
{
    real z;
    // z = sign(y) * lambert_w0(d * y^2) / d; // doesn't work as lambert_w0 not found
}
model
{
    d ~ gamma(1, 2);
    z ~ normal(0, 1);
    // todo jacobian adjustment since i'm putting a prior on a transformed parameter
}