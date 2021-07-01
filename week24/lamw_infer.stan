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
    y = z/2 exp(d * z^2)
}
model
{
    d ~ gamma(1, 2);
    z ~ normal(0, 1);
}