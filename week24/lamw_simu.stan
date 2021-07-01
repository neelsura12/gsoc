data
{
    int N;
    real<lower=0> d;
    real m;
    real<lower=0> s;

}
generated quantities
{
    real y[N];
    real z[N];
    for (i in 1:N)
    {
        z[i] = normal_rng(m, s);
        y[i] = z[i] * exp(d/2*z[i]^2);
    }
}