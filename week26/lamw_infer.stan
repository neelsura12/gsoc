functions
{
    int signum(real y)
    {
        return y > 0 ? 1 : -1;
    }
    real lambertw_delta(real y, real delta)
    {
	real x;
    	x = signum(y) * sqrt(lambert_w0(delta * square(y)) / delta);
	return x;	
    }
}
data
{
    int N;
    real y[N];
}
parameters
{
    real<lower=0> delta;
    real mu;
    real<lower=0> sigma;
}
model
{
    delta ~ normal(0, 1); 
    mu ~ normal(0, 1);
    sigma ~ normal(0, 1); 
    for (i in 1:N) 
    {
	real y_hat = (y[i] - mu)/sigma;
	real x = lambertw_delta(y_hat, delta);
        target += normal_lpdf(x*sigma + mu | mu, sigma);
        target += log(x/(y_hat*(1 + delta*square(x))));
    }
}
