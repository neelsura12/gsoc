functions {
    int signum(real y) {
        return y > 0 ? 1 : -1;
    }
    
    vector standardize(vector y, real mu, real sigma) {
        return (y - mu) / sigma;
    }
    
    vector unstandardize(vector y_std, real mu, real sigma) {
        return y_std * sigma + mu;
    }

    vector lambertw_delta_lp(vector y_std, real delta) {
        int N = num_elements(y_std);
        vector[N] x;

        for (n in 1:N)
            x[n] = signum(y_std[n]) * sqrt(lambert_w0(delta * square(y_std[n])) / delta);
        // target += log(fabs(x)) - log(fabs(y_std)) - log(1 + delta * square(x));
        target += -delta/2 * square(x) - log(1 + delta * square(x));
        return x;
    }
}
data {
    int N;
    vector[N] y;
    real<lower=2.5> nu;
    real mu;
    real<lower=0> sigma;
}
parameters {
    real<lower=0.01> delta;
}
model {
    vector[N] x = unstandardize(lambertw_delta_lp(standardize(y, mu, sigma), delta), mu, sigma);

    delta ~ exponential(1);
    
    x ~ student_t(nu, mu, sigma);
    // target += -log(sigma);
}