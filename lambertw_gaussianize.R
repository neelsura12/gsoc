library(LambertW)

# simu data y = x^2+e with iid cauchy errors
x = seq(0, 10, 0.1)
x_sq = sapply(x, function(x) x^2)
e = rcauchy(n=length(x))
y = x_sq + e

# gaussianize the error (y - x_sq)
tmp = Gaussianize(y - x_sq, return.tau.mat=TRUE)
e_hat = tmp$input
tau_mat = tmp$tau.mat # to ungaussianize

# model assumes normal error, but reality is cauchy, expect bad fit
model1 = lm(y ~ x_sq)
summary(model1)

# model2, expect closer match to reality as error normalized
model2 = lm((e_hat + x_sq) ~ x_sq)
summary(model2)

# NS: this part shows that MSE is about the same in either case 

# calculate residual
y_pred1 = coef(model1)[1] + coef(model1)[2]*x_sq

# can undo the Gaussianization
e_pred = coef(model2)[1] + coef(model2)[2]*x_sq - x_sq
y_pred2 = x_sq + Gaussianize(e_pred, inverse = TRUE, tau.mat=tau_mat)

print(mean((y_pred1-y)^2))
print(mean((y_pred2-y)^2))