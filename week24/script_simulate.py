import numpy as np
from eth_tools.utils.stan_utility import compile_model
from matplotlib import pyplot as plt

model = compile_model('/week24/lamw_simu.stan')
model.sampling(data=dict(N=100, m=0, s=1, d=1), chains=1, iter=100, algorithm='Fixed_param')
out = model.sampling(data=dict(N=100, m=0, s=1, d=1), chains=1, iter=100, algorithm='Fixed_param')

fig, axes = plt.subplots(2, 1)
axes[0].hist(out.extract('y')['y'].flatten(), bins=np.arange(-10, 10, 0.5))
axes[1].hist(out.extract('z')['z'].flatten(), bins=np.arange(-10, 10, 0.5))
plt.show()
