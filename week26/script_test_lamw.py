from cmdstanpy import CmdStanModel

N = 1000
simu_model = CmdStanModel(stan_file='/Users/nshah/work/gsoc/week26/lamw_simu.stan')
simu_fit = simu_model.sample(data={'mu': 0.5, 'sigma': 2, 'delta': 1}, chains=1, iter_sampling=N, iter_warmup=0, fixed_param=True);
y = simu_fit.draws_pd()['y'].values

infer_model = CmdStanModel(stan_file='/Users/nshah/work/gsoc/week26/lamw_infer.stan')
infer_fit = infer_model.sample(data={'y': y, 'N': N}, chains=4);
infer_fit.diagnose()
infer_fit.summary()
