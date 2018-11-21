This file contains code for the Semi-Conjugate SABRE method from \cite{Davies17_COST}.
The method is contained in SC_SABRE.R and that will call the the other functions, 

@article{Davies17_COST,
	author = {Vinny Davies and Richard Reeve and William T. Harvey and Francois F Maree and Dirk Husmeier},
	title = {A sparse hierarchical {B}ayesian model for detecting relevant antigenic sites in virus evolution},
	journal = {Computational Statistics},
	year = {2017},
	volume = {32},
	number = {3},
	pages = {803--843},
}

Inputs into the model are as follows

y - Response vector
X - Explanatory variables
Z - List of random effect vectors (one vector per group)
iter - Number of iterations
it_count - Method prints iterations every it_count MCMC samples
mu_hyp, var_mu_w - hyper-parameters for \mu_w (mu_w)
alpha_e, beta_e - hyper-parameter for \simga_\varepsilon^2 (sigsq_e) 
alpha_w, beta_w - hyper-parameter for \simga_w^2 (sigsq_w) 
alpha_b, beta_b - hyper-parameter for \simga_b^2 (sigsq_b) 
pi_a, pi_b - hyper-parameter for \pi (pi)
gamma0, w0, sigsq_e0, sigsq_w0, b0, sigsq_b0, mu_w0, pi0 - are parameter initialisations 
f - the number of \gamma parameters to propose simultaneously
pi_prop - MCMC tuning parameter. Probabiliy of proposing \gamma=1
prop_mu_mu, prop_var_mu, prop_mu_sig, prop_var_sig - MCMC tuning parameters need when all \gamma=0
sigsq_w_change - allows you to change the group structure of regression coefficients (not described in paper and not tested - stick to default advised)

Model returns

t - vector of total times at each iteration
wf - intercept parameter and regression coefficients (where 0 if \gamma=0) (\textbf{w}^*)
sigsq_e - model error variance (\sigma_\varepsilon^2)
sigsq_w - regression coefficients error variance (\sigma_w^2)
b - random-effect coefficients
sigsq_b - random-effect coefficients error variance (\sigma_b^2)
gamma - binary selection parameters (\boldsymbol\gamma)
pi - probability of regression coefficient inclusion parameter (\pi)
iters - iterations
mu_w - mean of regression coefficients parameter (\mu_w)
accept - gamma proposal acceptance rate
accept2 - adjusted gamma proposal acceptance rate (correct for proposals that are the same as current state