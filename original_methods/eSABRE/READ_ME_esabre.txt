This file contains code for the eSABRE method from \cite{Davies18}.
The method is contained in SC_SABRE.R and that will call the the other functions, 

@article{Davies18,
	author = {Vinny Davies and William T. Harvey and Richard Reeve and Dirk Husmeier},
	title = {Improving the identification of antigenic sites in the {H1N1} Influenza virus through accounting for the experimental structure in a sparse hierarchical {B}ayesian model},
	journal = {In Submission},
	year = {2018},
	volume = { },
	number = { },
	pages = { },
}

Inputs into the model are as follows

y - Response vector
X - Explanatory variables
Z - List of random effect vectors (one vector per group)
Challenge - vector of challenge strains
Protective - vector of challenge strains
iter - Number of iterations
it_count - Method prints iterations every it_count MCMC samples
alpha_e, beta_e - hyper-parameter for \simga_\varepsilon^2 (sigsq_e) 
alpha_e, beta_e - hyper-parameter for \simga_y^2 (sigsq_y) 
alpha_w, beta_w - hyper-parameter for \simga_w^2 (sigsq_w) 
alpha_b, beta_b - hyper-parameter for \simga_b^2 (sigsq_b) 
pi_a, pi_b - hyper-parameter for \pi (pi)
gamma0, w0, sigsq_e0, sigsq_w0, b0, sigsq_b0, mu_w0, pi0 - are parameter initialisations 
f - the number of \gamma parameters to propose simultaneously
pi_prop - MCMC tuning parameter. Probabiliy of proposing \gamma=1

Model returns

t - vector of total times at each iteration
Zname - random effect factor levels (not checked, advise not to use)
X - explanatory variables for mu_y
y - response variable
M - design matrix for virus pair specification
m_prot - protective strains
m_chal - challenge strains
levelcolZ - random effect matrix summary statistics
Z - random effects matrix
Mvec - internal parameter for eSABRE virus pair specification
wf - intercept parameter and regression coefficients (where 0 if \gamma=0) (\textbf{w}^*)
sigsq_e - latent variable fit error variance (\sigma_\varepsilon^2)
sigsq_w - regression coefficients error variance (\sigma_w^2)
b - random-effect coefficients
sigsq_b - random-effect coefficients error variance (\sigma_b^2)
gamma - binary selection parameters (\boldsymbol\gamma)
pi - probability of regression coefficient inclusion parameter (\pi)
mu_w - mean of regression coefficients parameter (\mu_w)
mu_y - samples of estimated true HI assay values (\mu_y)
sigsq_y - model fit error variable (\sigma_y^2)