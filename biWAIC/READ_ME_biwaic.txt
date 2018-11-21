This file contains the code for nWAIC and biWAIC applied to the SABRE methods.
Code does not currently generalised to other method, but to do this would be trivial.
Code is used via iWAIC_eSABRE function, which then calls the other functions.

@article{Davies18,
	author = {Vinny Davies and William T. Harvey and Richard Reeve and Dirk Husmeier},
	title = {Improving the identification of antigenic sites in the {H1N1} Influenza virus through accounting for the experimental structure in a sparse hierarchical {B}ayesian model},
	journal = {In Submission},
	year = {2018},
	volume = { },
	number = { },
	pages = { },
}

Code input ares

model - eSABRE model list
it_seq - iterations on which we wish to calculate nWAIC / biWAIC

Outputs are

biWAIC - biWAIC value
nWAIC - nWAIC value
biLIKE - (marginal) likelihoods on which biWAIC was evaluated
nLIKE - likelihoods on which nWAIC was evaluated