SABRE Methods 
---------
Understanding how viruses offer protection against closely related emerging strains is vital for creating effective vaccines. 
For many viruses, multiple serotypes often co-circulate and testing large numbers of vaccines can be infeasible. 
Therefore the development of an in silico predictor of cross-protection between strains is important to help optimise vaccine choice. 
Here we present a sparse hierarchical Bayesian model for detecting relevant antigenic sites in virus evolution (SABRE) 
which can account for the experimental variability in the data and predict antigenic variability. 
The method uses spike and slab priors to identify sites in the viral protein which are important for the neutralisation of the virus.
Additionally, we propose an extended sparse hierarchical Bayesian model that can deal with the pairwise structure and inherent experimental variability
in the antigenic variability data through the introduction of latent variables (eSABRE).

Content 
---------

This repository contains the original code (R) for the SABRE and eSABRE methods, which are available in the folder 
['Original Methods'](https://github.com/vinnydavies/sabre_methods/original_methods).
Additionally, there is ongoing work looking at improving the speed and performance of the SABRE methods available in the folder 
['Updated Methods'](https://github.com/vinnydavies/sabre_methods/updated_methods).


Examples 
---------

Examples can be found in the folder ['Example'](https://github.com/vinnydavies/sabre_methods/examples)

Publications 
---------

To reference the original SABRE methods, please cite our papers:

Davies, V., Reeve, R., Harvey, W., Maree, F., & Husmeier, D. (2014, April). [Sparse Bayesian variable selection for the identification of antigenic variability
in the foot-and-mouth disease virus.](http://proceedings.mlr.press/v33/davies14.pdf) In Artificial Intelligence and Statistics (pp. 149-158). 

Davies, V., Reeve, R., Harvey, W. T., Maree, F. F., & Husmeier, D. (2017). 
[A sparse hierarchical Bayesian model for detecting relevant antigenic sites in virus evolution.](https://link.springer.com/article/10.1007/s00180-017-0730-6)
Computational Statistics, 32(3), 803-843.

To reference the eSABRE method, please also cite the following paper:

Davies, V., Harvey, W. T., Reeve, R., & Husmeier, D. (2019). 
[Improving the identification of antigenic sites in the H1N1 influenza virus through accounting for the experimental structure in a sparse hierarchical Bayesian model.](
https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssc.12338?af=R)
Journal of the Royal Statistical Society: Series C (Applied Statistics), 68(4), 859-885.
