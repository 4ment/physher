Here are some examples of configuration files with the short description of the model and algorithm.

Every analysis use an alignment made of 69 influenza A haemagglutinin DNA sequences.

# JC69-ML.json
- Model:
  - Fixed unrooted topology
  - JC69 substitution model
- Algorithm: maximum likelihood

# JC69-NNI-ML.json
- Model:
  - Unrooted topology optimized using NNIs
  - JC69 substitution model
- Algorithm: maximum likelihood

# GTR-G4-ML.json
- Model:
  - Fixed unrooted topology
  - GTR substitution model
  - Rate heterogeneity across site using discretized gamma distribution (4 categories)
- Algorithm: maximum likelihood

# JC69-time-ELBO.json
- Model:
  - Time tree
  - Fixed rooted topology
  - JC69 substitution model
  - Strict clock
  - Constant population size coalescent
- Algorithm:
  - variational inference
  - ELBO is maximized
  - 1 sample for calculating the gradient and 100 samples for the ELBO