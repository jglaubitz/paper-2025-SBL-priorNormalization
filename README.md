# Efficient sampling for sparse Bayesian learning using hierarchical prior normalization

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)

This repository contains code to reproduce the results presented in the
article
```bibtex
@article{glaubitz2025efficient,
  title={Efficient sampling for sparse Bayesian learning using hierarchical prior normalization},
  author={Glaubitz, Jan and Marzouk, Youssef},
  journal={Journal of Computational Physics},
  year={2025},
  month={05},
}
```

If you find these results useful, please cite the article mentioned above. 

## Abstract

We introduce an approach for efficient Markov chain Monte Carlo (MCMC) sampling for challenging high-dimensional distributions in sparse Bayesian learning (SBL). 
The core innovation involves using hierarchical prior-normalizing transport maps (TMs), which are deterministic couplings that transform the sparsity-promoting SBL prior into a standard normal one. 
We analytically derive these prior-normalizing TMs by leveraging the product-like form of SBL priors and Knothe--Rosenblatt (KR) rearrangements.
These transform the complex target posterior into a simpler reference distribution equipped with a standard normal prior that can be sampled more efficiently. 
Specifically, one can leverage the standard normal prior by using more efficient, structure-exploiting samplers.  
Our numerical experiments on various inverse problems---including signal deblurring, inverting the non-linear inviscid Burgers equation, and recovering an impulse image---demonstrate significant performance improvements for standard MCMC techniques.


## Numerical experiments

This repository contains all source code required to reproduce the numerical
experiments presented in the paper. 
It is developed for [Julia](https://julialang.org/) v1.10.0 in Jupyter Notebook.

To reproduce the numerical experiments, clone this repository and start Julia with the project set to the local directory:
```shell
git clone https://github.com/jglaubitz/paper-2025-SBL-priorNormalization.git
cd paper-2025-SBL-priorNormalization
julia --project=.
```
Since some of the experiments involve generating up to six MCMC chains, running all experiments can take some time on a single core. 
We thus recommend starting Julia with multiple threads by appending `--threads N` to the Julia
command, with `N` being the number of threads.


## Authors

- [Jan Glaubitz](https://www.janglaubitz.com) (Linköping University, Sweden)
- [Youssef Marzouk](https://uqgroup.mit.edu/people) (MIT, USA)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
