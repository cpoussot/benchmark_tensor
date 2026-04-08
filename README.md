# Introduction

## Overview

This repository accompany the article by C. Poussot-Vassal, I-V. Gosea, P. Vuillemin and A.C. Antoulas in ["Tensor-based multivariate function approximation: methods benchmarking and comparison"](https://arxiv.org/abs/2506.04791) currently under reviewing, and its extendded [arXiv](https://arxiv.org/abs/2506.04791) version (regularly updated). The functions and script in this repository allow evaluating different codes for $n$-dimensional tensor approximation. More specifically, the current version evaluates the following codes:
- ["mlf1" and "mlf2"](https://github.com/cpoussot/mLF), implementing the multivariate Loewner Framework (Alg. 1 & 2) in Matlab 
- ["mdspack"](https://mordigitalsystems.fr/), implementing the  multivariate Loewner Framework in Fortran (developped by MOR Digital Systems)
- ["kan1"](https://github.com/andrewpolar), implementing a Kolmogorov Arnold Network in Matlab
- ["paaa" and "paaaalr"](https://github.com/lbalicki/parametric-AAA), implementing the parametric AAA and its low rank version in Matlab
- ["tensorflow"](https://www.tensorflow.org/?hl=fr), implementing the Multi Layer Perceptron in Python (not supported in this pacakge)

## Contributions claim

- To suggest a comprehensive benchmark collection together with a methodology for tensor approximation with a surrogate model and, 
- To provide a plug-n-play manner to report the results.

## Main reference

```
@article{PVGVA:2026,
	Author	= {C. Poussot-Vassal and I-V. Gosea and P. Vuillemin A.C. Antoulas},
	Title 	= {Tensor-based multivariate function approximation: methods benchmarking and comparison},
	Doi 	= {https://doi.org/10.48550/arXiv.2506.04791},
	Journal = {arXiv},
	Volume 	= {},
	Number 	= {},
	Pages 	= {},
	Year 	= {2026},
	URL     = {https://arxiv.org/abs/2506.04791}, 
}
```


# The "benchmark_tensor" MATLAB package 

The package contains
- `+run`: a set of functions


## Dependencies

- MATLAB R2023b or later (tested on this version)
- Toolboxes: the one listed in the introduction/overview section
- It is strongly recommended to dpwnload at least ["mlf"](https://github.com/cpoussot/mLF) since some printing function are embdedde there.

## A simple MATLAB code example

We provide two examples to test the approach. Refer to [SIAM Review paper](https://doi.org/10.1137/24M1656657) / [arXiv paper](https://arxiv.org/abs/2405.00495) for notations and related equations. 

First add the path where the `+mlf` package is.

```Matlab
addpath("location_of_mlf") % Add the location of the +mlf package
```

## Feedbacks

Please send any comment to C. Poussot-Vassal (charles.poussot-vassal@onera.fr) if you want to report any bug or user experience issues.

## Disclaimer

This deposit consitutes a research code that accompany the paper mentionned above. It is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

Notice also that pathological cases may appear. A more advanced and professional code, to deal with practical and theoretical issues/limitations is currently under development by the authors.
