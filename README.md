# Introduction

## Overview

This repository accompany the article by C. Poussot-Vassal, I-V. Gosea, P. Vuillemin and A.C. Antoulas in "Tensor-based multivariate function approximation: methods benchmarking and comparison", currently under reviewing, and its extendded [arXiv](https://arxiv.org/abs/2506.04791) version (regularly updated). The functions and script in this repository allow evaluating different codes for $n$-dimensional tensor approximation. More specifically, the current version evaluates the following codes:
- ["mlf1" and "mlf2"](https://github.com/cpoussot/mLF), implementing the multivariate Loewner Framework (Alg. 1 & 2), in Matlab 
- ["mdspack"](https://mordigitalsystems.fr/), implementing the  multivariate Loewner Framework, in Fortran (developped by MOR Digital Systems)
- ["kan1"](https://github.com/andrewpolar), implementing a Kolmogorov Arnold Network, in Matlab
- ["paaa" and "paaaalr"](https://github.com/lbalicki/parametric-AAA), implementing the parametric AAA and its low rank version, in Matlab
- ["tensorflow"](https://www.tensorflow.org/?hl=fr), implementing the Multi Layer Perceptron, in Python (not supported yet in this package)

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

The package contains the necessary material to constuct an evaluation report as the one given in ["Tensor-based multivariate function approximation: methods benchmarking and comparison"](https://arxiv.org/abs/2506.04791).

## Dependencies

- MATLAB R2023b or later (tested on this version)
- It is strongly recommended to download at least ["mlf"](https://github.com/cpoussot/mLF) since some printing and evaluation functions are embeded there.

## Detailed description

The package contains:
- `+run`: a set of functions for evaluation.
- `start_init`: a script that initialize the main variables,
    - it declares the main variables:
        - `SPACE_CAS=1:50`, the considered examples list (see `mlf.examples`)
        - `NTEST=500`, the number of random draw for evaluation of the model mismatch
        - `RESULT_PATH`, the path where you want to save the results (in the actual version, the path are my local ones, you should update them with yours)
        - `TEX_PATH`, the path where you want to save the LaTeX code generated and used for the report (in the actual version, the path are my local ones, you should update them with yours)
    - it adds to the Matlab path of the third parties software (in the actual version, the path are my local ones, you should update them with yours);
- `start_compare_step1`: this script applies the different approximation methods using different parametrizations.
- `start_compare_step2`: this script evaluates the approximation quality of the different methods with different tuning parameters w.r.t. the true function and keeps the best candidate per method (over all parametrizations).
- `start_compare_step3`: this script evaluates the best candidate for each method and report some statistics and figures.


## A simple procedure

- Set the variables in `start_init`
    - Set e.g. `SPACE_CAS=1:2`, `NTEST=500`
    - Chose the `RESULT_PATH` and `TEX_PATH` where you want to save your reults
    - Set the path for each method
- Run `start_compare_step1`
    - In `RESULT_PATH`, folders for the different methods are created, and in each folder, Matlab files with the computed models is created (e.g. `RESULT_PATH/mlf1/cas_1_mlf1.mat`).
    - Line 14: `METHOD_LIST` gathers all the method to be tested, but of course, to can restrict to one or few, or even add your own.
    - Line 17 to 47: you may try different parametrizations and change the possible parameters.
- Run `start_compare_step2`
    - In `RESULT_PATH`, in the same folder, a Matlab file containnig the best model is saved (e.g. `RESULT_PATH/mlf1/cas_1_mlf1_best.mat`).
- Run `start_compare_step3`
    - In `TEX_PATH/figures`, a folder with the results related to each case are saved  (e.g. `TEX_PATH/figures/case_1/all_stat.pdf`, `TEX_PATH/figures/case_1/eval_scaled.pdf`, `TEX_PATH/figures/case_1/table_main.tex`, `TEX_PATH/figures/case_1/text_loe.tex`, `TEX_PATH/figures/case_1/text_main.tex` and `TEX_PATH/figures/case_1/text_slide.tex`). These figures and LaTeX code will be used for reporting.
- Copy the files contained in `tex_pdf` in `TEX_PATH`
- Open `main.tex` and check line 31: originally it is `\def\CAS{2}`, where the `2` is related to the number of computed examples chosen in `SPACE_CAS`.
- Then compile the LaTeX file and obtain a PDf with the report.

## Feedbacks

Please send any comment to C. Poussot-Vassal (charles.poussot-vassal@onera.fr) if you want to report any bug or user experience issues.

## Disclaimer

This deposit consitutes a research code that accompany the paper mentionned above. It is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

Notice also that pathological cases may appear. A more advanced and professional code, to deal with practical and theoretical issues/limitations is currently under development by the authors.
