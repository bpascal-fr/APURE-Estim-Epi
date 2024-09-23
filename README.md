# Risk Estimate under a Nonstationary Generalized Autoregressive Model 



**Contributors:** B. Pascal (1) and S. Vaiter (2).  
(1) Nantes Université, École Centrale Nantes, CNRS, LS2N, UMR 6004, F-44000 Nantes, France.  
(2) CNRS, Université Côte d’Azur, LJAD, Nice,  France

Fundings: ANR-23-CE48-0009 [OptiMoCSI](https://optimocsi.cnrs.fr/) and ANR-18-CE40-0005 [GraVa](https://samuelvaiter.com/grava/).



This project contains the `Matlab` codes associated to the preprint:

> Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a Nonstationary Autoregressive Model for Data-Driven Reproduction Number  Estimation. Preprint. [arXiv:]()

### Data-driven piecewise linear estimation under nonstationary autoregressive Poisson model

Two demonstration scripts are provided:
- [`Demo_Synthetic_Poisson`](https://github.com/bpascal-fr/APURE-Estim-Epi/blob/main/Demo_Synthetic_Poisson.m), simple to run, all options set to the default values;
- [`Advanced_Synthetic_Example`](https://github.com/bpascal-fr/APURE-Estim-Epi/blob/main/Advanced_Synthetic_Poisson.m), all possible tunings unrolled for complete customization, enabling for example to reduce the computational load or to increase the precision of the explored parameter grid.

These scripts are accompanied with the ground truth used to generate synthetic observations, displayed in Figure 1 in [Pascal & Vaiter, 2024](), stored in [data](https://github.com/bpascal-fr/APURE-Estim-Epi/tree/main/data), so that the user can reproduce the numerical experiments presented in Section 4 of the associated preprint.

### Application to Data-Driven COVID-19 Reproduction Number Estimation

The demonstration script [`Demo_Epi_Covid`](https://github.com/bpascal-fr/APURE-Estim-Epi/blob/main/Demo_Epi_Covid.m) enables to download COVID-19 new infection counts reported in 200+ countries during 3 years of pandemic collected from National Health Services and made publicly available on the Johns Hopkins University [repository](https://coronavirus.jhu.edu/) and to estimate automatically the time-varying *weekly* instantaneous reproduction number under the extended epidemiological model proposed in Section 5 of [Pascal & Vaiter, 2024]().

## Project description

In [Pascal & Vaiter, 2024](), a novel *generalized nonstationary autoregressive* model has been proposed [5], encompassing, but not reducing to, one of the most popular model for the propagation of viral epidemics [1], a topic recently brought to the fore due to the COVID-19 pandemic [2,3].
This model is driven by a nonnegative time-varying *reproduction coefficient* $\mathsf{X}_t$ so that:
- when $\mathsf{X}_t>1$ the observed time series is exponentially increasing with time;
- while if $\mathsf{X}_t < 1$ the time series is shrinking exponentially fast.

The purpose is to accurately estimate $\{\mathsf{X}_t,  t = 1, ..., T\}$ from observations  
> $\mathsf{Y}_t = \mathcal{B} (\mathsf{X}_t \Psi_t(\mathsf{Y}))$  
>
corrupted by measurement noise $\mathcal{B}$ and involving a memory term $\Psi_t(\mathsf{Y})$, while using as little expert knowledge as possible. In particular, no ground truth of any sort is assumed available.

In epidemiology, observations consists in new infection counts, the measurement noise follows a Poisson distribution and the time-varying parameter to be estimated is the *reproduction number* $\mathsf{R}_t$, quantifying the intensity of the virus spread.
The major challenge to estimate $\{ \mathsf{R}_t, t=1,...,T\}$ is the low quality of infection counts reported during an ongoing epidemic in a crisis context, corrupted by significant administrative noise consisting of missing counts on week-ends and days off, cumulated counts and reporting errors.

To obtain accurate estimates of the time-varying reproduction number a variational framework have been introduced in [2,3], enforcing temporal smoothness of $\mathsf{R}_t$, which is expected for an indicator describing the propagation of a pathogen in a large susceptible population.
In [Pascal & Vaiter, 2024](), the administrative noise is smoothed out by considering *weekly* aggregated counts inspiring from [4], and a *scaled* Poisson model is introduced to better account for intrinsic fluctuations of infection counts, extending the state-of-the-art model proposed in [1].

In practice, the main bottleneck to the use of the use of variational estimators is that the accuracy of the estimate strongly depends on *hyperparameters* tuning.
Without available ground truth, hyperparameters are selected by minimizing specifically designed data-driven oracles, used as proxy for the estimation error.
Focusing on the nonstationary autoregressive *Poisson model*, [Pascal & Vaiter, 2024]() generalized the Stein's Unbiased Risk Estimate formalism to construct asymptotically unbiased risk estimators based on the derivation of an original autoregressive counterpart of Stein's lemma, referred to as prediction and estimation Autoregressive Poisson Unbiased Risk Estimates (APURE).

The `APURE-Estim-Epi`toolbox implements the evaluation of robustified prediction and estimation APURE, leveraging Finite Difference and Monte Carlo strategies and the minimization of these data-driven oracles selecting optimal hyperparameters.
Accurate piecewise linear estimates of the reproduction coefficient are provided.
The overall procedure is exemplified both on synthetic observations and on real COVID-19 infection counts.


[1] Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013). A new framework and software to estimate time-varying reproduction numbers during epidemics. *American Journal of Epidemiology*, 178(9), 1505-1512.

[2] Abry, P., Pustelnik, N., Roux, S., Jensen, P., Flandrin, P., Gribonval, R., Lucas, C.-G., Guichard, É., Borgnat, P., & Garnier, N. (2020). Spatial and temporal regularization to estimate COVID-19 reproduction number R(t): Promoting piecewise smoothness via convex optimization. *PlosOne*, 15(8), e0237901.

[3] Pascal, B., Abry, P., Pustelnik, N., Roux, S., Gribonval, R., & Flandrin, P. (2022). Nonsmooth convex optimization to estimate the Covid-19 reproduction number space-time evolution with robustness against low quality data. *IEEE Transactions on Signal Processing*, 70, 2859–2868.

[4] Nash, R. K., Bhatt, S., Cori, A., & Nouvellet, P. (2023). Estimating the epidemic reproduction number from temporally aggregated incidence data: A statistical modelling approach and software tool. *PLOS Computational Biology*, 19(8), e1011439.

[5] Pascal, B., Vaiter, S. (2024, September). Risk Estimate under a Nonstationary Autoregressive Model for Data-Driven Reproduction Number  Estimation. *Preprint*. [arXiv:]()

## Installation and dependencies

To download and install the toolbox:  

> - open a Terminal in your working folder,
> - execute `git clone https://github.com/bpascal-fr/APURE-Estim-Epi.git`.

The Matlab toolbox [stein-piecewise-filtering](https://github.com/bpascal-fr/stein-piecewise-filtering) should be downloaded and placed in the folder containing the toolbox:

> - `cd APURE-Estim-Epi`,
> - `git clone https://github.com/bpascal-fr/stein-piecewise-filtering.git`.
>