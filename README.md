
<!-- README.md is generated from README.Rmd. Please edit that file -->

# epiwave.mapping

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

Prototype R package and example code for computationall-yefficient,
fully-Bayesian semi-mechanistic spatio-temporal mapping of disease
infection incidence simultaneously from spatially-aggregated (e.g. at
health facility level) clinical case count timeseries, and targeted
infection prevalence surveys at specific point locations. This model and
software is based on ideas in the epiwave R package and uses some
functions implemented there, but it currently relies on the user to
develop their own model using the greta R package, rather than providing
wrapper functions to set up the model for the user.

## Installation

You can install the package from github, using the `remotes` package:

``` r
remotes::install_github("idem-lab/epiwave.mapping",
                        dependencies = TRUE)
```

## Model

Our main aim is to model, and generate spatio-temporal predictive maps
of, variation in *infection incidence*: the average number of new
infections per time period, per member of the population. We model the
logarithm of infection incidence with a space-time Gaussian process.
Since infection incidence is rarely perfectly observed, we infer it from
two independent and complementary datastreams:

- **Clinical case timeseries**: Counts of the numbers of cases reported
  over time, but potentially spatially aggregated across e.g. the
  catchment areas of health facilities, or the entire region, rather
  than at spatially precise coordinates

- **Infection prevalence survey points**: Random samples of individuals
  at a specific point location, with the number tested and number
  positive recorded. Typically these data are available only very
  infrequently.

Since neither datastream provides a direct estimate of infection
incidence, we model the generative observation process yielding each
datastream, using informative priors for key parameters wherever
possible. A number of options are available for
computationally-efficient space-time Gaussian process (GP) modelling for
applications of this type. Below we implement an approach well-suited to
our health facility catchment observation model and fully Bayesian
inference using Hamiltonian Monte Carlo: a separable space-time GP
employing a block-circulant embedding for the spatial process. This is
similar to the approach used in the `lgcp` R package, and detailed in
Diggle et al. (2013). We note that the clinical incidence observation
model we employ is a particular case of a Log-Gaussian Cox process. Our
contribution is accounting for delays in reporting, and linking model
this to prevalence survey data, in a practical and easily-extensible
framework.

### Clinical incidence model

We model the observed number of clinical cases
![C\_{h,t}](https://latex.codecogs.com/png.latex?C_%7Bh%2Ct%7D "C_{h,t}")
of our disease of interest in health facility
![h](https://latex.codecogs.com/png.latex?h "h") during discrete time
period ![t](https://latex.codecogs.com/png.latex?t "t") as a Poisson
random variable:

![C\_{h,t} \sim \text{Poisson}(\hat{C}\_{h,t})](https://latex.codecogs.com/png.latex?C_%7Bh%2Ct%7D%20%5Csim%20%5Ctext%7BPoisson%7D%28%5Chat%7BC%7D_%7Bh%2Ct%7D%29 "C_{h,t} \sim \text{Poisson}(\hat{C}_{h,t})")

The expectation of this Poisson random variable (the modelled/expected
number of clinical cases) is given by a weighted sum of (unobserved but
modelled) expected pixel-level clinical case counts
![\hat{C}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7BC%7D_%7Bl%2Ct%7D "\hat{C}_{l,t}")
at each of ![L](https://latex.codecogs.com/png.latex?L "L") pixel
locations ![l](https://latex.codecogs.com/png.latex?l "l"):

![\hat{C}\_{h,t} = \sum\_{l=1}^{L}{\hat{C}\_{l,t}w\_{l,h}}](https://latex.codecogs.com/png.latex?%5Chat%7BC%7D_%7Bh%2Ct%7D%20%3D%20%5Csum_%7Bl%3D1%7D%5E%7BL%7D%7B%5Chat%7BC%7D_%7Bl%2Ct%7Dw_%7Bl%2Ch%7D%7D "\hat{C}_{h,t} = \sum_{l=1}^{L}{\hat{C}_{l,t}w_{l,h}}")

where weights
![w\_{l,h}](https://latex.codecogs.com/png.latex?w_%7Bl%2Ch%7D "w_{l,h}")
give the ‘membership’ of the population in each pixel location to each
health facility, where each case at a given location
![l](https://latex.codecogs.com/png.latex?l "l") has probability
![w\_{l,h}](https://latex.codecogs.com/png.latex?w_%7Bl%2Ch%7D "w_{l,h}")
of reporting at health facility
![h](https://latex.codecogs.com/png.latex?h "h"), and therefore
![\sum\_{l=1}^L w\_{l,h} = 1](https://latex.codecogs.com/png.latex?%5Csum_%7Bl%3D1%7D%5EL%20w_%7Bl%2Ch%7D%20%3D%201 "\sum_{l=1}^L w_{l,h} = 1").
In practice this could be either a proportional (fractions of the
population attend different health facilities) or a discrete (the
population of each pixel location attends only one nearby health
facility) mapping.

At each pixel location, we model the unobserved clinical case count
![\hat{C}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7BC%7D_%7Bl%2Ct%7D "\hat{C}_{l,t}")
from the modelled number of new infections
![\hat{I}\_{l,t'}](https://latex.codecogs.com/png.latex?%5Chat%7BI%7D_%7Bl%2Ct%27%7D "\hat{I}_{l,t'}")
at the same location, during previous time periods
![t'](https://latex.codecogs.com/png.latex?t%27 "t'"), from the fraction
of infections
![\gamma\_{l,t'}](https://latex.codecogs.com/png.latex?%5Cgamma_%7Bl%2Ct%27%7D "\gamma_{l,t'}")
at that location and previous time period that would result in a
recorded clinical case, and a probability distribution
![\pi(t-t')](https://latex.codecogs.com/png.latex?%5Cpi%28t-t%27%29 "\pi(t-t')")
over delays ![t-t'](https://latex.codecogs.com/png.latex?t-t%27 "t-t'")
from infection to diagnosis and reporting:

![\hat{C}\_{l,t} = \sum\_{t-t' = 0}^{\tau\_\pi}{\hat{I}\_{l,t'} \\ \gamma\_{l,t'} \\ \pi(t-t')}](https://latex.codecogs.com/png.latex?%5Chat%7BC%7D_%7Bl%2Ct%7D%20%3D%20%5Csum_%7Bt-t%27%20%3D%200%7D%5E%7B%5Ctau_%5Cpi%7D%7B%5Chat%7BI%7D_%7Bl%2Ct%27%7D%20%5C%2C%20%5Cgamma_%7Bl%2Ct%27%7D%20%5C%2C%20%5Cpi%28t-t%27%29%7D "\hat{C}_{l,t} = \sum_{t-t' = 0}^{\tau_\pi}{\hat{I}_{l,t'} \, \gamma_{l,t'} \, \pi(t-t')}")

The distribution
![\pi(\Delta_t)](https://latex.codecogs.com/png.latex?%5Cpi%28%5CDelta_t%29 "\pi(\Delta_t)")
gives the discrete and finite (support on
![\Delta_t \in (0, \tau\_\pi)](https://latex.codecogs.com/png.latex?%5CDelta_t%20%5Cin%20%280%2C%20%5Ctau_%5Cpi%29 "\Delta_t \in (0, \tau_\pi)"))
probability distribution over delays from infection to reporting,
indexed on the time-periods considered in the model. This temporal
reweighting to account for a distribution over possible delays can be
considered as a ‘discrete-time convolution’ with
![\pi(\Delta_t)](https://latex.codecogs.com/png.latex?%5Cpi%28%5CDelta_t%29 "\pi(\Delta_t)")
the ‘kernel’. Below we discuss efficient methods for computing these
convolutions in greta.

### Infection prevalence model

We model the observed number of individuals who test positive for
infection
![N^+\_{l,t}](https://latex.codecogs.com/png.latex?N%5E%2B_%7Bl%2Ct%7D "N^+_{l,t}")
in an infection prevalence survey at location
![l](https://latex.codecogs.com/png.latex?l "l") at time
![t](https://latex.codecogs.com/png.latex?t "t") as a binomial sample,
given the number of individuals tested
![N\_{l,t}](https://latex.codecogs.com/png.latex?N_%7Bl%2Ct%7D "N_{l,t}"),
and the modelled prevalence of infections
![\hat{p}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7Bp%7D_%7Bl%2Ct%7D "\hat{p}_{l,t}")
across the population at that time/place:

![N^+\_{l,t} \sim \text{Binomial}(N\_{l,t}, \hat{p}\_{l,t})](https://latex.codecogs.com/png.latex?N%5E%2B_%7Bl%2Ct%7D%20%5Csim%20%5Ctext%7BBinomial%7D%28N_%7Bl%2Ct%7D%2C%20%5Chat%7Bp%7D_%7Bl%2Ct%7D%29 "N^+_{l,t} \sim \text{Binomial}(N_{l,t}, \hat{p}_{l,t})")

Similarly to clinical incidence, we model the infection prevalence at a
given location and time as a discrete-time convolution over previous
infection counts, divided by the total population in the pixel
![M_l](https://latex.codecogs.com/png.latex?M_l "M_l"), and the duration
of the time-period in days
![d](https://latex.codecogs.com/png.latex?d "d"):

![\hat{p}\_{l,t} = \frac{1}{M_l}\sum\_{t-t' = 0}^{\tau_q}{q(t-t') \\ \hat{I}\_{l,t'} \\ d^{-1}}](https://latex.codecogs.com/png.latex?%5Chat%7Bp%7D_%7Bl%2Ct%7D%20%3D%20%5Cfrac%7B1%7D%7BM_l%7D%5Csum_%7Bt-t%27%20%3D%200%7D%5E%7B%5Ctau_q%7D%7Bq%28t-t%27%29%20%5C%2C%20%5Chat%7BI%7D_%7Bl%2Ct%27%7D%20%5C%2C%20d%5E%7B-1%7D%7D "\hat{p}_{l,t} = \frac{1}{M_l}\sum_{t-t' = 0}^{\tau_q}{q(t-t') \, \hat{I}_{l,t'} \, d^{-1}}")

In this case, the kernel
![q(t-t')](https://latex.codecogs.com/png.latex?q%28t-t%27%29 "q(t-t')")
is not a probability distribution, but the kernel of a convolution
operation that maps the daily average number of new infections in
previous time periods
![t'](https://latex.codecogs.com/png.latex?t%27 "t'") to the number of
individuals who would test positive (to the diagnostic used) on a
randomly-selected day in the timeperiod `t`. This kernel, which is
specific to the tim eperiod chosen for modelling the temporal process,
can be computed from a related function:
![q\_{\text{daily}}(\Delta_t)](https://latex.codecogs.com/png.latex?q_%7B%5Ctext%7Bdaily%7D%7D%28%5CDelta_t%29 "q_{\text{daily}}(\Delta_t)")
which gives the probability that an individual will test positive
![\Delta_t](https://latex.codecogs.com/png.latex?%5CDelta_t "\Delta_t")
days after infection. The integral of this function is the average
number of days an infected person would test positive for. The function
![q\_{\text{daily}}()](https://latex.codecogs.com/png.latex?q_%7B%5Ctext%7Bdaily%7D%7D%28%29 "q_{\text{daily}}()")
can be estimated from empirical data on how the test sensitivity and
duration of infection vary over time since infection. This package
provides functionality: `transform_convolution_kernel()` to construct
timeperiod-specific kernels such as
![q()](https://latex.codecogs.com/png.latex?q%28%29 "q()") from daily
kernels like
![q\_{\text{daily}}()](https://latex.codecogs.com/png.latex?q_%7B%5Ctext%7Bdaily%7D%7D%28%29 "q_{\text{daily}}()")
for an arbitrary modelling timeperiod.

By convolving of the number of new infections *per day*
![\hat{I}\_{l,t'} \\ d^{-1}](https://latex.codecogs.com/png.latex?%5Chat%7BI%7D_%7Bl%2Ct%27%7D%20%5C%2C%20d%5E%7B-1%7D "\hat{I}_{l,t'} \, d^{-1}")
in all previous time periods with the fraction of those that would test
positive in a survey in time period
![t](https://latex.codecogs.com/png.latex?t "t"), gives an estimate of
the *number* of positive-testing people in the population at location
![l](https://latex.codecogs.com/png.latex?l "l"), on any given day
within the time period ![t](https://latex.codecogs.com/png.latex?t "t").
Dividing this by the population of that location,
![M_l](https://latex.codecogs.com/png.latex?M_l "M_l"), therefore yields
an estimate of the population proportion testing positive on any given
day; the parameter of the binomial distribution that captures the
prevalence survey data.

Note that our definition of
![\hat{p}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7Bp%7D_%7Bl%2Ct%7D "\hat{p}_{l,t}")
is the population prevalence of infections *that would test positive
using that diagnostic method* rather than the true fraction
infected/infectious at any one time. We also assume here that tests have
perfect specificity, though the model can easily be adapted to
situations where that is not the case.

### Infection incidence model

The expected number of new infections
![\hat{I}\_{l,t'}](https://latex.codecogs.com/png.latex?%5Chat%7BI%7D_%7Bl%2Ct%27%7D "\hat{I}_{l,t'}")
in location ![l](https://latex.codecogs.com/png.latex?l "l") during time
period ![t](https://latex.codecogs.com/png.latex?t "t") is modelled as
the product of the population
![M_l](https://latex.codecogs.com/png.latex?M_l "M_l") at that location,
the length of the timeperiod in days
![d](https://latex.codecogs.com/png.latex?d "d"), and the daily
infection incidence
![f\_{l,t}](https://latex.codecogs.com/png.latex?f_%7Bl%2Ct%7D "f_{l,t}")
at that location and time:

![\hat{I}\_{l,t'} = d \\ M_l \\ \hat{f}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7BI%7D_%7Bl%2Ct%27%7D%20%3D%20d%20%5C%2C%20M_l%20%5C%2C%20%5Chat%7Bf%7D_%7Bl%2Ct%7D "\hat{I}_{l,t'} = d \, M_l \, \hat{f}_{l,t}")

Whilst we mechanistically model the observation processes yielding our
data types, we employ a geostatistical approach to modelling
spatio-temporal variation in infection incidence, with spatio-temporal
covariates
![\mathbf{X}\_{l,t}](https://latex.codecogs.com/png.latex?%5Cmathbf%7BX%7D_%7Bl%2Ct%7D "\mathbf{X}_{l,t}")
and a space-time random effect
![\epsilon\_{l,t}](https://latex.codecogs.com/png.latex?%5Cepsilon_%7Bl%2Ct%7D "\epsilon_{l,t}")
with zero-mean Gaussian process prior:

![\text{log}(f\_{l,t}) = \alpha +\mathbf{X}\_{l,t} \beta + \epsilon\_{l,t}](https://latex.codecogs.com/png.latex?%5Ctext%7Blog%7D%28f_%7Bl%2Ct%7D%29%20%3D%20%5Calpha%20%2B%5Cmathbf%7BX%7D_%7Bl%2Ct%7D%20%5Cbeta%20%2B%20%5Cepsilon_%7Bl%2Ct%7D "\text{log}(f_{l,t}) = \alpha +\mathbf{X}_{l,t} \beta + \epsilon_{l,t}")

![\epsilon \sim GP(0, \mathbf{K})](https://latex.codecogs.com/png.latex?%5Cepsilon%20%5Csim%20GP%280%2C%20%5Cmathbf%7BK%7D%29 "\epsilon \sim GP(0, \mathbf{K})")

where ![\alpha](https://latex.codecogs.com/png.latex?%5Calpha "\alpha")
is a scalar intercept term,
![\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\beta") is a
vector of regression coefficients against the environmental covariates,
and
![\mathbf{K}](https://latex.codecogs.com/png.latex?%5Cmathbf%7BK%7D "\mathbf{K}")
is the space-time covariance function of the Gaussian process over
random effects
![\epsilon\_{l,t}](https://latex.codecogs.com/png.latex?%5Cepsilon_%7Bl%2Ct%7D "\epsilon_{l,t}").

There are many choices of space-time covariance structure for
![\mathbf{K}](https://latex.codecogs.com/png.latex?%5Cmathbf%7BK%7D "\mathbf{K}"),
though we use a separable combination of an isotropic spatial covariance
function with a temporal covariance function, to enable the use of a
range of computationally efficient simulation and calculation methods:

![K\_{l,t,l',t'} = \sigma^2 \\ K\_{\text{space}}(\|\|l-l'\|\|; \phi) \\ K\_{\text{time}}(\|t-t'\|; \theta)](https://latex.codecogs.com/png.latex?K_%7Bl%2Ct%2Cl%27%2Ct%27%7D%20%3D%20%5Csigma%5E2%20%5C%2C%20K_%7B%5Ctext%7Bspace%7D%7D%28%7C%7Cl-l%27%7C%7C%3B%20%5Cphi%29%20%5C%2C%20K_%7B%5Ctext%7Btime%7D%7D%28%7Ct-t%27%7C%3B%20%5Ctheta%29 "K_{l,t,l',t'} = \sigma^2 \, K_{\text{space}}(||l-l'||; \phi) \, K_{\text{time}}(|t-t'|; \theta)")

where
![\sigma^2](https://latex.codecogs.com/png.latex?%5Csigma%5E2 "\sigma^2")
is the marginal variance (amplitude) of the resultant Gaussian process,
![K\_{\text{space}}(.; \phi)](https://latex.codecogs.com/png.latex?K_%7B%5Ctext%7Bspace%7D%7D%28.%3B%20%5Cphi%29 "K_{\text{space}}(.; \phi)")
is an (unit variance) isotropic spatial kernel on euclidean distances
![\|\|l-l'\|\|](https://latex.codecogs.com/png.latex?%7C%7Cl-l%27%7C%7C "||l-l'||")
with parameter
![\phi \> 0](https://latex.codecogs.com/png.latex?%5Cphi%20%3E%200 "\phi > 0")
controlling the range of spatial correlation, and
![K\_{\text{time}}(\|t-t'\|; \theta)](https://latex.codecogs.com/png.latex?K_%7B%5Ctext%7Btime%7D%7D%28%7Ct-t%27%7C%3B%20%5Ctheta%29 "K_{\text{time}}(|t-t'|; \theta)")
is a a temporal kernel on time differences
![\|t-t'\|](https://latex.codecogs.com/png.latex?%7Ct-t%27%7C "|t-t'|"),
with parameter
![\theta](https://latex.codecogs.com/png.latex?%5Ctheta "\theta")
controlling the range of temporal correlation. Again for computational
reasons, we prefer Markovian kernels for
![K\_{\text{time}}](https://latex.codecogs.com/png.latex?K_%7B%5Ctext%7Btime%7D%7D "K_{\text{time}}").

## Computational approaches

### Inference

We fit the model via fully Bayesian inference, using MCMC (Hamiltonian
Monte Carlo) in the greta R package. Whilst a number of alternative
inference algorithms have been proposed and widely adopted for Bayesian
and non-Bayesian estimation of spatio-temporal Gaussian process models
(e.g. INLA), in out case, the aggregation of clinical cases across
catchments (a feature required by the model structure) leads to a
non-factorising likelihood, which nullifies many of the computational
advantages of that method. Further, the mapping from the Gaussian
process to the expectations of the two sampling distributions is
non-linear (due to sums of exponents in various places), requiring
computationally costly linearisation of the non-linearities. By
employing a fully-Bayesian MCMC approach, and focussing on
computationally-efficient model and implementation choices, rather than
approximations, we are able to estimate the model parameters relatively
quickly, with full treatment, propagation, and quantification of model
uncertainty, and retaining the ability to easily modify and interrogate
the model.

### Discrete-time convolution

This discrete convolutions (weighted sums over previous time points)
used to compute
![\hat{C}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7BC%7D_%7Bl%2Ct%7D "\hat{C}_{l,t}")
and
![\hat{p}\_{l,t}](https://latex.codecogs.com/png.latex?%5Chat%7Bp%7D_%7Bl%2Ct%7D "\hat{p}_{l,t}")
can be computationally intensive, depending on the number of time
periods modelled. A number of efficient computational approaches exist
to overcome these, and the optimal approach to use in greta depends on
the size of
![\tau^{max}](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D "\tau^{max}"):
the maximum duration of the delay in terms of the number of timeperiods
being considered. Where
![\tau^{max}](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D "\tau^{max}")
is relatively large (say,
![\tau^{max} \> 3](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D%20%3E%203 "\tau^{max} > 3")),
implementation of the convolution as a matrix multiply is likely to be
the most efficient in greta. If
![\tau^{max} \> 3](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D%20%3E%203 "\tau^{max} > 3")
and the number of time periods being modelled is much larger than
![\tau^{max}](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D "\tau^{max}"),
implementing this as a sparse matrix multiply (skipping computation on
zero elements of the convolution matrix) is likely to be most efficient.
Where
![\tau^{max}](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D "\tau^{max}")
is small (say,
![\tau^{max} \leq 3](https://latex.codecogs.com/png.latex?%5Ctau%5E%7Bmax%7D%20%5Cleq%203 "\tau^{max} \leq 3"))
the convolution can instead be computed with a sum of dense vectorised
additions and subtractions. These convolution approaches are implemented
in this package, and demonstrated below.

### Gaussian process simulation

Naive implementation of the full Gaussian process (GP) model is
typically very computationally intensive, due to the fact that the most
expensive step (inverting a dense covariance matrix) scales cubically
with the number of evaluation points. In the model we propose (and as in
a Log-Gaussian Cox process), we need to evaluate the GP at every pixel
location and time period in the study frame, so that we can aggregate
the expected clinical case count to compute the likelihood. This results
in a computationally impractical algorithm for all but the smallest
study frames.

One solution would be to employ an approximation to the full GP,
evaluated at only a subset of locations and times in the study frame,
and approximating the clinical case count calculation with some smaller
finite sum. Candidate GP approximation approaches include the SPDE
approximation to a Matern-type spatial kernel, on a computational mesh
(as used in INLA, Lindgren *et al.*, 2011), a sparse GP method over a
limited set of inducing points (Quinonera-Candela *et al.*, ?; see
`greta.gp`), or one of the closely-related penalised spline methods (see
`greta.gam`).

An appealing alternative in this case is to compute the full spatial GP
for every pixel in a regular grid (the same we use to record spatial
covariate values), by exploiting the block-circulant structure in the
resulting spatial covariance matrix. This requires using projected
coordinates, and expanding the spatial study frame (by a factor of 2 in
each dimension) to map the projection onto a torus in such a way that
the distances between pairs of pixels are preserved. This approach
enables simulation of the GP across all pixels, for a given time period,
via the fast fourier transform (FFT) - which scales only linearly with
the number of locations considered (? check when have wifi). This
enables us to simulate values of the spatial GP across all pixels very
cheaply, with no need to approximate the GP. In separable combination
with a Markovian temporl kernel, this enables very rapid inference.

## Example application

We demonstrate the model with application to mapping the (simulated)
infection incidence of malaria in Kenya from simulated clinical
incidence and infection prevalence data. Using the `sim_data()` function
from this package, we use a real grid of environmental covariates, but
simulate the locations of health facilities, prevalence survey
locations, and all datasets to which we will then fit the model.

First we load the bioclim covariates using the `terra` and `geodata` R
package:

``` r
library(terra)
library(geodata)
# download at a half-degree resolution. See ?geodata_path for how to save the data between sessions
bioclim_kenya <- geodata::worldclim_country(
  country = "KEN",
  var = "bio",
  res = 0.5,
  path = file.path(tempdir(), "bioclim")
)

pop <- geodata::population(
  year = 2020,
  res = 0.5,
  path = file.path(tempdir(), "population")
)

# crop and mask to Kenya
pop_kenya <- mask(pop, bioclim_kenya)
```

    #> Warning: package 'terra' was built under R version 4.2.3
    #> terra 1.7.71

Now we lower the spatial resolution. a little (to make the example run
faster) and simulate some fake data, using the first 5 covariates

``` r
bioclim_kenya <- terra::aggregate(bioclim_kenya, 10)
pop_kenya <- terra::aggregate(pop_kenya, 10)

library(epiwave.mapping)
set.seed(1)
data <- sim_data(
  covariates_rast = bioclim_kenya[[1:5]],
  population_rast = pop_kenya, 
  years = 2020:2024,
  n_health_facilities = 100,
  n_prev_surveys = 30)
```
