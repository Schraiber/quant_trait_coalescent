# quant_trait_coalescent

This R package simulates quantitative traits under simple scenarios, allowing the coalescent process, the mutational kernel that determines the effect sizes from mutation events, and mutational targets (loci).

Additionally, the package provides tools to characterize the resulting quantitative trait distributions, with particular interest paid to robust moment estimation (variance, skewness, and kurtosis estimators), non-normality tests, and non-unimodality tests.

# Use

## Mutational kernels

Upon mutation, the process draws quantitative trait effect sizes from a provided mutational kernel. Supported kernels:

```rnorm``` is the symmetric normal distribution

```rsn``` is the skewed normal distribution

```rlaplace``` is the Laplace (the symmetric exponential) exponential distribution

```rstable`` is the alpha-stable distribution

## Simulation

The function ```sim_qts``` simulates quantitative traits. This function relies on having ```ms``` installed (http://home.uchicago.edu/rhudson1/source/mksamples.html). Descriptions and default parameter values for ```sim_qts``` are given below. The function returns a list containing the simulation results.

```
# default values
num_replicates  = 50
num_loci        = 1:8                # Simulate for all values of log2(num_loci)
num_samples     = 1:8                #          for all values of log2(num_samples)
num_individuals = 2000               #          sampled from a population of size num_individuals
theta           = 0.1                #          for the eff. mutational rate (see paper)
kernel          = rnorm              #          drawing mutational effects from kernel
scale           = 1                  #          satisfying an eff. process scale (see paper)
ncores          = 16                 #          using mclapply with ncores processors
ms_exec         = "~/apps/msdir/ms"  #          calling "ms" using ms_exec
results = sim_qts(nrep=num_replicates,
                  loc=num_loci,
                  sam=num_samples,
                  npop=num_individuals,
                  kernel=kernel,
                  scale=scale,
                  ms_exec=ms_exec)
```

## Moments

The ```get_moment``` function computes the specified sample moment from ```sim_qts``` results.

```
m           = 4     # 2:variance, 3:skewness, 4:kurtosis
use_moments = TRUE  # TRUE:  use the moment function from the moments package
                    # FALSE: use the default R moment functions (var, skewness, kurtosis)
central     = TRUE  # TRUE:  use central moments (only works if use_moments enabled)
                    # FALSE: use moments about zero 

moment_results = get_moment(results,
                            loc=num_loci,
                            sam=num_samples,
                            moment=m,
                            use_moments=use_moments,
                            central=central)
```

The ```get_expected_moment``` function compues the expected moment given the process parameters

```
moment_expectation= get_expected_moment(kernel=kernel,
                                        loc=num_loci,
                                        sam=num_samples,
                                        theta=theta,
                                        sigma=1,
                                        skew=0,
                                        moment=4)
```


# Citation

For more details about the underlying theory, please see

Schraiber, J. G. and Landis M. L. Sensitivity of quantitative traits to mutational effects and number of loci.
(link when available)

