# fastEmu 

<!-- badges: start -->
[![R-CMD-check](https://github.com/statdivlab/fastEmu/workflows/R-CMD-check/badge.svg)](https://github.com/statdivlab/fastEmu/actions)
[![codecov](https://codecov.io/github/statdivlab/fastEmu/coverage.svg?branch=main)](https://app.codecov.io/github/statdivlab/fastEmu)
<!-- badges: end -->

`fastEmu` is an `R` package for estimating changes in the abundance of microbial categories (taxa, genes, etc.) associated with covariates using high throughput sequencing data.
`fastEmu` is a wrapper for the package [`radEmu`](https://github.com/statdivlab/radEmu) that implements a fast version of the hypothesis testing procedure implemented in `radEmu`. 

Check out [`radEmu`](https://github.com/statdivlab/radEmu) for a full list of reasons to use this method for your differential abundance analysis. Some highlights include:

- `radEmu` uses high throughput sequencing data to estimate changes in the "absolute abundance" of categories
- `radEmu` doesn't require data transformations, a reference taxon, or pseudocounts 
- `radEmu` is robust to differential sampling depths and differential detection of categories (taxa, genes, etc.)
- `radEmu` has great error control, including in small samples and data generated from different distributions

However, hypothesis testing in `radEmu` can be slow, especially when there is a large number of categories. This is why we created `fastEmu`! `fastEmu` uses a simplified model to test the same parameters\* tested in `radEmu`, but much faster. 

### What parameters are we testing exactly?

\*We don't estimate exactly the same parameters as in `radEmu`. In `radEmu`, we estimate 
the log fold difference in abundance associated with the covariate for each category *relative* to the typical log fold difference in abundance associated with the covariate across all categories. This lets us pick out the categories that are the most differentially abundant compared to the overall set of categories. In `fastEmu`, we cannot compare to the typical log fold difference all categories, because we need to reduce the size of our model in order to run it quickly. Instead, we compare to the typical log fold difference in a reference set of categories. The smaller the set, the faster the test will run.

You can either choose the reference set yourself (if you have a set of categories whose abundance you think varies little across covariate levels) or `fastEmu` can choose a data-driven reference set. This data-driven reference set is made up of the categories with the smallest estimated expected log fold difference across covariate levels relative to the typical log fold difference across the entire dataset. This means that the approximate median log fold difference within the reference set is very similar to the approximate median log fold difference within the full set of categories.

## Installation

To download `fastEmu`, use the code below.

``` r
# install.packages("devtools")
devtools::install_github("statdivlab/fastEmu")
library(fastEmu)
```

## Use

The vignettes demonstrate example usage of the main functions. Please [file an issue](https://github.com/statdivlab/fastEmu/issues) if you have a request for a tutorial that is not currently included. 

The following code will run `fastEmu` to estimate log fold differences in abundance for all categories in `Y` in comparison to the typical log fold difference in a data-driven reference set of categories (this parameter is nearly the same as the one estimated in `radEmu`). This code will then quickly run a robust score test for the hypothesis that the $5\text{th}$ category has a log fold difference in abundance associated with the treatment that is different than the typical log fold difference associated with the treatment across categories in the reference set. 

``` r
emu_test <- fastEmuFit(formula = ~ treatment,
                       data = sample_data, 
                       Y = count_data,
                       test_kj = data.frame(k = 2, j = 5),
                       reference_set = "data_driven")
```

## Documentation 

We additionally have a `pkgdown` [website](https://statdivlab.github.io/fastEmu/) that contains pre-built versions of our function [documentation](https://statdivlab.github.io/fastEmu/reference/index.html) and a [vignette](https://statdivlab.github.io/fastEmu/articles/intro_to_fastEmu.html).

## Citation

If you use `fastEmu` for your analysis, check back in soon for a preprint! 

## Bug Reports / Change Requests

If you encounter a bug or would like make a change request, please file it as an issue [here](https://github.com/statdivlab/fastEmu/intro_to_fastEmu).

If you're a developer, we would love to review your pull requests. 
