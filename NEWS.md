# fastEmu 2.0.0

* Major update, changes the way that the data-driven reference sets are determined, lets there be a reference set for each covariate instead of one reference set per analysis. This typically means the results will align closer with `radEmu`, although computation will take longer for more complex models with more covariates.
