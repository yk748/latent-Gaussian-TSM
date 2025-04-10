# Latent Gaussian models for multivariate and high-dimensional discrete-valued time series
Latent Gaussian models are designed for modeling count time series in multivariate, possibly high-dimensional settings. The key idea is to link the second-order properties of the discrete-valued observations to those of the latent Gaussian process through so-called link functions. This repository primarily reproduces the results in [Kim et al. (2024)](https://arxiv.org/pdf/2307.10454), where the Gaussian dynamic factor model is assumed as the latent process. It also includes code for related models, such as those in [Kim (2023)](https://www.proquest.com/openview/b7d1eae2131e518bc3af4ca4f2816513/1?pq-origsite=gscholar&cbl=18750&diss=y) and [D&uuml;ker et al. (2024)](https://arxiv.org/pdf/2301.00491), where the latent Gaussian processes are assumed to be vector autoregressive models. We are developing a fully functioning package, which will be released soon.
