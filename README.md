# Latent Gaussian models for multivariate and high-dimensional discrete-valued time series

Latent Gaussian models are designed for modeling count time series in multivariate, possibly high-dimensional settings. The key idea is to link the second-order properties of the discrete-valued observations to those of the latent Gaussian process through so-called link functions. This repository primarily reproduces the results in Kim, Düker, Fisher, and Pipiras ([2024](#ref-latent_dfm)), where the Gaussian dynamic factor model is assumed as the latent process. It also includes code for related models, such as those in Kim ([2023](#ref-latent_thesis)) and Düker, Lund, Pipiras ([2024](#ref-latent_var)), where the latent Gaussian processes are assumed to be vector autoregressive models. We are developing a fully functioning package, which will be released soon.

## References

<div id="ref-latent_var" class="references">
Marie-Christine Düker, Robert Lund, and Vladas Pipiras. 2024. "High-dimensional latent Gaussian count time series: Concentration results for autocovariances and applications." <em>Electronic Journal of Statistics</em>, 18(2), 5484–5562. [https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-18/issue-2/High-dimensional-latent-Gaussian-count-time-series--Concentration-results/10.1214/24-EJS3125.full](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-18/issue-2/High-dimensional-latent-Gaussian-count-time-series--Concentration-results/10.1214/24-EJS3125.full)

</div>
  
<div id="ref-latent_dfm" class="references">
Younghoon Kim, Marie-Christine Düker, Zachary F. Fisher, and Vladas Pipiras. 2023. "Latent Gaussian dynamic factor modeling and forecasting for multivariate count time series." <em>arXiv preprint arXiv:2307.10454.</em> <a href="https://arxiv.org/abs/2307.10454">https://arxiv.org/abs/2307.10454</a>

</div>

<div id="ref-latent_thesis" class="references">
Younghoon Kim. 2023. <em>Modeling Multiple-Subject and Discrete-Valued High-Dimensional Time Series.</em> PhD thesis, The University of North Carolina at Chapel Hill. <a href="https://www.proquest.com/openview/b7d1eae2131e518bc3af4ca4f2816513/1?pq-origsite=gscholar&cbl=18750&diss=y">https://www.proquest.com/openview/b7d1eae2131e518bc3af4ca4f2816513/1?pq-origsite=gscholar&cbl=18750&diss=y</a>

